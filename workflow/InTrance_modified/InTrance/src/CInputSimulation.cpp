/*!
\file CInputSimulation.cpp
\version 1.7
\brief CInputSimulation class implementation.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 08.11.2006 (AO)

\b Purpose:
   This file provides the interface of the class CInputSimulation,
   an auxiliary class that contains the InTrance simulation input. It parses
   those information from a given input file.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.

\b Copyright:
   (c) 2001-2006 Bernd Dachwald, Andreas Ohndorf
*/

// InTrance
#include <configInTrance.h>

// STL
#include <algorithm>
#include <sstream>

// InTrance
#include <CInputSimulation.h>

using namespace std;
using namespace SpaceLib;

extern TrajectoryType tt;         // global in CInTranceFitness.cpp
extern OptimizationGoalType goal; // global in CInTranceFitness.cpp

const std::string S_NOT_APPLICABLE = "n/a";


//------------------------------------------------------------------------------
void CInputSimulation::init(void)
{
  setSolarSystem(0);
  setInitialStateType(DFLT_INIT_STATE_TYPE);
  setInitCentralBodyName(DFLT_INIT_CENTRAL_BODY_NAME);
  setInitCentralBody(DFLT_INIT_CENTRAL_BODY);
  if (getInitCentralBody() != SUN) setHeliocentricSimulation(false);
  else setHeliocentricSimulation(true);

  setACTDistanceFactor(DFLT_ACT_DISTANCE_FACTOR);
  setACTExponent(DFLT_ACT_EXPONENT);

  min_start_time = DFLT_MIN_START_TIME;
  max_start_time = DFLT_MAX_START_TIME;

  // integration interval
  integration_interval = DFLT_ITGR_INTVL;
  ITGR_INTVL_REF_UNIT = SpaceLib::U_DAY;

  dynamic_integration_interval_reduction = DFLT_ITGR_INTVL_RED;
  max_steps = DFLT_MAX_CTRL_STEPS;

  min_output_points = DFLT_MIN_OUTPUT_POINTS;
  resolution_factor = DFLT_RESOLUTION_FACTOR;

  setInitialStateType(BODY_STATE);
  setInitialBody(SS_BARY);

  initInitialOrbitStateFields();  
  initInitialVectorStateFields();

  launch_v_inf.resize(N_VINF_PARAMS, 0.0);
  launch_v_inf_input_unit.resize(N_VINF_PARAMS, U_NO_UNIT);
  DFLT_V_INF_INPUT_UNIT.resize(N_VINF_PARAMS, U_NO_UNIT);
  DFLT_V_INF_INPUT_UNIT[VINF_MAX] = U_KM_PER_SEC;
  DFLT_V_INF_INPUT_UNIT[VINF_MIN] = U_KM_PER_SEC;
  DFLT_V_INF_INPUT_UNIT[VINF_AZI_MAX] = U_RADIAN;
  DFLT_V_INF_INPUT_UNIT[VINF_AZI_MIN] = U_RADIAN;
  DFLT_V_INF_INPUT_UNIT[VINF_ELE_MAX] = U_RADIAN;
  DFLT_V_INF_INPUT_UNIT[VINF_ELE_MIN] = U_RADIAN;

  INT_V_INF_UNIT.resize(N_VINF_PARAMS, U_NO_UNIT);
  INT_V_INF_UNIT[VINF_MAX] = U_AU_PER_DAY;
  INT_V_INF_UNIT[VINF_MIN] = U_AU_PER_DAY;
  INT_V_INF_UNIT[VINF_AZI_MAX] = U_RADIAN;
  INT_V_INF_UNIT[VINF_AZI_MIN] = U_RADIAN;
  INT_V_INF_UNIT[VINF_ELE_MAX] = U_RADIAN;
  INT_V_INF_UNIT[VINF_ELE_MIN] = U_RADIAN;
  
  setTrajectoryType(RENDEZVOUS);
  setOptGoalType(DFLT_OPT_GOAL_TYPE);
  setTargetStateType(BODY_STATE);
  setTargetBodyId(SS_BARY);

  number_of_ga_bodies = 0;
  ga_body.resize(0);

  initTargetOrbitalElementsFields();

  // maximum allowed relative target distance
  setMaxFinalRelDistance(DFLT_MAX_ABS_REL_DISTANCE);
  max_delta_r_input_unit = U_NO_UNIT;
  DFLT_MAX_DELTA_R_UNIT = U_KILOMETER;
  INT_MAX_DELTA_R_UNIT = U_AU;

  // maximum allowed relative velocity
  setMaxFinalRelVelocity(DFLT_MAX_ABS_REL_VELOCITY);
  max_delta_v_input_unit = U_NO_UNIT;
  DFLT_MAX_DELTA_V_UNIT = U_KM_PER_SEC;
  INT_MAX_DELTA_V_UNIT = U_AU_PER_DAY;

  prox_weights = 
    std::valarray<double>(DFLT_WEIGHT_VECTOR, N_KEPLERIAN_ELEMENTS);

  min_solar_distance = DFLT_MIN_SOLAR_DISTANCE;
  min_solar_distance_unit = U_NO_UNIT;
  DFLT_MIN_SOLAR_DISTANCE_UNIT = U_AU;
  INT_MIN_SOLAR_DISTANCE_UNIT = U_AU;

  des_solver = _NOT_SPECIFIED;
  max_error_rel = 0.0;
  max_error_abs = 0.0;
  impact_time = 0.0;
  impact_efficiency = 0.0;
  max_clock_jump = DFLT_MAX_CLOCK_JUMP;
  resolution_factor = DFLT_RESOLUTION_FACTOR;

  max_itgr_step_size = DFLT_MAX_ITGR_STEP_SIZE;
  min_itgr_step_size = DFLT_MIN_ITGR_STEP_SIZE;

  // initialize dynamic accuracy control (DAC) default parameters
  dac_usage = DFLT_DAC_USAGE;
  dac_factor_max = DFLT_DAC_FACTOR_MAX;
  dac_hc_start_size = DFLT_DAC_HC_START_SIZE;
  dac_hc_final_size = DFLT_DAC_HC_FINAL_SIZE;

  initDSSC();

} // void CInputSimulation::init(void)


//------------------------------------------------------------------------------
void CInputSimulation::initTargetOrbitalElementsFields(void)
{
  std::valarray<bool> t(DFLT_USE_TARGET_ORBITAL_ELEMENTS[0]
  , N_KEPLERIAN_ELEMENTS);
  use_target_orbital_elements.resize(N_KEPLERIAN_ELEMENTS, false);
  use_target_orbital_elements = t;


  target_orbital_elements.resize(N_KEPLERIAN_ELEMENTS);
  target_orbital_elements =
    std::valarray<double>(DFLT_TARGET_ORBITAL_ELEMENTS
    , N_KEPLERIAN_ELEMENTS);

  const int sz = target_orbital_elements.size();
  target_orbital_elements_input_unit.resize(sz, U_NO_UNIT);
  DFLT_TARGET_ORBITAL_ELEMENTS_UNIT.resize(sz, U_NO_UNIT);
  DFLT_TARGET_ORBITAL_ELEMENTS_UNIT[SEMIMAJORAXIS] = U_AU;
  DFLT_TARGET_ORBITAL_ELEMENTS_UNIT[ECCENTRICITY] = U_NO_UNIT;
  DFLT_TARGET_ORBITAL_ELEMENTS_UNIT[INCLINATION] = U_RADIAN;
  DFLT_TARGET_ORBITAL_ELEMENTS_UNIT[NODE] = U_RADIAN;
  DFLT_TARGET_ORBITAL_ELEMENTS_UNIT[PERIAPSIS] = U_RADIAN;
  DFLT_TARGET_ORBITAL_ELEMENTS_UNIT[ANOMALY] = U_RADIAN;

  // set internal units to standard input units
  INT_TARGET_ORBITAL_ELEMENTS_UNIT = DFLT_TARGET_ORBITAL_ELEMENTS_UNIT;


  max_target_orbital_element_deviation = 
    std::valarray<double>(DFLT_MAX_TARGET_ORBITAL_ELEMENT_DEVIATION
    , N_KEPLERIAN_ELEMENTS);

  const int sz2 = max_target_orbital_element_deviation.size();
  target_orbital_elements_dev_input_unit.resize(sz2, U_NO_UNIT);
  DFLT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT.resize(sz2, U_NO_UNIT);
  DFLT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT[SEMIMAJORAXIS] = U_AU;
  DFLT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT[ECCENTRICITY] = U_NO_UNIT;
  DFLT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT[INCLINATION] = U_RADIAN;
  DFLT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT[NODE] = U_RADIAN;
  DFLT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT[PERIAPSIS] = U_RADIAN;
  DFLT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT[ANOMALY] = U_RADIAN;

  // set internal units to standard input units
  INT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT = DFLT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT;

} // void initTargetOrbitalElementsFields(void)


//------------------------------------------------------------------------------
void CInputSimulation::initInitialOrbitStateFields(void)
{
  // initial orbital elements
  elements_of_initial_orbit.Set(DFLT_INIT_ORBITAL_ELEMENTS);
  elements_of_initial_orbit_input_unit.resize(N_KEPLERIAN_ELEMENTS, U_NO_UNIT);
  DFLT_ELEMENTS_OF_INITIAL_ORBIT_UNIT.resize(N_KEPLERIAN_ELEMENTS,U_NO_UNIT);
  DFLT_ELEMENTS_OF_INITIAL_ORBIT_UNIT[SEMIMAJORAXIS] = U_AU;
  DFLT_ELEMENTS_OF_INITIAL_ORBIT_UNIT[ECCENTRICITY] = U_NO_UNIT;
  DFLT_ELEMENTS_OF_INITIAL_ORBIT_UNIT[INCLINATION] = U_RADIAN;
  DFLT_ELEMENTS_OF_INITIAL_ORBIT_UNIT[NODE] = U_RADIAN;
  DFLT_ELEMENTS_OF_INITIAL_ORBIT_UNIT[PERIAPSIS] = U_RADIAN;
  DFLT_ELEMENTS_OF_INITIAL_ORBIT_UNIT[ANOMALY] = U_RADIAN;
} // void initInitialOrbitStateFields(void)


//------------------------------------------------------------------------------
void CInputSimulation::initInitialVectorStateFields(void)
{
  // initial vector state position vector
  setInitPosVec(0.0,0.0,0.0);
  const int sz1 = getInitPosVec().size();
  int i;
  initial_position_input_units.resize(sz1, U_NO_UNIT);
  DFLT_INITIAL_POSITION_UNITS.resize(sz1, U_NO_UNIT);
  for (i = 0; i < sz1; ++i) DFLT_INITIAL_POSITION_UNITS[i] = U_AU;

  // initial vector state velocity vector
  setInitVelVec(0.0,0.0,0.0);
  const int sz2 = getInitVelVec().size();
  initial_velocity_input_units.resize(sz2, U_NO_UNIT);
  DFLT_INITIAL_VELOCITY_UNITS.resize(sz2, U_NO_UNIT);
  for (i = 0; i < sz1; ++i) DFLT_INITIAL_VELOCITY_UNITS[i] = U_AU_PER_DAY;
} // void initInitialVectorStateFields(void)


//------------------------------------------------------------------------------
void CInputSimulation::initDSSC(void)
{
  // initialize dynamic step size control (DSSC) default parameters
  dssc_usage = DFLT_DSSC_USAGE;
  dssc_use_approach_control = DFLT_DSSC_USE_APPROACH_CONTROL;
  dssc_use_step_angle_control = DFLT_DSSC_USE_STEP_ANGLE_CONTROL;
  dssc_use_step_distance_control = DFLT_DSSC_USE_STEP_DISTANCE_CONTROL;
  dssc_approach_step_size_factor = DFLT_DSSC_APPROACH_STEP_SIZE_FACTOR;

  // max. allowed control step angle
  dssc_max_step_angle = DFLT_DSSC_MAX_STEP_ANGLE;
  DFLT_DSSC_MAX_STEP_ANGLE_UNIT = SpaceLib::U_RADIAN;
  INT_DSSC_MAX_STEP_ANGLE_UNIT = SpaceLib::U_RADIAN;

  // max. allowed control step linear distance
  dssc_max_step_distance = DFLT_DSSC_MAX_STEP_DISTANCE;
  INT_DSSC_MAX_STEP_DISTANCE_UNIT = SpaceLib::U_AU;
  DFLT_DSSC_MAX_STEP_DISTANCE_UNIT = SpaceLib::U_AU;

  // max. allowed time span of control step
  dssc_max_step_size = DFLT_DSSC_MAX_STEP_SIZE;
  DFLT_DSSC_MAX_STEP_SIZE_UNIT = SpaceLib::U_DAY;
  INT_DSSC_MAX_STEP_SIZE_UNIT = SpaceLib::U_DAY;

  // min. allowed time span of control step
  dssc_min_step_size = DFLT_DSSC_MIN_STEP_SIZE;
  DFLT_DSSC_MIN_STEP_SIZE_UNIT = SpaceLib::U_DAY;
  INT_DSSC_MIN_STEP_SIZE_UNIT = SpaceLib::U_DAY;

  dssc_approach_step_size_factor = DFLT_DSSC_APPROACH_STEP_SIZE_FACTOR;
} // void initDSSC(void)


//------------------------------------------------------------------------------
void CInputSimulation::initGTOC2(void)
{
  gtoc2_bodies.resize(4);
  for (unsigned int i=0;i<gtoc2_bodies.size();++i) gtoc2_bodies[i].reserve(350);
} // void CInputSimulation::initGTOC2(void)

//------------------------------------------------------------------------------
void CInputSimulation::readDynIntegrationReductionSetting(Global::ParamList& p)
throw (std::invalid_argument)
{
  const std::string method_name = 
"void CInputSimulation::readDynIntegrationReductionSetting(Global::ParamList&)";

  try
  {
    InputBoolParameter("DYN_INTEGRATION_INTERVAL"
      , dynamic_integration_interval_reduction, p, method_name
      , TRAILING_CHARS, true);
  }
  catch (std::invalid_argument&)
  {
    Warning("parameter \"DYN_INTEGRATION_INTERVAL\" not found. Using default setting.\n");
  }
} // void CInputSimulation::readDynIntegrationReductionSetting(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputSimulation::readFileIdentifier(Global::ParamList& p) 
throw (std::invalid_argument&)
{
  const std::string method_name = 
    "void CInputSimulation::readFileIdentifier(Global::ParamList& p)";

  const std::string err_msg1 = method_name + 
    std::string("not an InTrance 1.7 simulation parameter file\n");

  std::string s_identifier;

  InputStringParameter("FILE_TYPE", s_identifier, p, method_name,
    TRAILING_CHARS);
  removeChars(s_identifier, string(" \""));          // remove quotes and spaces
  upper(s_identifier);

  if (s_identifier!="INTRANCE1.7SIMULATIONPARAMETERFILE") throw
    std::invalid_argument(err_msg1);

} // void CInputInTrance::readFileIdentifier(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputSimulation::readInitStateBodyParams(Global::ParamList& p)
throw (std::invalid_argument&)
{
  const std::string method_name = 
    "void CInputSimulation::readInitStateBodyParams(Global::ParamList& p)";

  const std::string err_msg1 = method_name + 
    ": missing parameter (INTIAL_BODY_NAME)\n";

  string s_init_body;
  InputStringParameter("INITIAL_BODY_NAME", s_init_body, p, method_name
    , TRAILING_CHARS);

  // is the desired target a major body (planet, earth's moon)
  initial_body = m_solsys->getPlanetIndex(s_init_body);

  if (initial_body == 0)
  {
    // second, try to find an asteroid with the given name
    initial_body = m_solsys->insertAsteroid(s_init_body);
    if (initial_body == 0)
    {
      // and at last try comet names
      initial_body = m_solsys->insertComet(s_init_body);

      if (initial_body == 0) // error handling (no body found)
      {
        const std::string err_msg2 = method_name + 
          std::string(": couldn't find initial body (\"") + s_init_body + 
          std::string("\")\n");
        throw std::invalid_argument(err_msg2);
      } // if (initial_body == 0)
    } // if (initial_body == 0)
  } // if (initial_body == 0)
} // void CInputSimulation::readInitStateBodyParams(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputSimulation::readInitStateOrbitParams(Global::ParamList& p)
throw (std::invalid_argument&)
{
  const std::string method_name =
    "void CInputSimulation::readInitStateOrbitParams(Global::ParamList&)";

  try
  {
    InputStringParameter("INITIAL_CENTRAL_BODY", getInitCentralBodyName()
      , p, method_name, TRAILING_CHARS);
    upper(getInitCentralBodyName());
    setInitCentralBodyName(getInitCentralBodyName());
    // the following call to setInitCentralBodyName must be done to
    // concurrently set the ID of the central body, which must be obtained
    // during run time from the stored CSolarSystem object
  } catch (std::invalid_argument) {}

  for (int i = SEMIMAJORAXIS;i<N_KEPLERIAN_ELEMENTS;++i)
  {
    const KEPLERIAN_ELEMENT e = KEPLERIAN_ELEMENT(i);
    elements_of_initial_orbit_input_unit[e] = 
      InputParameter(INIT_ELEMENTS_PARAM_NAME[e]
    , elements_of_initial_orbit[e], p, method_name
      , TRAILING_CHARS, true);
    if (elements_of_initial_orbit_input_unit[e] != U_NO_UNIT)
    {
      elements_of_initial_orbit[e] = 
        elements_of_initial_orbit_input_unit[e]->convertTo( elements_of_initial_orbit[e], DFLT_ELEMENTS_OF_INITIAL_ORBIT_UNIT[e]);
    }
    else
    {
      elements_of_initial_orbit_input_unit[e] =
        DFLT_ELEMENTS_OF_INITIAL_ORBIT_UNIT[e];
    }
  }

} // void readInitStateOrbitParams(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputSimulation::readInitStateVectorParams(Global::ParamList& p)
throw (std::invalid_argument&)
{
  const std::string method_name =
    "void CInputSimulation::readInitStateVectorParams(Global::ParamList&)";

  try                                                   // input of central body
  {
    InputStringParameter("INITIAL_CENTRAL_BODY", initial_central_body_name
      , p, method_name, TRAILING_CHARS);
    // the following call to setInitCentralBodyName must be done to
    // concurrently set the ID of the central body, which must be obtained
    // during run time from the stored CSolarSystem object
    setInitCentralBodyName(getInitCentralBodyName());
  } catch (std::invalid_argument) {}

  unsigned int i;
  for (i = 0; i<3;++i)                    // input of position vector components
  {
   initial_position_input_units[i] =
     InputParameter(INIT_STATE_POS_VECTOR_PARAM_NAMES[i]
     , initial_position[i], p, method_name, TRAILING_CHARS, true);
   if (initial_position_input_units[i] != U_NO_UNIT)
   {
     initial_position[i] = initial_position_input_units[i]->convertTo(
       initial_position[i], DFLT_INITIAL_POSITION_UNITS[i] );
   }
   else
   {
     initial_position_input_units[i] = DFLT_INITIAL_POSITION_UNITS[i];
   }
  } // for(i=0;i<3;++i)


  for (i = 0; i<3; ++i)                  // input of velocity vector components
  {
   initial_velocity_input_units[i] = 
     InputParameter(INIT_STATE_VEL_VECTOR_PARAM_NAMES[i]
     , initial_velocity[i], p, method_name, TRAILING_CHARS, true);
   if (initial_velocity_input_units[i] != U_NO_UNIT)
   {
     initial_velocity[i] = initial_velocity_input_units[i]->convertTo(
       initial_velocity[i], DFLT_INITIAL_VELOCITY_UNITS[i] );
   }
   else
   {
     initial_velocity_input_units[i] = DFLT_INITIAL_VELOCITY_UNITS[i];
   }  
  } // for(i=0;i<3;++i)

} // void readInitStateVectorParams(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputSimulation::readInitialStateParameters(Global::ParamList& p)
throw (std::invalid_argument&)
{
  const std::string method_name =
    "void CInputSimulation::readInitialStateParameters(Global::ParamList&)";
  const std::string err_msg1 = method_name + 
    "invalid initial state type";

  switch (getInitialStateType())
  {

  case BODY_STATE:
    readInitStateBodyParams(p);
    break;

  case FREE_ORBIT_STATE:
    readInitStateOrbitParams(p);
    break;

  case VECTOR_STATE:
    readInitStateVectorParams(p);
    break;

  default:
    throw std::invalid_argument(err_msg1);
    break;

  } // switch (...)
} // void CInputSimulation::readInitialStateParameters(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputSimulation::readInitialStateType(Global::ParamList& p)
throw (std::invalid_argument)
{
  const std::string method_name = 
    "void CInputSimulation::readInitialStateType(Global::ParamList& p)";

  const std::string err_msg1 = method_name + 
    std::string(": couldn't find parameter \"INITIAL_STATE\"!\n");

  string s_init_state;
  InputStringParameter("INITIAL_STATE", s_init_state, p, err_msg1
    , TRAILING_CHARS);
  removeBlanks(s_init_state);upper(s_init_state);

  const int sz1 = sizeof(START_STATE_TYPES)/sizeof(START_STATE_TYPES[0]);
  const int sz2 = sizeof(START_STATE_TYPES_PARAM_VALUE)
    / sizeof(START_STATE_TYPES_PARAM_VALUE[0]);

  assert(sz1==sz2);

  int i = -1;
  while ( (++i) < sz1)
  {
    if (START_STATE_TYPES_PARAM_VALUE[i] == s_init_state)
    {
      setInitialStateType(START_STATE_TYPES[i]);
      break;
    } // if ...
  } // while ...
  if (i == sz1) 
  {
    const std::string err_msg2 = method_name + 
      std::string(": couldn't parse initial state type \"") +
      s_init_state + std::string("\"\n");

    throw std::invalid_argument(err_msg2);
  } // if (i == sz1)
} // void CInputSimulation::readInitialStateType(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputSimulation::readIntegrationInterval(Global::ParamList& p)
throw (std::invalid_argument&)
{
  const std::string method_name = 
    "void CInputSimulation::readIntegrationInterval(Global::ParamList&)";

  const std::string err_msg1 = method_name + 
    std::string(": missing parameter (INTEGRATION_INTERVAL)\n");

  integration_interval_input_unit = InputParameter("INTEGRATION_INTERVAL"
    , integration_interval, p, err_msg1, TRAILING_CHARS);
  if (integration_interval_input_unit != SpaceLib::U_NO_UNIT)
  {
  integration_interval = 
    integration_interval_input_unit->convertTo(integration_interval
    ,ITGR_INTVL_REF_UNIT);
  } else
  {
    integration_interval_input_unit = ITGR_INTVL_REF_UNIT;
  }
} // void CInputSimulation::readIntegrationInterval(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputSimulation::readIntegrationMaxAbsError(Global::ParamList& p)
  throw (std::invalid_argument&)
{
  const std::string method_name =
    "void InputSimulation::readIntegrationMaxAbsError(Global::ParamList& p)";

  InputParameter("MAX_ABSOLUTE_ERROR", max_error_abs, p, method_name,
    TRAILING_CHARS);

} // void InputSimulation::readIntegrationMaxAbsError(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputSimulation::readIntegrationMaxRelError(Global::ParamList& p)
  throw (std::invalid_argument&)
{
  const std::string method_name =
    "void InputSimulation::readIntegrationMaxRelError(Global::ParamList& p)";

  InputParameter("MAX_RELATIVE_ERROR", max_error_rel, p, method_name,
    TRAILING_CHARS);

} // void InputSimulation::readIntegrationMaxRelError(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputSimulation::readIntegrator(Global::ParamList& p)
throw (std::invalid_argument)
{
  const std::string method_name =
    "void CInputSimulation::readIntegrator(Global::ParamList&)";

  std::string s_des_solver;
  InputStringParameter("INTEGRATOR", s_des_solver,p,method_name,TRAILING_CHARS);
  removeChars(s_des_solver, " \"");
  upper(s_des_solver);

  const int sz1=sizeof(INTEGRATOR_TYPES)/sizeof(INTEGRATOR_TYPES[0]);
  const int sz2=sizeof(INTEGRATOR_TYPE_NAMES)/sizeof(INTEGRATOR_TYPE_NAMES[0]);

  // check that all known types have a literal parameter name and vice versa
  assert(sz1==sz2); 

  int i = -1;
  while ( (++i) < sz1 ) // loop over all known integrator types
  {
    if (s_des_solver == INTEGRATOR_TYPE_NAMES[i])  // we found the integrator
    {
      des_solver = INTEGRATOR_TYPES[i];
      break;
    } // if ...
  } // while ...
  if (i == sz1)     // error handling
  {
    const std::string err_msg1 = method_name +
      std::string(": couldn't determine integrator type (\"")+
      s_des_solver + std::string("\")\n");
    throw std::invalid_argument(err_msg1);
  } // if (i == sz1)

} // void CInputSimulation::readIntegrator(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputSimulation::readMaxIntegrationSteps(Global::ParamList& p)
throw (std::invalid_argument&)
{
  const std::string method_name = 
    "void CInputSimulation::readMaxIntegrationSteps(Global::ParamList&)";

  const std::string err_msg1 = method_name + 
    std::string(": missing parameter (INTEGRATION_STEPS)\n");

  InputParameter("INTEGRATION_STEPS", max_steps, p, err_msg1,
    TRAILING_CHARS);
} // void CInputSimulation::readMaxIntegrationSteps(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputSimulation::readMaxStartTime(Global::ParamList& p)
throw (std::invalid_argument&)
{
  const std::string method_name = 
    "void CInputSimulation::readMaxStartTime(Global::ParamList&)";

  const std::string err_msg1 = method_name + 
    std::string(": missing parameter (SIM_START_TIME_MAX)\n");

  InputParameter("SIM_START_TIME_MAX", max_start_time, p, err_msg1,
    TRAILING_CHARS);
} // void CInputSimulation::readStartTime(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputSimulation::readMinOutputPoints(Global::ParamList& p) 
throw (std::invalid_argument)
{
  const std::string method_name = 
    "void CInputSimulation::readMinOutputPoints(Global::ParamList&)";

  const std::string err_msg1 = method_name + 
    std::string(": missing parameter (MIN_OUTPUT_POINTS)\n");

  InputParameter("MIN_OUTPUT_POINTS", min_output_points, p, err_msg1,
    TRAILING_CHARS);
} // void CInputSimulation::readMinOutputPoints(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputSimulation::readMinSolarDistance(Global::ParamList& p) 
throw (std::invalid_argument)
{
  const std::string method_name =
    "void CInputSimulation::readMinSolarDistance(Global::ParamList&)";

  min_solar_distance_unit = InputParameter("MIN_SOLAR_DISTANCE"
    , min_solar_distance, p, method_name, TRAILING_CHARS);
  if (min_solar_distance_unit == U_NO_UNIT)
  {
    min_solar_distance_unit = DFLT_MIN_SOLAR_DISTANCE_UNIT;
  }
  if (min_solar_distance_unit != INT_MIN_SOLAR_DISTANCE_UNIT)
  {
    min_solar_distance = min_solar_distance_unit->convertTo(min_solar_distance,
     INT_MIN_SOLAR_DISTANCE_UNIT);
  }
} // void CInputSimulation::readMinSolarDistance(Global::ParamList&) 


//------------------------------------------------------------------------------
void CInputSimulation::readMinStartTime(Global::ParamList& p)
throw (std::invalid_argument&)
{
  const std::string method_name = 
    "void CInputSimulation::readMinStartTime(Global::ParamList&)";

  const std::string err_msg1 = method_name + 
    std::string(": missing parameter (SIM_START_TIME_MIN)\n");

  InputParameter("SIM_START_TIME_MIN", min_start_time, p, err_msg1,
    TRAILING_CHARS);
} // void CInputSimulation::readStartTime(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputSimulation::readAccuracyBoundaries(Global::ParamList& p)
throw (std::invalid_argument)
{

  const std::string method_name = 
    "void CInputSimulation::readAccuracyBoundaries(Global::Parameters& p)";
  const std::string err_msg1 = method_name + std::string("unkown trajectory type\n");

  // maximum final relative distance
  switch (trajectory)
  {
  case RENDEZVOUS:
  case FLYBY:
  case RIMPACT:
  case ACTPROBLEM:
  case SPIPROBLEM:         
  case GTOC2TRAJ:
    readMaxFinalRelDistParam(p);
    break;
  default:
    throw invalid_argument(err_msg1);
    break;
  }

  // maximum final relative velocity
  switch (trajectory)
  {
  case RENDEZVOUS:
  case ACTPROBLEM:
    readMaxFinalRelVelParam(p);
    break;
  case FLYBY:
  case RIMPACT:
  case SPIPROBLEM:
    max_dv = HUGE_VAL;
    break;
  case GTOC2TRAJ:
    readMaxFinalRelVelParam(p);
    break;
  default:
    throw invalid_argument(err_msg1);
    break;
  }
}


//------------------------------------------------------------------------------
void CInputSimulation::readMaxFinalRelDistParam(Global::ParamList& p)
{
  const std::string method_name = 
    "void CInputSimulation::readMaxFinalRelDistParam(Global::ParamList& p)";

  max_delta_r_input_unit = InputParameter("TARGET_DIST_MAX", max_dr, p
    , method_name, TRAILING_CHARS, true);
 
  if (max_delta_r_input_unit == U_NO_UNIT)
  {
    max_delta_r_input_unit = DFLT_MAX_DELTA_R_UNIT;
  }
  if (max_delta_r_input_unit != INT_MAX_DELTA_R_UNIT)
  {
    max_dr = max_delta_r_input_unit->convertTo(max_dr,
      INT_MAX_DELTA_R_UNIT);
  }
} // void readMaxFinalRelDistParam(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputSimulation::readMaxFinalRelVelParam(Global::ParamList& p)
{
  const std::string method_name =
    "void CInputSimulation::readMaxFinalRelVelParam(Global::ParamList& p)";

  max_delta_v_input_unit = InputParameter("TARGET_RELVEL_MAX", max_dv, p
    , method_name, TRAILING_CHARS, true);
  if (max_delta_v_input_unit == U_NO_UNIT)
  {
    max_delta_v_input_unit = DFLT_MAX_DELTA_V_UNIT;
  }
  if (max_delta_v_input_unit != INT_MAX_DELTA_V_UNIT)
  {
    max_dv = max_delta_v_input_unit->convertTo(max_dv,
      INT_MAX_DELTA_V_UNIT);
  }
} // void readMaxFinalRelVelParam(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputSimulation::readActProblemParameters(Global::ParamList& p)
  throw (std::invalid_argument&)
{
  const std::string method_name = 
    "void CInputSimulation::inputActProblemParameters(Global::ParamList&)";

  InputParameter("ACT_DISTANCE_FACTOR", act_distance_factor, p, method_name,
    TRAILING_CHARS, false);

  InputParameter("ACT_EXPONENT", act_exponent, p, method_name,
    TRAILING_CHARS, false);
} // void inputActProblemParameters(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputSimulation::readGTOCOMP2Parameters(Global::ParamList& p)
throw (std::invalid_argument&)
{
  initGTOC2();

  using namespace std;
  string fn,name;
  double elements[6], epoch, mean_motion;
  const double _gm = m_solsys->GM_AUd(m_solsys->getPlanetIndex(SUN_NAME));
  int group,i, n_bodies=0;
  const int bufsize = 128;
  char buf[bufsize];
  CBody* asteroid = 0;
  CBody* b_sun = m_solsys->getBody(m_solsys->getPlanetIndex(SUN_NAME));
  BodyID id_temp = 0;

  InputFilenameParameter("TARGET_FILE",fn, p, "",TRAILING_CHARS, "",true,true);
  ifstream is(fn.c_str()); // so we should not need error checks here    
  is.getline(buf,bufsize, '\n');
  Output("\nGTOC2: Reading orbit parameters of target bodies");
  while (!is.eof())
  {
    is >> name;name = string("B")+name;
    for (i = SEMIMAJORAXIS; i<N_KEPLERIAN_ELEMENTS; ++i) is >> elements[i];
    for (i = INCLINATION; i<N_KEPLERIAN_ELEMENTS; ++i) elements[i] *= DEG2RAD;
    is >> epoch;
    is >> group;
    mean_motion = CalcMeanMotion(elements[SEMIMAJORAXIS], _gm);

    asteroid = new CBody(name, 1.0, epoch
     ,elements[SEMIMAJORAXIS],elements[ECCENTRICITY],elements[INCLINATION]
     ,elements[NODE]
     ,elements[NODE]+elements[PERIAPSIS]                   // curlypi
     ,elements[NODE]+elements[PERIAPSIS]+elements[ANOMALY] // longitude L
     ,0,0,0,0,0, mean_motion, KEPLERM, b_sun);
    if (asteroid!=0)
    {
     asteroid->setCenterBodyId(SUN);                             // heliocentric

     // insert asteroid into solar system object
     id_temp = m_solsys->insertBody(asteroid, true
       , m_solsys->NumberOfPlanets()+1);

     // and in problem specific data structure
     (gtoc2_bodies[group-1]).push_back(asteroid);
    } else //if (asteroid!=0)
    {
      ostringstream os;
      Output("\ncouldn't insert all bodies!\n");
      os << n_bodies;
      Output(os.str() + " inserted!\n");
    } // else
    if ( ((++n_bodies)%50) == 0) Output(".");
  } // while (!is.eof())
  Output("\n");
} // void CInputSimulation::readGTOCOMP2Parameters(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputSimulation::readDACParameters(ParamList& p) 
throw (std::invalid_argument)
{
  const std::string s = "void CInputSimulation::inputDACParameters(\
                        ParamList& p) throw(std::invalid_argument)";

  InputBoolParameter("USE_DAC", 
    dac_usage, p, s, TRAILING_CHARS, false);

  if (dac_usage==true)
  {
    InputParameter("DAC_MAX_FACTOR", dac_factor_max,
      p, s, TRAILING_CHARS, false);

    InputParameter("DAC_HC_START_SIZE", dac_hc_start_size,
      p, s, TRAILING_CHARS, false);

    InputParameter("DAC_HC_FINAL_SIZE", dac_hc_final_size,
      p, s, TRAILING_CHARS, false);
  }
} // void CInputSimulation::inputDACParameters(ParamList&) 


//------------------------------------------------------------------------------
void CInputSimulation::readDeflectionParameters(Global::ParamList& p) 
throw (std::invalid_argument)
{
  const std::string method_name = 
    "void CInputSimulation::inputDeflectionParameters(Global::ParamList&)";

  // read needed parameter of time of undisturbed asteroid's impact on target body
  InputParameter("IMPACT_TIME", impact_time, p, method_name, TRAILING_CHARS);

  // read needed parameter of hitting s/c's impact efficiency
  InputParameter("IMPACT_EFFICIENCY", impact_efficiency, p, method_name
    , TRAILING_CHARS);

} // void CInputSimulation::inputDeflectionParameters(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputSimulation::readDSSCParameters(ParamList& p)
throw(std::invalid_argument)
{
  const std::string s = "void CInputSimulation::inputDSSCParameters(\
                        ParamList& p) throw(std::invalid_argument)";

  InputBoolParameter("USE_DSSC", 
    dssc_usage, p, s, TRAILING_CHARS, false);

  if (dssc_usage == true)         // read needed parameter values
  {
    // total velocity step size control ?
    InputBoolParameter("DSSC_STEP_DISTANCE_CONTROL"
      , dssc_use_step_distance_control, p, s, TRAILING_CHARS, false);

    if (dssc_use_step_distance_control == true)
    {
      dssc_max_step_distance_input_unit = 
        InputParameter("DSSC_MAX_STEP_DISTANCE", dssc_max_step_distance,
        p, s, TRAILING_CHARS, true);
      if (dssc_max_step_distance_input_unit == SpaceLib::U_NO_UNIT)
      {
        dssc_max_step_distance_input_unit = DFLT_DSSC_MAX_STEP_DISTANCE_UNIT;
      }
      if (dssc_max_step_distance_input_unit!=INT_DSSC_MAX_STEP_DISTANCE_UNIT)
      {
        dssc_max_step_distance = 
          dssc_max_step_distance_input_unit->convertTo(dssc_max_step_distance,
          INT_DSSC_MAX_STEP_DISTANCE_UNIT);
      }
    } // if (dssc_use_step_distance_control == true)

    // and/or orbital angular velocity control ?
    InputBoolParameter("DSSC_STEP_ANGLE_CONTROL"
      , dssc_use_step_angle_control, p, s, TRAILING_CHARS, false);

    if (dssc_use_step_angle_control == true)
    {
      dssc_max_step_angle_input_unit = InputParameter("DSSC_MAX_STEP_ANGLE"
        ,dssc_max_step_angle, p, s, TRAILING_CHARS, true);
      if (dssc_max_step_angle_input_unit == SpaceLib::U_NO_UNIT)
      {
        dssc_max_step_angle_input_unit=DFLT_DSSC_MAX_STEP_ANGLE_UNIT;
      }
      if (dssc_max_step_angle_input_unit != INT_DSSC_MAX_STEP_ANGLE_UNIT)
      {
        dssc_max_step_angle = dssc_max_step_angle_input_unit->convertTo(
          dssc_max_step_angle, INT_DSSC_MAX_STEP_ANGLE_UNIT);
      }
    }

    // target approach step size control
    InputBoolParameter("DSSC_APPROACH_CONTROL"
      , dssc_use_approach_control, p, s, TRAILING_CHARS, false);

    if (dssc_use_approach_control == true)
    {
      InputParameter("DSSC_APPROACH_STEP_SIZE_FACTOR",
        dssc_approach_step_size_factor,p,s,TRAILING_CHARS, true);
    }

    // maximum dynamic step sizes
    dssc_max_step_size_input_unit = 
      InputParameter("DSSC_MAX_STEP_SIZE", dssc_max_step_size,
      p, s, TRAILING_CHARS, true);
    if (dssc_max_step_size_input_unit == SpaceLib::U_NO_UNIT)
    {
      dssc_max_step_size_input_unit=DFLT_DSSC_MAX_STEP_SIZE_UNIT;
    }
    if (dssc_max_step_size_input_unit != INT_DSSC_MAX_STEP_SIZE_UNIT)
    {
      dssc_max_step_size = dssc_max_step_size_input_unit->convertTo(
        dssc_max_step_size, INT_DSSC_MAX_STEP_SIZE_UNIT);
    }

    // minimum dynamic step sizes
    dssc_min_step_size_input_unit = 
      InputParameter("DSSC_MIN_STEP_SIZE", dssc_min_step_size,
      p, s, TRAILING_CHARS, true);
    if (dssc_min_step_size_input_unit == SpaceLib::U_NO_UNIT)
    {
      dssc_min_step_size_input_unit=DFLT_DSSC_MIN_STEP_SIZE_UNIT;
    }
    if (dssc_min_step_size_input_unit != INT_DSSC_MIN_STEP_SIZE_UNIT)
    {
      dssc_min_step_size = dssc_min_step_size_input_unit->convertTo(
        dssc_min_step_size, INT_DSSC_MIN_STEP_SIZE_UNIT);
    }
  }
} // void CInputSimulation::inputDSSCParameters(ParamList& p)

//------------------------------------------------------------------------------
void CInputSimulation::modifyDSSCParameters(void)
{
  const string method_name = "void CInputSimulation::modifyDSSCParameters(void)";
  if (getMaxFinalRelDistance()>0.0)
  {
    const double val = (getDSSCMinStepSize()/(getMaxFinalRelDistance()/2.0));
    switch (getTrajectoryType())
    {
    case RENDEZVOUS:
    case FLYBY:
    case RIMPACT:
    case ACTPROBLEM:
      if (getDSSCApproachFactor() > val)
      {
        setDSSCApproachFactor(val);
        Warning(method_name + string(": approach factor had to be modified \
                                     because of too large value!\n"));
      }
    default:break;
    }
  }
} // void CInputSimulation::modifyDSSCParameters(void)

//------------------------------------------------------------------------------
void CInputSimulation::readHyperbolicExcessVelocityParameters(
  Global::ParamList& p) throw (std::invalid_argument)
{
  const std::string method_name = "void CInputSimulation::\readHyperbolicExcessVelocityParameters(\
Global::ParamList& p)";

  int i;
  for (i=0;i<N_VINF_PARAMS;++i)
  {
    launch_v_inf_input_unit[i] = InputParameter(INIT_VINF_PARAM_NAMES[i],
      launch_v_inf[i], p, method_name, TRAILING_CHARS, false);
    if (launch_v_inf_input_unit[i] == U_NO_UNIT)
    {
      launch_v_inf_input_unit[i] = DFLT_V_INF_INPUT_UNIT[i];
    }
    if (launch_v_inf_input_unit[i] != INT_V_INF_UNIT[i])
    {
      launch_v_inf[i] = launch_v_inf_input_unit[i]->convertTo(
        launch_v_inf[i], INT_V_INF_UNIT[i]);
    }
  } // for (i=0;i<N_VINF_PARAMS;++i)

  // assert azimuth boundaries
  launch_v_inf[VINF_AZI_MAX] = std::min(PI-EPS, launch_v_inf[VINF_AZI_MAX]);
  launch_v_inf[VINF_AZI_MIN] = std::max(-PI+EPS, launch_v_inf[VINF_AZI_MIN]);

  // assert elevation boundaries
  launch_v_inf[VINF_ELE_MAX]= std::min(HALFPI-EPS,launch_v_inf[VINF_ELE_MAX]);
  launch_v_inf[VINF_ELE_MIN]= std::max(-HALFPI+EPS,launch_v_inf[VINF_ELE_MIN]);

} // void CInputSimulation::inputHyperbolicExcessVelocityParameters(ParamList&)


//------------------------------------------------------------------------------
///< read possible gravity assist bodies
void CInputSimulation::readGravitiyAssistBodies(Global::ParamList& p)
throw (std::invalid_argument)
{
  ParamList::iterator i;
  ga_body.resize(NEPTUNE+1);
  for (int loop=MERCURY; loop<=NEPTUNE; loop++) ga_body[loop]=false;
  ostringstream t;
  string st2; int idx;
  for (idx=0; idx<number_of_ga_bodies; idx++)
  {
    st2="GRAVITY_ASSIST_BODY"; t << (idx+1);
    st2 += t.str();
    i = find_if(p.begin(), p.end(), name_equal<Param>(st2));
    if (i!=p.end()) 
    {
      removeTrailingChars((*i).second, TRAILING_CHARS); // remove not needed chars
      removeChars((*i).second, string(" \""));          // remove possible quotes
      ga_body[m_solsys->getPlanetIndex((*i).second)]=true;
    }
    else throw invalid_argument(string("CInputSimulation::CInputSimulation(): missing parameter (")+st2+string(")\n"));
    st2.erase();t.str(string(""));
    p.erase(i);
  }
} // void CInputSimulation::inputGravitiyAssistBodies(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputSimulation::readOptimizationGoal(Global::ParamList& p)
    throw (std::invalid_argument)
{
  const std::string method_name =
    "void CInputSimulation::inputOptimizationGoal(Global::ParamList& p)";
  const std::string err_msg1 = method_name + 
    std::string("not all optimization goals do have parameters\n");
  const std::string err_msg2 = method_name + 
    std::string("unknown optimization goal\n");

  string og;
  InputStringParameter("OPTIMIZATION_GOAL", og, p, method_name, TRAILING_CHARS);
  removeChar(og, ' ');                                 // remove possible spaces
  upper(og);                              // convert all charaters to upper case

  int sz1 = sizeof(OPT_GOAL_TYPES)/sizeof(OPT_GOAL_TYPES[0]);
  if (sz1!=N_OPT_GOAL_TYPES) throw std::invalid_argument(err_msg1);

  int i = -1;
  while ((++i)<sz1)
  {
    if (og==OPT_GOAL_TYPE_PARAM_VALUE[i])
    {
      opt_goal = OPT_GOAL_TYPES[i];
      break;
    }
  }
  if (i==N_OPT_GOAL_TYPES) throw std::invalid_argument(err_msg2);

  switch (opt_goal)
  {
  case TRANSFER_TIME:
  case PROPELLANT_MASS:
    break;

  case MAXIMUM_DEFLECTION:
    readDeflectionParameters(p);
    break;

  case RETROGRADE_IMPACT_CRANKING:
  case ACT_DEFLECTION:
    readActProblemParameters(p);
    break;

  case GTOCOMP2:
    readGTOCOMP2Parameters(p);
    break;

  case SOLAR_POLAR_IMAGER:
    break;

  default:
    throw std::invalid_argument(err_msg1);
    break;
  }
} // CInputSimulation::inputOptimizationGoal(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputSimulation::readTargetBodySettings(Global::ParamList& p)
  throw (std::invalid_argument)
{
  const std::string method_name = 
    "void CInputSimulation::inputTargetBodySettings(Global::ParamList&)";
  const std::string err_msg1 = method_name + 
    std::string(": Target body not found!\n");
  const std::string err_msg2 = method_name + 
    std::string(": Lagrange point object could not be inserted into solar \
                system.\n");
  const std::string err_msg3 = method_name +
    std::string(": Error during creation of Lagrange - point object.\n");
  const std::string err_msg4 = method_name +
    std::string(": Wrong Lagrange point syntax (bodies not in solar \
                system).\n");
  const std::string err_msg5 = method_name +
    std::string(": Wrong Lagrange point syntax (to few seperators).\n");

  string targetBody;
  InputStringParameter("TARGET_BODY_NAME", targetBody, p, method_name,
    TRAILING_CHARS);
  removeChars(targetBody, " ");

  // Determines whether the target is a Langrangian equilibrium point and
  // if it's not, try to insert the desired "real" target body into the
  // solar system object.
  string ts(targetBody.substr(0,2));// cut first two characters from target name
  if ((ts != "L1") && (ts != "L2") && (ts != "L3") && (ts != "L4") 
    && (ts != "L5")) // not a Lagrange point
  {
    target_body=m_solsys->getPlanetIndex(targetBody);
    if (target_body == 0)
    {
      target_body = m_solsys->insertAsteroid(targetBody);
      if (target_body == 0)
      {
        target_body = m_solsys->insertComet(targetBody);
        if (target_body == 0) throw invalid_argument(err_msg1);
      }
    }
  }
  else                                     // it is a Lagrange equilibrium point
  {
    // Now we have to parse the given L-point because it contains two 
    // system determining body names.

    // count number of seperating characters
    string::difference_type c=count(targetBody.begin(), targetBody.end(), '_');

    if (c == 2)
    {
      int which;
      istringstream(targetBody.substr(1,1)) >> which;
      targetBody.erase(0, targetBody.find('_')+1);
      string b1(targetBody.substr(0,targetBody.find('_')));
      targetBody.erase(0,targetBody.find('_')+1);
      string b2(targetBody);

      BodyID bID1 = m_solsys->getPlanetIndex(b1);  // get body IDs
      BodyID bID2 = m_solsys->getPlanetIndex(b2);

      if ( (bID1!=0) && (bID2!=0)) // are those IDs valid?
      {
        CBody* bp1=m_solsys->getBody(bID1);  // get pointers to CBody objects
        CBody* bp2=m_solsys->getBody(bID2);

        CLagrange* L;
        L = new CLagrange(bp1, bp2, which); // make new object

        if (L!=0)
        {
          target_body=m_solsys->insertBody(L, true, 100); // insert it

          if (target_body==0) throw invalid_argument(err_msg2);

        } else throw invalid_argument(err_msg3);
      } else throw invalid_argument(err_msg4);
    } else throw invalid_argument(err_msg5);

  } // end (if a Lagrange point)
} // void CInputSimulation::inputTargetBodySettings(Global::ParamList&)


////------------------------------------------------------------------------------
//void CInputSimulation::readGTOC2TargetBodySettings(Global::ParamList& p)
//throw (std::invalid_argument&)
//{
//  using namespace std;
//  string fn,name;
//  double elements[6], epoch, mean_motion;
//  const double _gm = m_solsys->GM_AUd(m_solsys->getPlanetIndex(SUN_NAME));
//  int group,i, n_bodies=0;
//  const int bufsize = 128;
//  char buf[bufsize];
//  CBody* asteroid = 0;
//  CBody* b_sun = m_solsys->getBody(m_solsys->getPlanetIndex(SUN_NAME));
//
//  InputFilenameParameter("TARGET_FILE",fn, p, "",TRAILING_CHARS, "",true,true);
//  ifstream is(fn.c_str()); // so we should not need error checks here    
//  is.getline(buf,bufsize, '\n');
//  Output("\nGTOC2: Reading orbit parameters of target bodies");
//  while (!is.eof())
//  {
//    is >> name;
//    for (i = SEMIMAJORAXIS; i<N_KEPLERIAN_ELEMENTS; ++i) is >> elements[i];
//    for (i = INCLINATION; i<N_KEPLERIAN_ELEMENTS; ++i) elements[i] *= DEG2RAD;
//    is >> epoch;
//    is >> group;
//    mean_motion = CalcMeanMotion(elements[SEMIMAJORAXIS], _gm);
//
//    asteroid = new CBody(name, 1.0, epoch
//     ,elements[SEMIMAJORAXIS],elements[ECCENTRICITY],elements[INCLINATION]
//     ,elements[NODE]
//     ,elements[NODE]+elements[PERIAPSIS]                   // curlypi
//     ,elements[NODE]+elements[PERIAPSIS]+elements[ANOMALY] // longitude L
//     ,0,0,0,0,0, mean_motion, KEPLERM, b_sun);
//    if (asteroid!=0)
//    {
//     asteroid->setCenterBodyId(SUN);
//     m_solsys->insertBody(asteroid, true, 2000);
//    } else //if (asteroid!=0)
//    {
//      ostringstream os;
//      Output("\ncouldn't insert all bodies!\n");
//      os << n_bodies;
//      Output(os.str() + " inserted!\n");
//    } // else
//    if ( ((++n_bodies)%50) == 0) Output(".");
//  } // while (!is.eof())
//  Output("\n");
//} // void CSimulation::readGTOC2TargetBodySettings(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputSimulation::readTargetBodyOrbitSettings(Global::ParamList& p)
  throw (std::invalid_argument)
{
  const std::string method_name = 
    "void CInputSimulation::readTargetBodyOrbitSettings(Global::ParamList&)";
  const std::string err_msg1 = method_name + 
    std::string(": Target body not found!\n");

  string targetBody;
  InputStringParameter("TARGET_ORBIT_NAME", targetBody, p, method_name,
    TRAILING_CHARS);
  removeChars(targetBody, " ");

  setTargetBodyId(getSolarSystem()->getPlanetIndex(targetBody));
  if (getTargetBodyId() == 0)
  {
    setTargetBodyId(m_solsys->insertAsteroid(targetBody));
    if (getTargetBodyId() == 0)
    {
      setTargetBodyId(m_solsys->insertComet(targetBody));
      if (getTargetBodyId() == 0) throw invalid_argument(err_msg1);
    } // if (getTargetBodyId() == 0)
  } // if (getTargetBodyId() == 0)
} // void readTargetBodyOrbitSettings(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputSimulation::readTargetFreeOrbitSettings(Global::ParamList& p)
  throw (std::invalid_argument)
{
  const std::string method_name =
    "void CInputSimulation::inputTargetFreeOrbitSettings(Global::ParamList&) ";

  for ( int i = SEMIMAJORAXIS ; i < N_KEPLERIAN_ELEMENTS ; ++i )
  {
    try
    {
      // read desired value of target's semimajor axis
      target_orbital_elements_input_unit[i] =
        InputParameter(TARGET_ELEMENT_PARAM_NAME[i], target_orbital_elements[i], 
        p, method_name, TRAILING_CHARS, true);
      if (target_orbital_elements_input_unit[i] != U_NO_UNIT)
      {
        target_orbital_elements[i] =
          target_orbital_elements_input_unit[i]->convertTo(
          target_orbital_elements[i],DFLT_TARGET_ORBITAL_ELEMENTS_UNIT[i]);
      } // if (target_orbital_elements_input_unit[i] != U_NO_UNIT)
      else
      {
        target_orbital_elements_input_unit[i] =
          DFLT_TARGET_ORBITAL_ELEMENTS_UNIT[i];
      } // else
      use_target_orbital_elements[i] = true; // set flag
    } // try
    catch (std::invalid_argument)
    {
       use_target_orbital_elements[i] = false;   
    }
  } // for ( int i = SEMIMAJORAXIS ; i < N_KEPLERIAN_ELEMENTS ; ++i )

} // void CInputSimulation::inputTargetFreeOrbitSettings(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputSimulation::readTargetTypeSettings(Global::ParamList& p)
throw (std::invalid_argument)
{
  switch (target_state)
  {

  case BODY_STATE:
    if (getTrajectoryType()!=GTOC2TRAJ) readTargetBodySettings(p);
    break;

  case BODY_ORBIT_STATE:
    readTargetBodyOrbitSettings(p); 
    break;
 
  case FREE_ORBIT_STATE:
    readTargetFreeOrbitSettings(p);
    break;

  case ORBITAL_ELEMENTS_STATE:
    readTargetOrbitalElementsSettings(p);
    break;

  case RESERVED_01:break;

  default:
    // target state type unknown, so throw an exception
    throw invalid_argument(string("CInputSimulation::CInputSimulation():\
                                  unknown target state type\n"));
    break;
  }
} // void CInputSimulation::inputTargetTypeSettings(ParamList& p)


//------------------------------------------------------------------------------
void CInputSimulation::readTargetType(Global::ParamList& p) 
throw (std::invalid_argument)
{
  const std::string method_name = 
    "CInputSimulation::inputTarget(Global::ParamList& p)";
  const std::string err_msg_0 = method_name +
    std::string(": missing parameter (TARGET_STATE)\n");
  const std::string err_msg_1 = method_name +
    std::string(": unknown target state type\n");

  string s_target;
  InputStringParameter("TARGET_STATE", s_target, p, err_msg_0, TRAILING_CHARS);
  upper(s_target);removeChars(s_target, " \"");
  
  int sz1 = 
    sizeof(OPTIMIZATION_PROBLEM_DESC)/sizeof(OPTIMIZATION_PROBLEM_DESC[0]);
  int sz2 = 
    sizeof(TARGET_STATE_TYPES)/sizeof(TARGET_STATE_TYPES[0]);
  int sz3 = 
    sizeof(TARGET_TRAJECTORIES)/sizeof(TARGET_TRAJECTORIES[0]);
  // assure that all constant arrays have the same dimension
  assert( (sz1==sz2) && (sz1==sz3) && (sz2==sz3) );

  int i = -1;                      // -1 due to following postincrement operator
  while ( (++i) < sz1 )            // loop over possible target type desciptions
  {
    if (s_target==OPTIMIZATION_PROBLEM_DESC[i])
    {
      target_state = TARGET_STATE_TYPES[i];   // set problem specific state type
      trajectory   = TARGET_TRAJECTORIES[i]; // problem specific trajectory type
      break;
    }
  } // while loop
  if (i==sz1) throw invalid_argument(err_msg_1); // error exception
} // void inputTarget(Global::ParamList&) throw (std::invalid_argument)


//------------------------------------------------------------------------------
void CInputSimulation::readParameter(
  Global::ParamList& p
 ,std::string name
 ,double& val
 ,PhysUnit& input_unit
 ,PhysUnit& default_unit
 ,PhysUnit& internal_unit
) throw (std::invalid_argument&, std::logic_error&)
{
  const std::string method_name = "CInputSimulation::readParameter(...)";

  input_unit = InputParameter(name, val, p, method_name, TRAILING_CHARS, true);

  if (input_unit == U_NO_UNIT)
  {
    input_unit = default_unit;
  } // if (input_unit == U_NO_UNIT)

  if (input_unit != internal_unit)
  {
    val = input_unit->convertTo(val, internal_unit);        
  } // if (input_unit != internal_unit)

} // void readParameter(...)


//------------------------------------------------------------------------------
void CInputSimulation::readTargetOrbitalElementsSettings(Global::ParamList& p)
throw (std::invalid_argument&, std::logic_error&)
{
  const std::string method_name = "void CInputSimulation::readTargetOrbital\
ElementsSettings(Global::ParamList& p)";

  const std::string err_msg1 = method_name + std::string(": \
one of the needed data structures is of incompatible size!\n");

  const std::string err_msg2 = method_name + std::string(": \
not a single target orbital element given!\n");

  // validity check 
  std::valarray<size_t> sz(size_t(0),11);
  sz[0]=sizeof(TARGET_ELEMENT_PARAM_NAME)/sizeof(TARGET_ELEMENT_PARAM_NAME[0]);
  sz[1]=target_orbital_elements.size();
  sz[2]=use_target_orbital_elements.size();
  sz[3]=max_target_orbital_element_deviation.size();
  sz[4]=sizeof(TARGET_ELEMENT_MAX_DEVIATION_PARAM_NAME)/
    sizeof(TARGET_ELEMENT_MAX_DEVIATION_PARAM_NAME[0]);
  sz[5]=target_orbital_elements_input_unit.size();
  sz[6]=DFLT_TARGET_ORBITAL_ELEMENTS_UNIT.size();
  sz[7]=INT_TARGET_ORBITAL_ELEMENTS_UNIT.size();
  sz[8]=target_orbital_elements_dev_input_unit.size();
  sz[9]=DFLT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT.size();
  sz[10]=INT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT.size();

  sz-=N_KEPLERIAN_ELEMENTS;
  if (sz.sum() != 0) throw std::logic_error(err_msg1);
  // end of validity check

  setTargetOrbitalElementMaxDeviation(SEMIMAJORAXIS,
    INT_MAX_DELTA_R_UNIT->convertTo(max_dr
    ,INT_TARGET_ORBITAL_ELEMENTS_UNIT[SEMIMAJORAXIS]));

  for ( int i = SEMIMAJORAXIS ; i < N_KEPLERIAN_ELEMENTS ; ++i )
  {
    try
    {
      readParameter(p
       , TARGET_ELEMENT_PARAM_NAME[i]
       , target_orbital_elements[i]
       , target_orbital_elements_input_unit[i]
       , DFLT_TARGET_ORBITAL_ELEMENTS_UNIT[i]
       , INT_TARGET_ORBITAL_ELEMENTS_UNIT[i]);

      readParameter(p
       , TARGET_ELEMENT_MAX_DEVIATION_PARAM_NAME[i]
       , max_target_orbital_element_deviation[i]
       , target_orbital_elements_dev_input_unit[i]
       , DFLT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT[i]
       , INT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT[i]);

      use_target_orbital_elements[i] = true; // set flag

      //// read weight factor
      //InputParameter(PROX_WEIGHT_PARAM_NAME[i], prox_weights[i], p
      //  , method_name, TRAILING_CHARS, false);
    } // try
    catch (std::invalid_argument)
    {
       use_target_orbital_elements[i] = false;   
    } // catch (std::invalid_argument)
  } // for ( int i = SEMIMAJORAXIS ; i < N_KEPLERIAN_ELEMENTS ; ++i )
  // check whether at least one target parameter is given and throw
  // an exception to indicate an error
  unsigned int u = 0;
  while ((u<N_KEPLERIAN_ELEMENTS) &&(use_target_orbital_elements[u]==false)) ++u;
  if (u == N_KEPLERIAN_ELEMENTS) throw std::invalid_argument(err_msg2);
} // void inputTargetOrbitalElementsSettings(Global::ParamList& p)


//------------------------------------------------------------------------------
CInputSimulation::CInputSimulation(ParamList& parameters,
                                   SpaceLib::CSolarSystem* solsys)
{
  const string prefix("CInputSimulation::\
                      CInputSimulation(ParamList&, CSolarSystem*):");
  init();                                                  // init member fields
  assert(solsys!=0);
  setSolarSystem(solsys);

  readFileIdentifier(parameters);
  readMinStartTime(parameters);
  readMaxStartTime(parameters);
  readIntegrationInterval(parameters);
  readMaxIntegrationSteps(parameters);
  readMinOutputPoints(parameters);
  readDynIntegrationReductionSetting(parameters);
  readInitialStateType(parameters);
  readInitialStateParameters(parameters);
  readHyperbolicExcessVelocityParameters(parameters);

  // read potential gravity assist bodies
  InputParameter("NUMBER_OF_GRAVITY_ASSISTS", number_of_ga_bodies,
    parameters, prefix, TRAILING_CHARS, false);
  if (number_of_ga_bodies == 0) ga_body.resize(0);
  else readGravitiyAssistBodies(parameters);

  // estimate needed target state and trajectory type
  readTargetType(parameters); 

  // read max. delta_r and max. delta_v (problem dependend)
  readAccuracyBoundaries(parameters);

  // read/parse target settings according to prior input of target state type
  readTargetTypeSettings(parameters); 

  readOptimizationGoal(parameters); // read desired optimization goal
  readMinSolarDistance(parameters); 
  readIntegrator(parameters);
  readIntegrationMaxAbsError(parameters);
  readIntegrationMaxRelError(parameters);
  readDACParameters(parameters);  // dynamic accuracy control parameters
  readDSSCParameters(parameters); // dynamic step size control parameters
  modifyDSSCParameters();

};


//------------------------------------------------------------------------------
const CInputSimulation& CInputSimulation::operator=(const CInputSimulation& rhs)
{
  gtoc2_bodies = rhs.gtoc2_bodies;

  m_solsys=rhs.m_solsys;
  setInitialStateType(rhs.getInitialStateType());

  // set the name of initial central body and, afterwards, the ID of this
  // central body. Therefore there's no need to call "setInitialCentralBody"
  // because it's is implicitely called.
  setInitCentralBodyName(rhs.getInitCentralBodyName());

  setHeliocentricSimulation(rhs.isHeliocentric());

  // initial position and velocity vector
  setInitPosVec(rhs.getInitPosVec());
  initial_position_input_units=rhs.initial_position_input_units;
  DFLT_INITIAL_POSITION_UNITS=rhs.DFLT_INITIAL_POSITION_UNITS;
  setInitVelVec(rhs.getInitVelVec());
  initial_velocity_input_units=rhs.initial_velocity_input_units;
  DFLT_INITIAL_VELOCITY_UNITS = rhs.DFLT_INITIAL_VELOCITY_UNITS;

  act_distance_factor = rhs.act_distance_factor;
  act_exponent = rhs.act_exponent;
  des_solver = rhs.des_solver;
  dynamic_integration_interval_reduction = 
    rhs.dynamic_integration_interval_reduction;
  ga_body = rhs.ga_body;
  impact_efficiency = rhs.impact_efficiency;
  impact_time = rhs.impact_time ;
  min_start_time = rhs.min_start_time;
  max_start_time = rhs.max_start_time;
  integration_interval = rhs.integration_interval;  
  max_steps = rhs.max_steps;
  min_output_points = rhs.min_output_points;

  setInitialStateType(rhs.getInitialStateType());
  initial_body = rhs.initial_body;
  elements_of_initial_orbit = rhs.elements_of_initial_orbit;
  elements_of_initial_orbit_input_unit =
    rhs.elements_of_initial_orbit_input_unit;
  DFLT_ELEMENTS_OF_INITIAL_ORBIT_UNIT =
    rhs.DFLT_ELEMENTS_OF_INITIAL_ORBIT_UNIT;

   // initial excess velocity parameters
  launch_v_inf = rhs.launch_v_inf;
  launch_v_inf_input_unit = rhs.launch_v_inf_input_unit;
  DFLT_V_INF_INPUT_UNIT = rhs.DFLT_V_INF_INPUT_UNIT;
  INT_V_INF_UNIT = rhs.INT_V_INF_UNIT;

  trajectory = rhs.trajectory;
  opt_goal = rhs.opt_goal;

  number_of_ga_bodies = rhs.number_of_ga_bodies;

  target_state = rhs.target_state;

  use_target_orbital_elements = rhs.use_target_orbital_elements;

  target_orbital_elements = rhs.target_orbital_elements;

  target_orbital_elements_input_unit=rhs.target_orbital_elements_input_unit;

  DFLT_TARGET_ORBITAL_ELEMENTS_UNIT=rhs.DFLT_TARGET_ORBITAL_ELEMENTS_UNIT;

  max_target_orbital_element_deviation = 
    rhs.max_target_orbital_element_deviation;

  target_orbital_elements_dev_input_unit =
    rhs.target_orbital_elements_dev_input_unit;

  DFLT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT =
    rhs.DFLT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT;

  target_body = rhs.target_body;
  prox_weights = rhs.prox_weights;

  max_dr = rhs.max_dr;
  max_delta_r_input_unit = rhs.max_delta_r_input_unit;
  DFLT_MAX_DELTA_R_UNIT = rhs.DFLT_MAX_DELTA_R_UNIT;
  INT_MAX_DELTA_R_UNIT = rhs.INT_MAX_DELTA_R_UNIT;
  max_dv = rhs.max_dv;
  max_delta_v_input_unit = rhs.max_delta_v_input_unit;
  DFLT_MAX_DELTA_V_UNIT = rhs.DFLT_MAX_DELTA_V_UNIT;
  INT_MAX_DELTA_V_UNIT = rhs.INT_MAX_DELTA_V_UNIT;

  min_solar_distance = rhs.min_solar_distance;
  min_solar_distance_unit = rhs.min_solar_distance_unit;
  DFLT_MIN_SOLAR_DISTANCE_UNIT = rhs.DFLT_MIN_SOLAR_DISTANCE_UNIT;
  INT_MIN_SOLAR_DISTANCE_UNIT = rhs.INT_MIN_SOLAR_DISTANCE_UNIT;

  max_error_rel = rhs.max_error_rel;
  max_error_abs = rhs.max_error_abs;
  max_clock_jump = rhs.max_clock_jump;
  resolution_factor = rhs.resolution_factor;
  
  max_itgr_step_size = rhs.max_itgr_step_size;
  min_itgr_step_size = rhs.min_itgr_step_size;

  // copy dynamic accuracy control (DAC) parameters
  dac_usage = rhs.dac_usage;
  dac_factor_max = rhs.dac_factor_max;
  dac_hc_start_size = rhs.dac_hc_start_size;
  dac_hc_final_size = rhs.dac_hc_final_size;

  // copy dynamic step size control (DSSC) parameters
  dssc_usage = rhs.dssc_usage;
  dssc_use_approach_control = rhs.dssc_use_approach_control;
  dssc_use_step_angle_control = rhs.dssc_use_step_angle_control;
  dssc_use_step_distance_control = rhs.dssc_use_step_distance_control;
  dssc_approach_step_size_factor = rhs.dssc_approach_step_size_factor;

  dssc_max_step_angle = rhs.dssc_max_step_angle;
  INT_DSSC_MAX_STEP_ANGLE_UNIT=rhs.INT_DSSC_MAX_STEP_ANGLE_UNIT;
  dssc_max_step_angle_input_unit=rhs.dssc_max_step_angle_input_unit;
  DFLT_DSSC_MAX_STEP_ANGLE_UNIT = rhs.DFLT_DSSC_MAX_STEP_ANGLE_UNIT;

  dssc_max_step_distance = rhs.dssc_max_step_distance;
  dssc_max_step_distance_input_unit = rhs.dssc_max_step_distance_input_unit;
  INT_DSSC_MAX_STEP_DISTANCE_UNIT = rhs.INT_DSSC_MAX_STEP_DISTANCE_UNIT;
  DFLT_DSSC_MAX_STEP_DISTANCE_UNIT = rhs.DFLT_DSSC_MAX_STEP_DISTANCE_UNIT;

  dssc_max_step_size = rhs.dssc_max_step_size;
  dssc_max_step_size_input_unit=rhs.dssc_max_step_size_input_unit;
  INT_DSSC_MAX_STEP_SIZE_UNIT=rhs.INT_DSSC_MAX_STEP_SIZE_UNIT;
  DFLT_DSSC_MAX_STEP_SIZE_UNIT =rhs.DFLT_DSSC_MAX_STEP_SIZE_UNIT;

  dssc_min_step_size = rhs.dssc_min_step_size;
  dssc_min_step_size_input_unit=rhs.dssc_min_step_size_input_unit;
  INT_DSSC_MIN_STEP_SIZE_UNIT=rhs.INT_DSSC_MIN_STEP_SIZE_UNIT;
  DFLT_DSSC_MIN_STEP_SIZE_UNIT=rhs.DFLT_DSSC_MIN_STEP_SIZE_UNIT;

  return (*this);
}

//------------------------------------------------------------------------------
void CInputSimulation::printDacSettings(const unsigned int w) const
{
  ostringstream ss;

  // maximum DAC accuracy factor
  ss.str("");
  ss << setw(w) << left << "SIM: max. DAC accuracy factor" << ": "
    << dac_factor_max;
  if (dac_factor_max==DFLT_DAC_FACTOR_MAX) ss << "(DEFAULT)";
  ss << endl;Output(ss.str());

  // DAC start hypercube size
  ss.str("");
  ss << setw(w) << left << "SIM: DAC hypercube start size" << ": "
    << dac_hc_start_size;
  if (dac_hc_start_size==DFLT_DAC_HC_START_SIZE) ss << "(DEFAULT)";
  ss << endl; Output(ss.str());

  // DAC final hypercube size
  ss.str("");
  ss << setw(w) << left << "SIM: DAC hypercube final size" << ": "
    << dac_hc_final_size;
  if (dac_hc_final_size==DFLT_DAC_HC_FINAL_SIZE) ss << "(DEFAULT)";
  ss << endl;Output(ss.str());
} // void CInputSimulation::printDacSettings(const unsigned int w) const


//------------------------------------------------------------------------------
void CInputSimulation::printDSSCSettings(const unsigned int w) const
{
  ostringstream ss;

  ss.str("");      // approach control
  ss << setw(w) << left << "SIM: APPROACH CONTROL"<<": "
    << (dssc_use_approach_control == true?"yes":"no") << endl;
  Output(ss.str());

  ss.str(""); // step distance control
  ss << setw(w) << left << "SIM: STEP DISTANCE CONTROL"<<": "
    << (dssc_use_step_distance_control == true?"yes":"no") << endl;
  Output(ss.str());

  ss.str(""); // step angle control
  ss << setw(w) << left << "SIM: STEP ANGLE CONTROL"<<": "
    << (dssc_use_step_angle_control == true?"yes":"no") << endl;
  Output(ss.str());

  if (dssc_use_approach_control == true)
  {
    ss.str("");   // approach step size factor
    ss << setw(w) << left << "SIM: APPROACH STEP SIZE FACTOR"<<": "
      << dssc_approach_step_size_factor << endl;
    Output(ss.str());
  }

  if (dssc_use_step_angle_control == true)
  {
    ss.str("");   // max orbit angle per step
    ss << setw(w) << left << "SIM: max. orbit angle per step"<<": ";
    ss << dssc_max_step_angle << INT_DSSC_MAX_STEP_ANGLE_UNIT;
    ss << " (";
    ss << " ";
    ss << INT_DSSC_MAX_STEP_ANGLE_UNIT->convertTo(dssc_max_step_angle
      ,dssc_max_step_angle_input_unit);
    ss << dssc_max_step_angle_input_unit;
    ss << ")";
    ss << endl;
    Output(ss.str());
  }

  if (dssc_use_step_distance_control == true)
  {
    ss.str("");   // max distance per step
    ss << setw(w) << left << "SIM: max. distance per step"<<": ";
    ss << dssc_max_step_distance << INT_DSSC_MAX_STEP_DISTANCE_UNIT;
    ss << " (";
    ss << " ";
    ss << INT_DSSC_MAX_STEP_DISTANCE_UNIT->convertTo(dssc_max_step_distance
      ,dssc_max_step_distance_input_unit);
    ss << dssc_max_step_distance_input_unit;
    ss << ")";
    ss << endl;
    Output(ss.str());
  }

  ss.str("");  // min. DSSC step size
  ss << setw(w) << left << "SIM: min. DSSC step size "<<": ";
  ss << dec << dssc_min_step_size << INT_DSSC_MIN_STEP_SIZE_UNIT;
  ss << " (";
  ss << " ";
  ss << INT_DSSC_MIN_STEP_SIZE_UNIT->convertTo(dssc_min_step_size
    ,dssc_min_step_size_input_unit);
  ss << dssc_min_step_size_input_unit;
  ss << ")";
  ss << endl;
  Output(ss.str());

  ss.str(""); // max. DSSC step size
  ss << setw(w) << left << "SIM: max. DSSC step size "<<": ";
  ss << dec << dssc_max_step_size << INT_DSSC_MAX_STEP_SIZE_UNIT;
  ss << " (";
  ss << " ";
  ss << INT_DSSC_MAX_STEP_SIZE_UNIT->convertTo(dssc_max_step_size
    ,dssc_max_step_size_input_unit);
  ss << dssc_max_step_size_input_unit;
  ss << ")";
  ss << endl;
  Output(ss.str());
} // void CInputSimulation::printDSSCSettings(const unsigned int w) const


//------------------------------------------------------------------------------
void CInputSimulation::printSettings(void) const
{
  ostringstream ss;
  const unsigned int w(37);

  // optimization goal
  ss.str("");
  ss << setw(w) << left<< "SIM: optimization goal"<<": " 
     << OPT_GOAL_DESC[opt_goal] << endl;
  Output(ss.str());

  // optimization goal
  ss.str("");
  ss << setw(w) << left<< "SIM: simulation type"<<": ";
  if (isHeliocentric()) ss << "heliocentric" << endl;
  else ss << "bodycentric" << endl;
  Output(ss.str());

  printInitStateSettings(w);

  // target state
  ss.str("");
  ss << setw(w) << left << "SIM: target state"<<": "
     << StateTypeDesc[target_state] << endl;
  Output(ss.str());

  printTargetOrbitalElementsSettings(w);

  // launch excess velocity vector parameters
  for (int i=0;i<N_VINF_PARAMS;++i)
  {
    ss.str("");
    ss << "SIM: " << setw(w-5) << left << INIT_VINF_PARAM_NAMES[i]
    << ": ";
    ss << launch_v_inf[i] << " " << INT_V_INF_UNIT[i];
    ss << " (" << INT_V_INF_UNIT[i]->convertTo(launch_v_inf[i],
      launch_v_inf_input_unit[i]);
    ss << " " << launch_v_inf_input_unit[i] << ")\n";
    Output(ss.str());
  }

  // resolution factor
  ss.str("");
  ss << setw(w) << left<< "SIM: max. input clock jump"<<": " 
     << max_clock_jump;
  if (max_clock_jump==DFLT_MAX_CLOCK_JUMP) ss<< " (default)";
  ss<<endl;Output(ss.str());

  // resolution factor
  ss.str("");
  ss << setw(w) << left<< "SIM: approx. resolution factor"<<": " 
     << resolution_factor;
  if (resolution_factor==DFLT_RESOLUTION_FACTOR) ss<< " (default)";
  ss<<endl;Output(ss.str());

  // allowed maximum final distance
  ss.str("");
  ss << setw(w) << left<< "SIM: maximum pos. error"<<": ";
  ss << max_dr << " " << INT_MAX_DELTA_R_UNIT;
  ss << " (";
  ss << INT_MAX_DELTA_R_UNIT->convertTo(max_dr, max_delta_r_input_unit);
  ss << " " << max_delta_r_input_unit;
  ss << ")"<<endl;
  Output(ss.str());

  // maximum velocity error
  ss.str("");
  ss << setw(w) << left<< "SIM: maximum vel. error"<<": ";
  ss << max_dv << " " << INT_MAX_DELTA_V_UNIT;
  ss << " (";
  ss << INT_MAX_DELTA_V_UNIT->convertTo(max_dv, max_delta_v_input_unit);
  ss << " " << max_delta_v_input_unit;
  ss << ")" << endl;
  Output(ss.str());

  // dynamic accuracy control
  ss.str("");
  ss << setw(w) << left << "SIM: dynamic accuracy control (DAC)" << ": "
    << (dac_usage==true?"yes":"no")<< endl;
  Output(ss.str());

  if (dac_usage==true) printDacSettings(w);

  // dynamic step size control
  ss.str("");
  ss << setw(w) << left << "SIM: dynamic step size control (DSSC)"<<": "
    << (dssc_usage==true?"yes":"no")<< endl;
  Output(ss.str());

  if (dssc_usage == true) printDSSCSettings(w);

  // maximum integration step size [s]
  ss.str("");
  ss << setw(w) << left<< "SIM: max. integration step size"<<": "
    << max_itgr_step_size <<" s";
  if (max_itgr_step_size==DFLT_MAX_ITGR_STEP_SIZE) ss << " (DEFAULT)";
  ss << endl;
  Output(ss.str());

  if (getTrajectoryType() == ACTPROBLEM) printACTParameterSettings(w);

}


//------------------------------------------------------------------------------
void CInputSimulation::printACTParameterSettings(
  const unsigned int w) const
{
  std::ostringstream ss("");
    // 1st ESA ACT competition parameter 1
  ss.str("");
  ss << setw(w) << left << "SIM: ACT distance factor" << ": ";
  ss << act_distance_factor << endl;
  Output(ss.str());

  // 1st ESA ACT competition parameter 2
  ss.str("");
  ss << setw(w) << left << "SIM: ACT exponent" << ": ";
  ss << act_exponent << endl;
  Output(ss.str());
} // 


//------------------------------------------------------------------------------
void CInputSimulation::printTargetOrbitalElementsSettings(
  const unsigned int w) const
{
  std::ostringstream ss("");
  size_t i = 0;
  while (i<use_target_orbital_elements.size())
  {
    if (use_target_orbital_elements[i] == true) break;
    ++i;
  }
  if (i<use_target_orbital_elements.size())
  {
    ss.str("");    
    for (unsigned int e = SEMIMAJORAXIS; e<N_KEPLERIAN_ELEMENTS;++e)
    {
      const KEPLERIAN_ELEMENT el= KEPLERIAN_ELEMENT(e);
      if (useTargetOrbitalElement(el) == true)
      {
        ss.str("");
        ss << scientific << setprecision(2);
        ss << setw(5) << left << "SIM: ";
        ss << setw(w>=5?w-5:5) << left 
          << TARGET_ELEMENT_PARAM_NAME[el] <<": ";

        ss << getTargetOrbitalElement(el);
        ss << " " << INT_TARGET_ORBITAL_ELEMENTS_UNIT[el];
        ss << " (";
        ss << INT_TARGET_ORBITAL_ELEMENTS_UNIT[el]->convertTo(
getTargetOrbitalElement(el), target_orbital_elements_input_unit[el]);
        ss << " " << target_orbital_elements_input_unit[el];

        ss << ") | ";
        ss << getTargetOrbitalElementMaxDeviation(el) 
          << INT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT[el];
        ss << " (";
        ss << INT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT[el]->convertTo(
getTargetOrbitalElementMaxDeviation(el)
, target_orbital_elements_dev_input_unit[el]);
        ss << " " << target_orbital_elements_dev_input_unit[el];
        ss << ")\n";
        Output(ss.str());   
      } // if (useTargetOrbitalElement(el) == true)
    } // for (unsigned int e = SEMIMAJORAXIS; e<N_KEPLERIAN_ELEMENTS;++e)
  } // if (i<use_target_orbital_elements.size())
} // void CInputSimulation::printTargetOrbitalElementsSettings(...) const

//----------------------------------------------------------------------------
void CInputSimulation::printInitStateSettings(
  const unsigned int w) const
{
  const std::string method_name = "void CInputSimulation::\
getInitialStateSettings(";
  const std::string warn_msg_01 = ": unknown initial state type!\n";
  std::ostringstream ss("");

  // initial state type
  ss << setw(w) << left << "SIM: initial state type" << ": ";
  ss << StateTypeDesc[getInitialStateType()];
  if (getInitialStateType()==DFLT_INIT_CENTRAL_BODY) ss << " (DEFAULT)";
  ss << endl;
  Output(ss.str());

  // initial central body
  ss.str("");
  ss << setw(w) << left << "SIM: initial central body" << ": ";
  ss << getInitCentralBodyName();
  if (getInitCentralBodyName()==DFLT_INIT_CENTRAL_BODY_NAME) ss << " (DEFAULT)";
  ss << endl;
  Output(ss.str());

  switch (getInitialStateType())
  {
  case BODY_STATE:
    printInitStateBodySettings(w);
    break;
  case FREE_ORBIT_STATE:
    printInitStateElementsSettings(w);
    break;
  case VECTOR_STATE:
    printInitStateVecSettings(w);
    break;
  default:
    Warning(method_name + warn_msg_01);
    break;
  } // switch (getInitialStateType())

} // void printInitialStateSettings(const unsigned int w) const


//----------------------------------------------------------------------------
void CInputSimulation::printInitStateBodySettings(const unsigned int w) const
{
  std::ostringstream ss("");

  // initial start body
  ss.str("");
  ss << setw(w) << left << "SIM: initial body" << ": ";
  ss << getInitialBodyName();
  if (getInitialBodyName()==DFLT_INIT_START_BODY_NAME) ss << " (DEFAULT)";
  ss << endl;
  Output(ss.str());
} // void printInitStateBodySettings(const unsigned int w) const


//----------------------------------------------------------------------------
void CInputSimulation::printInitStateElementsSettings(const unsigned int w)
const
{
  std::ostringstream ss("");
  unsigned int i;
  for (i=SEMIMAJORAXIS; i<N_KEPLERIAN_ELEMENTS;++i)
  {
    ss.str("");
    const KEPLERIAN_ELEMENT e = KEPLERIAN_ELEMENT(i);
    
    ss << "SIM: init. " << setw(w-11) << left << KEPLERIAN_ELEMENT_NAME[e];
    ss << ": ";
    ss << getInitialOrbitElement(e)<< " ";
    ss << *(DFLT_ELEMENTS_OF_INITIAL_ORBIT_UNIT[e]);
    if (getInitialOrbitElement(e) == DFLT_INIT_ORBITAL_ELEMENTS[e]) 
      ss << " (DEFAULT)";
    ss << endl;
    Output(ss.str());
  }
} // void printInitStateElementsSettings(const unsigned int w) const


//----------------------------------------------------------------------------
void CInputSimulation::printInitStateVecSettings(const unsigned int w)
const
{
  std::ostringstream ss("");

  ss << setw(w) << left << "SIM: initial pos. vec [AU]" << ": ";
  ss << getInitPosVec();
  ss << endl;Output(ss.str());

  ss.str("");
  ss << setw(w) << left << "SIM: initial vel. vec [AU/d]" << ": ";
  ss << getInitVelVec();
  ss << endl;Output(ss.str());

} // void printInitStateElementsSettings(const unsigned int w) const


//----------------------------------------------------------------------------
void CInputSimulation::setInitCentralBodyName(const std::string& name)
{
  const std::string method_name = "CInputSimulation::setInitCentralBodyName(\
                                  const std::string& name)";
  const std::string warn_msg_01 = ": couldn't find object '"+name+"'!\n";

  initial_central_body_name = name;
  if (getSolarSystem() != 0)
  {
    const BodyID id = 
      getSolarSystem()->getPlanetIndex(initial_central_body_name);
    if (id != 0) 
    {
      setInitCentralBody(id);
      if (getInitCentralBody() == SUN)
      {
        setHeliocentricSimulation(true);
      }
      else 
      {
        setHeliocentricSimulation(false);
      }
    }
    else
    {
      Warning(method_name + warn_msg_01);
    }
  }
} // void setInitCentralBodyName(const std::string& name)

//----------------------------------------------------------------------------
const std::string& CInputSimulation::getInitialBodyName(void) const
throw (std::logic_error&)
{
  const std::string method_name = "const std::string& CInputSimulation::\
getInitialBodyName(void) const";
  const std::string err_msg_01 = ": no valid CSolarSystem object!\n";

  if (getSolarSystem() != 0)
  {
    return getSolarSystem()->getBodyName(getInitialBody());
  }
  else 
  {
    throw std::logic_error(method_name + err_msg_01);
    return S_NOT_APPLICABLE;
  }
} // const std::string& getInitialBodyName(void) const