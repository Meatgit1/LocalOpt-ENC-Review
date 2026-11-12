/*!
\file CInputSimulation.h
\version 1.7
\brief CInputSimulation class interface.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 07.11.2006 (AO)

\b Purpose:
   This file provides the interface for the class CInputSimulation.
   CInputSimulation is an auxiliary class that contains the values of the InTrance simulation input file.

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

#ifndef CINPUTSIMULATION_H_INCLUDED
#define CINPUTSIMULATION_H_INCLUDED

// InTrance
#include <InTranceConst.h>

// SpaceLib
#include <Global.h>
#include <MLDE.h>
#include <SpaceLib.h>


/// default optimization criterion
const OptimizationGoalType DFLT_OPT_GOAL_TYPE = TRANSFER_TIME;

/// default max. final abs. relative distance to target state [km]
const double DFLT_MAX_ABS_REL_DISTANCE = 1.0e5;

/// default max. final abs. relative velocity w.r.t. target state [km/s]
const double DFLT_MAX_ABS_REL_VELOCITY = 1.0e-1;

/// default value of minimum start date [MJD]
const double DFLT_MIN_START_TIME = MJD_J2000;
/// default value of maximum start date [MJD]
const double DFLT_MAX_START_TIME = MJD_J2000;

// 1st ESA ACT Optimization Competition problem parameters
/// default value of param 1 of ACT fitness
const double DFLT_ACT_DISTANCE_FACTOR = 1.0e2;
/// default value of param 2 of ACT fitness
const double DFLT_ACT_EXPONENT = 3.0;

/// default time span of max. integration
const double DFLT_ITGR_INTVL = 0.0;

/// default value of integration interval reduction flag
const bool DFLT_ITGR_INTVL_RED = false;

/// default value of max. allowed control steps
const int DFLT_MAX_CTRL_STEPS = 0;

/// default value of minimum output points
const int DFLT_MIN_OUTPUT_POINTS = 0;

// by default the new trajectory should have the same step size as the old one
const double DFLT_MAX_CLOCK_JUMP(MathLib::TWOPI);
const double DFLT_RESOLUTION_FACTOR(1.0);
const double DFLT_MAX_ITGR_STEP_SIZE(5.0);
const double DFLT_MIN_ITGR_STEP_SIZE(0.0);

// By default, dynamic accuracy control (DAC) is switched off and all values
// are set to represent the behaviour without DAC.
const bool DFLT_DAC_USAGE = false; // no dynamic accuracy control
const double DFLT_DAC_FACTOR_MAX = 1.0; // no change = 1.0
const double DFLT_DAC_HC_START_SIZE = 1.0;
const double DFLT_DAC_HC_FINAL_SIZE = 1.0;

const bool DFLT_DSSC_USAGE = false; // do not use dynamic step size control
const bool DFLT_DSSC_USE_APPROACH_CONTROL = false;
const bool DFLT_DSSC_USE_STEP_ANGLE_CONTROL = false;
const bool DFLT_DSSC_USE_STEP_DISTANCE_CONTROL = false;
const double DFLT_DSSC_APPROACH_STEP_SIZE_FACTOR = 0.0;
const double DFLT_DSSC_MAX_STEP_ANGLE = (5.0*DEG2RAD); // [rad]
const double DFLT_DSSC_MAX_STEP_DISTANCE = 1.0e7*KM2AU; // [km]
const double DFLT_DSSC_MIN_STEP_SIZE = 1.0; // [d]
const double DFLT_DSSC_MAX_STEP_SIZE = 1.0e1; // [d]

const double DFLT_MIN_SOLAR_DISTANCE = 0.2;

const SpaceLib::StateType DFLT_INIT_STATE_TYPE = BODY_STATE;
const std::string DFLT_INIT_START_BODY_NAME = EARTH_NAME;
const BodyID DFLT_INIT_START_BODY = EARTH;
const std::string DFLT_INIT_CENTRAL_BODY_NAME = SUN_NAME;
const BodyID DFLT_INIT_CENTRAL_BODY = SUN;

enum INIT_VINF_PARAM
{
   VINF_MIN
  ,VINF_MAX
  ,VINF_AZI_MIN
  ,VINF_AZI_MAX
  ,VINF_ELE_MIN
  ,VINF_ELE_MAX
  ,N_VINF_PARAMS
};

/// parameter names of initial excess velocity vector description
const std::string INIT_VINF_PARAM_NAMES[] =
{
   "INITIAL_VINF_MIN"
  ,"INITIAL_VINF_MAX"
  ,"INITIAL_VINF_AZIMUTH_MIN"
  ,"INITIAL_VINF_AZIMUTH_MAX"
  ,"INITIAL_VINF_ELEV_MIN"
  ,"INITIAL_VINF_ELEV_MAX"
};

/// default target values of orbital elements target state
const MathLib::Vector DFLT_INIT_ORBITAL_ELEMENTS(
   1.0 // semimajor axis [AU]
  ,0.0 // eccentricity
  ,0.0 // inclination [rad]
  ,0.0 // longitude of ascending node [rad]
  ,0.0 // argument of pericenter [rad]
  ,0.0 // mean anomaly [rad]
);


/// default values for target state of orbital elements
const bool DFLT_USE_TARGET_ORBITAL_ELEMENTS[N_KEPLERIAN_ELEMENTS] =
{
   false // semimajor axis
  ,false // eccentricity
  ,false // inclination
  ,false // longitude of ascending node
  ,false // argument of pericenter
  ,false // mean anomaly
};

/// default target values of orbital elements target state
const double DFLT_TARGET_ORBITAL_ELEMENTS[N_KEPLERIAN_ELEMENTS] =
{
   1.0 // semimajor axis [AU]
  ,0.0 // eccentricity
  ,0.0 // inclination [rad]
  ,0.0 // longitude of ascending node [rad]
  ,0.0 // argument of pericenter [rad]
  ,0.0 // mean anomaly [rad
};

/// maximum allowed deviation from target orbital element value
const double DFLT_MAX_TARGET_ORBITAL_ELEMENT_DEVIATION[N_KEPLERIAN_ELEMENTS] =
{
   1.0e6*KM2AU   // semimajor axis [AU]
  ,5.0e-2        // eccentricity []
  ,5.0*DEG2RAD   // inclination [RAD]
  ,5.0*DEG2RAD   // longitude of ascending node [RAD]
  ,5.0*DEG2RAD   // argument of pericenter [RAD]
  ,5.0*DEG2RAD   // mean anomaly [RAD]
};


/// available differential equation system (DES) solver types
const MathLib::DESSolverType INTEGRATOR_TYPES[] = 
{
  MathLib::_RK4
 ,MathLib::_RKF54
 ,MathLib::_DOPRI8
 ,MathLib::_DESG
};

/// literal parameter values
const std::string INTEGRATOR_TYPE_NAMES[] =
{
  "RK4"
 ,"RKF54"
 ,"DOPRI8"
 ,"DESG"
};


/// literal descriptions of possible optimization problems to
/// solve with InTrance
const std::string OPTIMIZATION_PROBLEM_DESC[] =
{
  "BODYRENDEZVOUS",
  "BODYFLYBY",
  "RETROGRADEBODYIMPACT",
  "ACTBODYIMPACT",
  "ACTBODYORBITIMPACT",
  "SPIPROBLEM",
  "BODYORBITRENDEZVOUS",
  "BODYORBITFLYBY",
  "FREEORBITRENDEZVOUS",
  "FREEORBITFLYBY",
  "ORBITALELEMENTS",
  "GTOC2",
  "TEST1",
  "TEST2",
  "TEST3",
  "TEST4",
  "TEST5",
  "TEST6",
  "TEST7",
  "TEST8",
  "TEST9"
};

/// initial state types
const SpaceLib::StateType START_STATE_TYPES[] =
{
   BODY_STATE
  ,FREE_ORBIT_STATE
  ,VECTOR_STATE
};


/// literal parameter value
const std::string START_STATE_TYPES_PARAM_VALUE[] =
{
   "BODY"
  ,"ORBIT"
  ,"VECTOR"
};


/// target state types corresponding to desired optimization type
const SpaceLib::StateType TARGET_STATE_TYPES[] =
{
  BODY_STATE,                                   ///< "BODYRENDEZVOUS" - problem
  BODY_STATE,                                        ///< "BODYFLYBY" - problem
  BODY_STATE,                             ///< "RETROGRADEBODYIMPACT" - problem
  BODY_STATE,                                    ///< "ACTBODYIMPACT" - problem
  BODY_ORBIT_STATE,                         ///< "ACTBODYORBITIMPACT" - problem
  BODY_ORBIT_STATE,                                 ///< "SPIPROBLEM" - problem
  BODY_ORBIT_STATE,                        ///< "BODYORBITRENDEZVOUS" - problem
  BODY_ORBIT_STATE,                             ///< "BODYORBITFLYBY" - problem
  FREE_ORBIT_STATE,                        ///< "FREEORBITRENDEZVOUS" - problem
  FREE_ORBIT_STATE,                             ///< "FREEORBITFLYBY" - problem
  ORBITAL_ELEMENTS_STATE,                      ///< "ORBITALELEMENTS" - problem
  BODY_STATE,                                            ///< "GTOC2" - problem
  RESERVED_01,                                 ///< not yet used, test purposes
  RESERVED_01,                                 ///< not yet used, test purposes
  RESERVED_01,                                 ///< not yet used, test purposes
  RESERVED_01,                                 ///< not yet used, test purposes
  RESERVED_01,                                 ///< not yet used, test purposes
  RESERVED_01,                                 ///< not yet used, test purposes
  RESERVED_01,                                 ///< not yet used, test purposes
  RESERVED_01,                                 ///< not yet used, test purposes
  RESERVED_01                                  ///< not yet used, test purposes
};


/// trajectory types corresponding to desired optimization type
const TrajectoryType TARGET_TRAJECTORIES[] =
{
  RENDEZVOUS,                                    ///< "BODYRENDEZVOUS" - problem
  FLYBY,                                              ///< "BODYFLYBY" - problem
  RIMPACT,                                 ///< "RETROGRADEBODYIMPACT" - problem
  ACTPROBLEM,                                     ///< "ACTBODYIMPACT" - problem
  ACTPROBLEM,                                ///< "ACTBODYORBITIMPACT" - problem
  SPIPROBLEM,                                        ///< "SPIPROBLEM" - problem
  RENDEZVOUS,                               ///< "BODYORBITRENDEZVOUS" - problem
  FLYBY,                                         ///< "BODYORBITFLYBY" - problem
  RENDEZVOUS,                               ///< "FREEORBITRENDEZVOUS" - problem
  FLYBY,                                         ///< "FREEORBITFLYBY" - problem
  RENDEZVOUS,                                   ///< "ORBITALELEMENTS" - problem
  GTOC2TRAJ,                                                  ///< GTOC2 Problem
  TRAJ_TEST1,                                   ///< not yet used, test purposes
  TRAJ_TEST2,                                   ///< not yet used, test purposes
  TRAJ_TEST3,                                   ///< not yet used, test purposes
  TRAJ_TEST4,                                   ///< not yet used, test purposes
  TRAJ_TEST5,                                   ///< not yet used, test purposes
  TRAJ_TEST6,                                   ///< not yet used, test purposes
  TRAJ_TEST7,                                   ///< not yet used, test purposes
  TRAJ_TEST8,                                   ///< not yet used, test purposes
  TRAJ_TEST9                                    ///< not yet used, test purposes
};

/// name of initial orbital elements parameter
const std::string INIT_ELEMENTS_PARAM_NAME[] =
{
  "INITIAL_SEMIMAJOR_AXIS"
  ,"INITIAL_ECCENTRICITY"
  ,"INITIAL_INCLINATION"
  ,"INITIAL_LONGITUDE_OF_NODE"
  ,"INITIAL_ARGUMENT_OF_PERICENTER"
  ,"INITIAL_MEAN_ANOMALY"
};

/// name of parameter in input file
const std::string TARGET_ELEMENT_PARAM_NAME[] =
{
   "TARGET_SEMIMAJOR_AXIS"
  ,"TARGET_ECCENTRICITY"
  ,"TARGET_INCLINATION"
  ,"TARGET_LONGITUDE_OF_NODE"
  ,"TARGET_ARGUMENT_OF_PERICENTER"
  ,"TARGET_MEAN_ANOMALY"
};


/// name of max. deviation parameter value
const std::string TARGET_ELEMENT_MAX_DEVIATION_PARAM_NAME[] =
{
   "TARGET_MAX_ABS_DELTA_A"
  ,"TARGET_MAX_ABS_DELTA_E"
  ,"TARGET_MAX_ABS_DELTA_I"
  ,"TARGET_MAX_ABS_DELTA_LONGITUDE_OF_NODE"
  ,"TARGET_MAX_ABS_DELTA_ARGUMENT_OF_PERICENTER"
  ,"TARGET_MAX_ABS_DELTA_MEAN_ANOMALY"
};


/// default weight vector for target state = TARGET_ORBITAL_ELEMENTS problem
const double DFLT_WEIGHT_VECTOR[] = 
{
   1.0  // semimajor axis
  ,1.0  // eccentricity
  ,10.0 // inclination
  ,0.0  // longitude of nodes
  ,0.0  // argument of periapsis
  ,0.0  // mean anomaly
};

/// name of parameters associated with "target_state" = TARGET_ORBITAL_ELEMENTS
const std::string PROX_WEIGHT_PARAM_NAME[] =
{
  "WEIGHT_A",
  "WEIGHT_E",
  "WEIGHT_I",
  "WEIGHT_NODE",
  "WEIGHT_PERIAPSIS",
  "WEIGHT_MEAN_ANOMALY"
};

/// parameter description for proximity weights 
/// ("target_state" = TARGET_ORBITAL_ELEMENTS)
const std::string PROX_WEIGHT_PARAM_DESC[] =
{
  "semimajor axis change time weight"
 ,"eccentricity change time weight"
 ,"inclination axis change time weight"
 ,"longitude of node change time weight"
 ,"argument of periapsis change time weight"
 ,"mean anomaly change time weight"
};

/// parameter names of input position vector
const std::string INIT_STATE_POS_VECTOR_PARAM_NAMES[] =
{
   "INITIAL_POS_X"
  ,"INITIAL_POS_Y"
  ,"INITIAL_POS_Z"
};

/// parameter names of input velocity vector
const std::string INIT_STATE_VEL_VECTOR_PARAM_NAMES[] =
{
   "INITIAL_VEL_X"
  ,"INITIAL_VEL_Y"
  ,"INITIAL_VEL_Z"
};


/// InTrance simulation configuration input object
class CInputSimulation 
{

public:

  SpaceLib::CSolarSystem* m_solsys;                     ///< solar system object

  double act_distance_factor;         ///< special parameter of fitness function
                ///< for solving the problem of the 1st ACT optimization contest
  double act_exponent;                      ///< second param of 1st ACT problem

  double min_start_time;     ///< minimum sim. start time [modified julian date]
  double max_start_time;     ///< maximum sim. start time [modified julian date]

  double integration_interval;                  ///< integration interval [days]
  PhysUnit integration_interval_input_unit;  ///< input unit of itgr. interval
  PhysUnit ITGR_INTVL_REF_UNIT;  ///< expected standard unit of itgr. interval

  bool dynamic_integration_interval_reduction; 
                                  ///< dynamic reduction of integration interval

  int max_steps;                             ///< maximum number of output steps
  int min_output_points;                    ///< minimum number of output points

private:

  bool dac_usage;                       ///< dynamic accuracy control (DAC) flag
  double dac_factor_max;                   ///< maximum accuracy factor to apply
  double dac_hc_start_size;   ///< hypercube size, from where DAC should be used
  double dac_hc_final_size;           ///< hypercube sizeof final accuracy limit

  bool dssc_usage;              ///< dynamic step size control (DSSC) usage flag
  bool dssc_use_approach_control; 
                          ///< should target approach step size control be used?
  bool dssc_use_step_angle_control; 
                   ///< should orbit angular velocity step size control be used?
  bool dssc_use_step_distance_control;
                           ///< should total velocity step size control be used?
  double dssc_approach_step_size_factor; 
                         ///< factor to calculate step size from target distance
  double dssc_max_step_angle;    ///< DSSC max. orbit angle per control step
  PhysUnit dssc_max_step_angle_input_unit;
  ///< DSSC max. orbit angular velocity input unit
  PhysUnit INT_DSSC_MAX_STEP_ANGLE_UNIT;
  ///< internal unit of DSSC max. orbit angle per step
  PhysUnit DFLT_DSSC_MAX_STEP_ANGLE_UNIT;
  ///< DSSC max. orbit angular velocity reference (default) unit

                                   
  double dssc_max_step_distance;      ///< DSSC max. total distance per day [km]
  PhysUnit dssc_max_step_distance_input_unit;
  ///< DSSC max. distance input unit
  PhysUnit INT_DSSC_MAX_STEP_DISTANCE_UNIT;
  ///< internal unit of DSSC max. step distance parameter
  PhysUnit DFLT_DSSC_MAX_STEP_DISTANCE_UNIT;
  ///< DSSC max. distance reference (default) unit


  double dssc_min_step_size;         ///< DSSC minimum integration step size (d)
  PhysUnit dssc_min_step_size_input_unit; ///< [d]
  PhysUnit INT_DSSC_MIN_STEP_SIZE_UNIT;
  ///< internal unit of DSSC min. step size
  PhysUnit DFLT_DSSC_MIN_STEP_SIZE_UNIT;

  double dssc_max_step_size;         ///< DSSC maximum integration step size (d)
  PhysUnit dssc_max_step_size_input_unit; ///< [d]
  PhysUnit INT_DSSC_MAX_STEP_SIZE_UNIT;
  ///< internal unit of DSSC min. step size
  PhysUnit DFLT_DSSC_MAX_STEP_SIZE_UNIT;

public:

  /// pointer to used CSolarSystem object
  const CSolarSystem* getSolarSystem(void) const;

  /// pointer to used CSolarSystem object
  CSolarSystem* getSolarSystem(void);

  /// set a new CSolarSystem object
  void setSolarSystem(CSolarSystem* sol_sys);

  /// name of central body
  const std::string& getInitCentralBodyName(void) const;

  /// name of central body (non-const version)
  std::string& getInitCentralBodyName(void);

  /// get Id of central body
  BodyID getInitCentralBody(void) const;

  /// set new central body
  void setInitCentralBody(const BodyID);

  /// current distance factor for 1st ACT problem solution
  double getACTDistanceFactor(void) const;
  /// set new distance factor for 1st ACT problem solution
  void setACTDistanceFactor(const double factor);

  /// current exponent for 1st ACT problem solution
  double getACTExponent(void) const;
  /// set new exponent for 1st ACT problem solution
  void setACTExponent(const double exponent);

  /// return maximum allowed # of steering steps
  int getMaxStepsCount(void) const; 
  /// set max. allowed # of steering steps
  void setMaxStepsCount(const int count);

  /// min. # of final output solution steps
  int getMinOutputPoints(void) const;
  /// set new min. # of final output solution steps
  void setMinOutputPoints(const int n_min_points);

  /// return setting of whether dynamic accuracy should be used or not
  bool useDynamicAccuracyControl(void) const;
  /// enable/disable dynamic accuracy control
  void setDynamicAccuracyControl(const bool use_or_not);

  /// hypercube size from where to start dynamic accuracy control
  double getDACMaxHypercubeSize(void) const;

  /// set new max. hypercube size for accuracy control
  void setDACMaxHypercubeSize(const double val);

  /// min. hypercube size for dynamic accuracy control
  double getDACMinHypercubeSize(void) const;

  /// set new min. hypercube size for danymic accuracy control
  void setDACMinHypercubeSize(const double val);

  /// get max. accuracy factor
  double getDACMaxFactor(void) const;

  /// set new max. accuracy factor
  void setDACMaxFactor(const double val);


  /// return setting of whether the integration interval should be reduced to
  /// current best solutions flight time
  bool getDynamicIntegrationIntervalReduction(void) const;

  /// enable/disable dynamic integration interval reduction
  void setDynamicIntegrationIntervalReduction(const bool);

  /// return status of Dynamic Step Size Control
  bool getStatusDSSC(void) const;

  /// enable/disable Dynamic Step Size Control
  void setStatusDSSC(const bool status);

  /// enable Dynamic Step Size Control
  void enableDSSC(void);

  /// disable Dynamic Step Size Control
  void disableDSSC(void);

  /*! \brief return status of Step Angle (Orbital Angle) Control (DSSC)
   If set, this system detemines before every steering step the size of the
   next step dependent on the current orbital angular velocity. The current 
   orbit angular velocity is calculated from the rates of change of polar
   coordinate angles (orbital angular velocity = sqrt((dphi/dt)^2+(dtheta/dt)^2) ).
   Consequently, the step size (days, hours) is greater on lower orbital 
   angular velocities and lower on higher orbital angular velocities.
   \return currently set status
   \sa setStatusDSSCStepAngleControl, setDSSCMaxStepAngle, getDSSCMaxStepAngle
  */
  bool getStatusDSSCStepAngleControl(void) const;

  /// enable/disable DSSC step angle control
  void setStatusDSSCStepAngleControl(const bool status);

  ///  enable DSSC step angle control
  void enableDSSCStepAngleControl(void);

  ///  disable DSSC step angle control
  void disableDSSCStepAngleControl(void);

  /*! \brief set the max. allowed orbital step angle (DSSC)
    \param angle No steering step is allowed to exceed this value.
    \sa getDSSCMaxStepAngle,getStatusDSSCStepAngleControl
  */
  void setDSSCMaxStepAngle(const double angle);

  /*! \brief get the max. allowed orbital steering step angle (DSSC)
    \sa setDSSCMaxStepAngle,getStatusDSSCStepAngleControl
  */
  double getDSSCMaxStepAngle(void) const;

  /*! \brief return status of step size control (distance)
    If set, before everey steering step it is determined how big the next step
    should be. This is done using the current absolute spacecraft velocity.
    The new step size (time step) is calculated so that a given max.
    distance is not exceeded during control step.
   \return currently set status
   \sa setDSSCMaxStepDistance, getDSSCMaxStepDistance
  */
  bool getStatusDSSCStepDistanceControl(void) const;

  /// enable/disable steering step size control
  void setStatusDSSCStepDistanceControl(const bool status);

  /// enablesteering step size control
  void enableDSSCStepDistanceControl(void);

  /// disable steering step size control
  void disableDSSCStepDistanceControl(void);

  /*! \brief set steering step size control parameter of max. allowed step distance
    \param distance max. allowed distance of control step
    \ sa getDSSCMaxStepDistance,getStatusDSSCStepDistanceControl
  */
  void setDSSCMaxStepDistance(const double distance);

  /// return max. allowed distance per control step
  double getDSSCMaxStepDistance(void) const;

  /// return status (on/off) of target approach controlled steering step size
  bool getStatusDSSCApproachControl(void) const;

  /// enable/disable target approach controlled steering step size
  void setStatusDSSCApproachControl(const bool status);

  /// enable approach controlled steering step size (part of DSSC)
  void enableDSSCApproachControl(void);

  /// disable approach controlled steering step size (part of DSSC)
  void disableDSSCApproachControl(void);

  /*! \brief set DSSC approach control parameter
  The parameter given with factor will be multiplied with the current target 
  distance to determine the next steering step size (time step -> day, hour).
  So, the relative distance (AU) is converted by this factor to a time span.
  \param factor new DSSC approach control step size parameter
   \sa getDSSCApproachFactor,getStatusDSSCApproachControl
  */
  void setDSSCApproachFactor(const double factor);

  /// get DSSC approach control distance factor
  double getDSSCApproachFactor(void) const;

  /*! brief Minimum allowed DSSC step size (time step)
    \sa setDSSCMinStepSize
  */
  double getDSSCMinStepSize(void) const;

  /*! \brief Set minimum allowed DSSC step size (time step)
   \param step_size the new min step size
   \sa getDSSCMinStepSize
  */
  void setDSSCMinStepSize(const double step_size);

  /*! brief Maximum allowed DSSC step size (time step)
    \sa setDSSCMaxStepSize
  */
  double getDSSCMaxStepSize(void) const;

  /*! \brief Set maximum allowed DSSC step size (time step)
   \param step_size the new max step size
   \sa setDSSCMaxStepSize
  */
  void setDSSCMaxStepSize(const double step_size);


  /*! \brief initial condition of simulation
   \sa StateType
  */
  StateType getInitialStateType(void) const;

  /// set new initial simulation condition
  void setInitialStateType(const StateType state_type);

  /// initial body in case of a body initial state
  BodyID getInitialBody(void) const;

  /// return name of initial start body
  const std::string& getInitialBodyName(void) const throw (std::logic_error&);

  /*! \brief set Id of new initial body to use as start condition
   \param body_id Id of new body to start simulation from
  */
  void setInitialBody(const BodyID body_id);

  /*! \brief get a specific Keplerian element of initial orbit

    If a set of Keplerian elements, which mustn't be complete, is used
    as initial condition for orbit optimization, a particular element of the
    currently used set of orbit parameters can be obtained.

    \param element the desired Keplerian element (see also KEPLERIAN_ELEMENT)
    \return current value of desired Keplerian element
    \sa setInitialOrbitElement
  */
  double getInitialOrbitElement(const SpaceLib::KEPLERIAN_ELEMENT element) const;

  /*! \brief set a specific Keplerian element of initial orbit element set
  * 
  * \param element the element that is subject to change
  * \param value value to set for specified initial condition set's element
  * \sa getInitialOrbitElement
  */
  void setInitialOrbitElement(const SpaceLib::KEPLERIAN_ELEMENT element,
    const double value);

  /*! \brief status of specific Keplerian element usage for target orbit
   *
   * \param e the Keplerian element of which the usage status is needed
   * \return current usage setting of specified element
   * \sa getTargetOrbitalElement
   */
  bool useTargetOrbitalElement(const KEPLERIAN_ELEMENT e) const;

  /*! \brief Keplerian element of target orbit
   *
   * In case of an orbit used as target state, of which not all elements must
   * be specified, every element can be read by this method.
   *
   * \param e the desired Keplerian element of target orbit
   * \return value of desired target orbit element
   * \sa useTargetOrbitalElement
   */
  double getTargetOrbitalElement(const KEPLERIAN_ELEMENT e) const;

  /*! \brief maximum allowed deviation of current element from specified one
   *
   * If the target state is an orbit specified by a, possibly incomplete, set of 
   * Keplerian elements then a maximum allowed deviation of current orbits elements
   * with respect to the specified ones must be given. Otherwise it is not possible
   * to determine whether the target state (the specified orbit) is reached with
   * sufficient accuracy or not. This max. deviation may be specified for every
   * element of the target orbit. If not, default values will be used.
   *
   * \param e Keplerian element of which the max. allowed deviation is needed
   * \return current set maximum allowed element deviation
   * \sa setTargetOrbitalElementMaxDeviation, useTargetOrbitalElement
   * \sa DFLT_MAX_TARGET_ORBITAL_ELEMENT_DEVIATION
   */
  double getTargetOrbitalElementMaxDeviation(
    const KEPLERIAN_ELEMENT e) const;

  /*! \brief set maximum allowed deviation fo current element from specified one
   *
   * This method can, and should, be used to alter the built-in allowed max.
   * deviation of orbit elements of date to the specified orbit elements.
   *
   * \param e Keplerian element of which the max. allowed deviation is subject to
   change
   * \param val new max. allowed deviation of specified element
   * \sa setTargetOrbitalElementMaxDeviation, useTargetOrbitalElement
   * \sa DFLT_MAX_TARGET_ORBITAL_ELEMENT_DEVIATION
   */
  void setTargetOrbitalElementMaxDeviation(
     const KEPLERIAN_ELEMENT e
    ,const double val);

  /*! \brief Id of current target body
   * 
   * If the current optimization's target is a celestial body (e.g. in case of 
   * body rendezvous or body fly-by type mission), this particular body is
   * accessible in CSolarSystem objects by the Id returned by this method.
   *
   * \return Id of body currently set as target body
   * \sa setTargetBodyId
   */
  BodyID getTargetBodyId(void) const;

  /*! \brief set new target body by providing it's Id as parameter
   *
   * If the current optimization's target is a celestial body (e.g. in case of 
   * body rendezvous or body fly-by type mission), this particular body is
   * accessible in CSolarSystem objects by the Id provided to this method.
   * So, basically, you change the target of the current optimization.
   *
   * \param target_id Id of new target body
   * \sa getTargetBodyId
   */
  void setTargetBodyId(const BodyID target_id);

  /*! \brief current trajectory type
   *
   * \return currently set type of desired trajectory
   * \sa setTrajectoryType
   */
  TrajectoryType getTrajectoryType(void) const;

  /*! \brief set new type of trajectory
   *
   * \param traj_type type of trajectory to use from now on
   * \sa getTrajectoryType
   */
  void setTrajectoryType(const TrajectoryType traj_type);

  /*! \brief get type of optimization goal
   * \return currently set optimization goal type
   * \sa setOptGoalType
   */
  OptimizationGoalType getOptGoalType(void) const;

  /*! \brief set new optimization goal type
   * \param nog new type of optimization goal
   * \sa getOptGoalType
   */
  void setOptGoalType(const OptimizationGoalType nog);

  /*! \brief current target state type
   *
   * \return currently set target state type */
  StateType getTargetStateType(void) const;

  /*! \brief set new target state type
   * 
   * \param new_type new type of target state to set
   * \sa getTargetStateType */
  void setTargetStateType(const StateType new_type);

//> subject to possibly deletion
  double getTargetOrbitalElementProxWeight(
     const KEPLERIAN_ELEMENT e) const;

  void setTargetOrbitalElementProxWeight(
     const KEPLERIAN_ELEMENT e
    ,const double val);

  const std::valarray<double>& getTargetOrbitalElementProxWeights(void) const;
  std::valarray<double>& getTargetOrbitalElementProxWeights(void);
//> subject to possibly deletion

  /*! \brief get current integrator type (for integration of equations of motion)
   * 
   * \return currently set integrator object class type
   * \sa setIntegratorType*/
  DESSolverType getIntegratorType(void) const;

  /*! \brief set new integrator type (for integration of equations of motion)
   *
   * \param integrator_type integrator type to use from now on
   * \sa getIntegratorType */
  void setIntegratorType(const DESSolverType integrator_type);

  /*! \brief currently allowed maximum relative integration error
   *
   * The used integrator objects perform error checks after every internal
   * substep to assure that given max. values of absolute and relative errors
   * are not reached. Those error boundaries may be provided by the user to
   * adapt optimization speed because, naturally, higher error error boundaries
   * result in fewer integration steps and, therefore, higher integration speed.
   *
   * \sa setMaxRelIntegrationError*/
  double getMaxRelIntegrationError(void) const;

  /*! \brief set new allowed max. integration error (relative)
   *
   * For an explanation see getMaxRelIntegrationError.
   *
   * \param rel_error new value of max. error (rel.)
   * \sa getMaxRelIntegrationError */
  void setMaxRelIntegrationError(const double rel_error);

  /*! \brief currently allowed maximum absolute integration error
   *
   * For an explanation see getMaxRelIntegrationError.
   *
   * \sa setMaxAbsIntegrationError*/
  double getMaxAbsIntegrationError(void) const;

  /*! \brief set new allowed max. integration error (absolute)
   *
   * For an explanation see getMaxRelIntegrationError.
   *
   * \param abs_error newvalue of max. error (abs.)
   * \sa getMaxAbsIntegrationError*/
  void setMaxAbsIntegrationError(const double abs_error);

  /*! \brief minimum allowed distance from the Sun
   *
   * Due to the extrem conditions near the Sun there is a mininimum distance
   * that still allows to operate a spacecraft. This spacecraft/propulsion system-
   * dependend boundary value can be specified by the user to adapt the simulation
   * to the respective problem.
   *
   * \return currently set minimum value of allowed solar distance
   * \sa setMinSolarDistance*/
  double getMinSolarDistance(void) const;

  /*! \brief set minimum allowed distance from the Sun
   *
   * For an explanation see getMinSolarDistance.
   *
   * \param min_distance new value for min. distance from the Sun
   * \sa getMinSolarDistance*/
  void setMinSolarDistance(const double min_distance);

  /*! \brief Maximum final value of relative distance to target
   *
   * During the evaluation of possible problem solution by integrating
   * a trajectory, there often has a check to be performed whether certain
   * boundary conditions have already been met. For instance, if the task is
   * a rendezvous to a target body the final relative distance and velocity 
   * between this body and the spacecraft should be object to minimization.
   * How small the final differences of distance and/or velocity can be,
   * is determined by parameters which can be set and obtained by public methods.
   * After these boundary conditions are met, the primary optimization goal, 
   * e.g. the transfer time, is subject to optimization. This way it is assured to
   * get valid problem solutions.
   *
   * \return currently set value of max. rel. distance to target
   * \sa setMaxFinalRelDistance*/
  double getMaxFinalRelDistance(void) const;

  /*! \brief set new allowed maximum relative distance to target
   *
   * For an explanation see getMaxFinalRelDistance.
   *
   * \param max_final_rel_dist new allowed max. final distance
   * \sa getMaxFinalRelDistance*/
  void setMaxFinalRelDistance(const double max_final_rel_dist);

  /*! \brief Allowed maximum final relative velocity w.r.t target
   *
   * For an explanation see getMaxFinalRelDistance.
   *
   * \return current setting
   * \sa setMaxFinalRelVelocity*/
  double getMaxFinalRelVelocity(void) const;

  /*! \brief set new allowed maximum relative velocity w.r.t. target
   *
   * For an explanation see getMaxFinalRelDistance.
   *
   * \param max_final_rel_vel new allowed max. final rel. velocity
   * \sa getMaxFinalRelVelocity*/
  void setMaxFinalRelVelocity(const double max_final_rel_vel);

  const MathLib::Vector& getInitPosVec(void) const;
  void setInitPosVec(const MathLib::Vector& pos);
  void setInitPosVec(double p1, double p2, double p3);

  const MathLib::Vector& getInitVelVec(void) const;
  void setInitVelVec(const MathLib::Vector& vel);
  void setInitVelVec(double v1, double v2, double v3);

private:
  StateType initial_state;                ///< initial state type: body or orbit
  BodyID initial_body;                       ///< index of initial body, if body

  std::vector< std::vector<CBody*> > gtoc2_bodies;       ///< GTOC2-problem data

  MathLib::Vector elements_of_initial_orbit;
                                           ///< initial orbit elements (if used)
  std::vector<SpaceLib::CPhysUnit*> elements_of_initial_orbit_input_unit;
  std::vector<SpaceLib::CPhysUnit*> DFLT_ELEMENTS_OF_INITIAL_ORBIT_UNIT;

  bool body_centric; 
  ///< false if heliocentric simulation, true if planetocentric

  MathLib::Vector initial_position;             ///< s/c initial position vector
  std::vector<SpaceLib::CPhysUnit*> initial_position_input_units;
  std::vector<SpaceLib::CPhysUnit*> DFLT_INITIAL_POSITION_UNITS;

  MathLib::Vector initial_velocity;             ///< s/c initial velocity vector
  std::vector<SpaceLib::CPhysUnit*> initial_velocity_input_units;
  std::vector<SpaceLib::CPhysUnit*> DFLT_INITIAL_VELOCITY_UNITS;

  std::string initial_central_body_name; ///< name of initial orbit central body
  BodyID initial_central_body;                   ///< Id of initial central body

  StateType target_state; ///< target state type: body, body orbit or free orbit
  BodyID target_body;                                  ///< index of target body

  std::valarray<double> launch_v_inf;    ///< initial excess velocity parameters
  std::vector<SpaceLib::CPhysUnit*> launch_v_inf_input_unit;
  std::vector<SpaceLib::CPhysUnit*> DFLT_V_INF_INPUT_UNIT;
  std::vector<SpaceLib::CPhysUnit*> INT_V_INF_UNIT;

  TrajectoryType trajectory;   ///< trajectory type: rendezvous, flyby or impact
  OptimizationGoalType opt_goal;
                        ///< optimization goal: transfer time or propellant mass

  /// array of switches whether the specific orbital element should be used
  std::valarray<bool> use_target_orbital_elements;

  /// orbital element values of target orbit (if used)
  std::valarray<double> target_orbital_elements;
  std::vector<SpaceLib::PhysUnit> target_orbital_elements_input_unit;
  std::vector<SpaceLib::PhysUnit> DFLT_TARGET_ORBITAL_ELEMENTS_UNIT;
  std::vector<SpaceLib::PhysUnit> INT_TARGET_ORBITAL_ELEMENTS_UNIT;

  /// weights for computing proximity value in case of "target_state"=TARGET_ORBITAL_ELEMENTS
  std::valarray<double>  prox_weights;

  /// max. orbital element deviation
  std::valarray<double> max_target_orbital_element_deviation;
  std::vector<SpaceLib::PhysUnit> target_orbital_elements_dev_input_unit;
  std::vector<SpaceLib::PhysUnit> DFLT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT;
  std::vector<SpaceLib::PhysUnit> INT_TARGET_ORBITAL_ELEMENTS_DEV_UNIT;

  /// minimum sun-spacecraft distance [AU]
  double min_solar_distance;
  CPhysUnit* min_solar_distance_unit;
  CPhysUnit* DFLT_MIN_SOLAR_DISTANCE_UNIT;
  CPhysUnit* INT_MIN_SOLAR_DISTANCE_UNIT;


  DESSolverType      des_solver;    ///< DES solver: RK4, RKF54, DOPRI8, or DESG
  double             max_error_rel; 
    ///< maximum relative error for DES solver (only for RKF54, DOPRI8 and DESG)
  double             max_error_abs; 
    ///< maximum absolute error for DES solver (only for RKF54, DOPRI8 and DESG)

public:
  int number_of_ga_bodies;        ///< number of potential gravity assist bodies
  std::vector<bool>  ga_body;             
                          ///< potential gravity assist body usage array pointer
  double max_dr;                ///< maximum final distance to target state [km]
  PhysUnit max_delta_r_input_unit;
  PhysUnit DFLT_MAX_DELTA_R_UNIT;
  PhysUnit INT_MAX_DELTA_R_UNIT;
  double max_dv;       ///< maximum final relative velocity w.r.t. target [km/s]
  PhysUnit max_delta_v_input_unit;
  PhysUnit DFLT_MAX_DELTA_V_UNIT;
  PhysUnit INT_MAX_DELTA_V_UNIT;

  double impact_time;                 ///< asteroid's impact time on target body
  double impact_efficiency;
                  ///< measure for transfered linear momentum of impacting body
  double max_clock_jump;  
                ///< max jump in clock angle due to angle limitation to [-PI,PI]
  double resolution_factor;      ///< factor for calculating new number of steps

  double max_itgr_step_size;                  ///< maximum integration step size
  double min_itgr_step_size;                  ///< minimum integration step size

private:
  /// initialize member fields
  void init(void);

  /// initialize DSSC members
  void initDSSC(void);

  /// initialize GTOC 2 params
  void initGTOC2(void);

  /// initialize member fields for orbital elements type target
  void initTargetOrbitalElementsFields(void);

  /// initialize initial orbit state parameters
  void initInitialOrbitStateFields(void);

  /// initialize initial vector state parameters
  void initInitialVectorStateFields(void);

  /// set new name of central body
  void setInitCentralBodyName(const std::string& name);

  /// set internal flag
  void setHeliocentricSimulation(bool helio_centric);

  /// read setting whether the integration interval shall be reduced depending
  /// on the best solution found so far or not
  void readDynIntegrationReductionSetting(Global::ParamList&)
    throw (std::invalid_argument);

  /// read and check input file type
  void readFileIdentifier(Global::ParamList&) throw (std::invalid_argument&);

  /// basically read name of initial state body and insert it into
  /// internal solar system if necessary
  void readInitStateBodyParams(Global::ParamList&) 
    throw (std::invalid_argument&);

  /// read orbit parameters of initial state orbit
  void readInitStateOrbitParams(Global::ParamList&)
    throw (std::invalid_argument&);

  /// read parameters of vector initial state
  void readInitStateVectorParams(Global::ParamList&)
    throw (std::invalid_argument&);

  /// read initial state parameters depending on previously read initial
  /// state type
  void readInitialStateParameters(Global::ParamList&)
    throw (std::invalid_argument&);

  /// read initial (at launch) state type
  void readInitialStateType(Global::ParamList&) throw (std::invalid_argument);

  /// read longest integration interval in days
  void readIntegrationInterval(Global::ParamList&)
    throw (std::invalid_argument&);

  /// read maximum absolute integration error
  void readIntegrationMaxAbsError(Global::ParamList&)
    throw (std::invalid_argument&);

  /// read maximum relative integration error
  void readIntegrationMaxRelError(Global::ParamList&)
    throw (std::invalid_argument&);

  /// read type of integrator
  void readIntegrator(Global::ParamList&) throw (std::invalid_argument);

  /// read longest integration interval in days
  void readMaxIntegrationSteps(Global::ParamList&)
    throw (std::invalid_argument&);

  /// read simulation latest start time
  void readMaxStartTime(Global::ParamList&) throw (std::invalid_argument&);

  /// read minimum output points count (file output)
  void readMinOutputPoints(Global::ParamList&) throw (std::invalid_argument);

  /// read minimum allowed solar distance
  void readMinSolarDistance(Global::ParamList&) throw (std::invalid_argument);

  /// read simulation earliest start time
  void readMinStartTime(Global::ParamList&) throw (std::invalid_argument&);

  /// read and parse accuracy boundaries
  void readAccuracyBoundaries(Global::ParamList&)
    throw (std::invalid_argument&);

  /// read final relative distance parameter
  void readMaxFinalRelDistParam(Global::ParamList& p);

  /// read final relative velocity parameter
  void readMaxFinalRelVelParam(Global::ParamList& p);

  /// read special paramters for 1st ESA ACT competition
  void readActProblemParameters(Global::ParamList&)
    throw (std::invalid_argument&);

  /// read special parameters of 2nd Global Trajectory Optimization Competition
  void readGTOCOMP2Parameters(Global::ParamList&)
    throw (std::invalid_argument&);

  /// read angle approximation parameters
  void readApproximationParameters(Global::ParamList&);

  /// read an parse Dynamic Accuracy Control parameters
  void readDACParameters(Global::ParamList&) throw (std::invalid_argument&);

  /// read and parse deflection parameters
  void readDeflectionParameters(Global::ParamList&) 
    throw (std::invalid_argument&);

  /// read an parse Dynamic Step Size Control parameters
  void readDSSCParameters(Global::ParamList&) throw (std::invalid_argument&);

  /// adapt some DSSC parameters according to set boundary constraints
  void modifyDSSCParameters(void);

  /// read and parse possible gravity assist bodies
  void readGravitiyAssistBodies(Global::ParamList&) 
    throw (std::invalid_argument&);

  /// read parameters of possibly given hyperbolic excess velocity
  void readHyperbolicExcessVelocityParameters(Global::ParamList&)
    throw (std::invalid_argument&);

  /// read and parse optimization goal
  void readOptimizationGoal(Global::ParamList&)
    throw (std::invalid_argument&);

  /// read and parse hopefully known target body or orbit
  void readTargetBodySettings(Global::ParamList&)
    throw (std::invalid_argument&);

  /// read and parse hopefully known target body and use its orbit as target
  void readTargetBodyOrbitSettings(Global::ParamList& p)
  throw (std::invalid_argument);

  /// read and parse parameters of free orbit target type
  void readTargetFreeOrbitSettings(Global::ParamList&)
    throw (std::invalid_argument&);

  /*! \brief read simulation parameter
   * Tries to read a parameter with possibly given unit and stores read parameter
   * value in given internal units, so that a conversion might be needed. If
   * no unit is found the given default unit is used as input unit.
   * \param p list of all parameters
   * \param name name of desired parameter
   * \param val reference to variable where the parameter value will be stored
   * \param input_unit reference to input unit (will be changed)
   * \param default_unit default unit object
   * \param internal_unit internal unit object
   */
  void CInputSimulation::readParameter(Global::ParamList& p
 ,std::string name, double& val, PhysUnit& input_unit, PhysUnit& default_unit
 ,PhysUnit& internal_unit) 
 throw (std::invalid_argument&, std::logic_error&);

  /// read and parse parameters for target state type orbital elements
  void readTargetOrbitalElementsSettings(Global::ParamList&)
    throw (std::invalid_argument&, std::logic_error&);

  /// read and parse desired target state type
  void readTargetType(Global::ParamList&) throw (std::invalid_argument);

  /// read and parse target state type settings
  void readTargetTypeSettings(Global::ParamList&)
    throw (std::invalid_argument&);

  /// print initial vector state settings
  void printInitStateBodySettings(const unsigned int w) const;

  /// print initial vector state settings
  void printInitStateVecSettings(const unsigned int w) const;

  /// print initial orbital elements state settings
  void printInitStateElementsSettings(const unsigned int w) const;

public:

  CInputSimulation(void) {init();}; ///< standard constructor

  /// Constructs CInputSimulation object from parameter list and 
  /// pointer to solar system object
  CInputSimulation(ParamList& parameters, SpaceLib::CSolarSystem* solsys);

  /// heliocentric simulation ?
  bool isHeliocentric(void) const;

  /// planetocentric simulation ?
  bool isBodycentric(void) const;

  
  /*! \brief get a asteroid body of GTOC2 - problem
   * \param group asteroid group
   * \param no # of asteroid in group
   * \return CBody object of desired asteroid if existent, 0 otherwise
  */
  CBody* getGTOC2Body(const int group, const unsigned int no) const;

  /// assign operator
  const CInputSimulation& operator=(const CInputSimulation& rhs);

  /// print dynamic accuracy control (DAC) settings
  void printDacSettings(const unsigned int w) const;

  /// print dynamic dynamic step size control (DSSC) settings
  void printDSSCSettings(const unsigned int w) const;

  void printSettings(void) const;             ///< prints stored actual settings

  /// print initial state settings
  void printInitStateSettings(const unsigned int w) const;

  /// print special ACT competition parameters
  void printACTParameterSettings(const unsigned int w) const;

  /// print settings for target state type "ORBITAL ELEMENTS"
  void printTargetOrbitalElementsSettings(const unsigned int w) const;

  ~CInputSimulation() {};                                        ///< destructor

  inline void reset(void) {init();};       ///< initialize with default settings

  /// get minimum allowed start/launch date as Modified Julian Date [MJD]
  double getMinStartTime(void) const {return min_start_time;};
  /// set minimum allowed start/launch date as Modified Julian Date [MJD]
  void setMinStartTime(const double val) {min_start_time=val;};

  /// get maximum allowed start/launch date as Modified Julian Date [MJD]
  double getMaxStartTime(void) const {return max_start_time;};
  /// set maximum allowed start/launch date as Modified Julian Date [MJD]
  void setMaxStartTime(const double val) {max_start_time=val;};

  /// get maximum integration interval in days [d]
  double getIntegrationInterval(void) const {return integration_interval;};
  /// set maximum integration interval in days [d]
  void setIntegrationInterval(const double val) {integration_interval=val;};

  /*! \brief Minimum launch (optimization start) absolute s/c velocity.
   *
   * \return Setting of minimum s/c velocity at the beginning of flight.
   * \sa setMinLaunchVelocity, getMaxLaunchVelocity */
  double getMinLaunchVelocity(void) const;

  /*! \brief Set new value of minimum launch (optimization start) velocity.
   *
   * \param vel new setting of minimum s/c launch velocity
   * \sa getMinLaunchVelocity, setMaxLaunchVelocity*/
  void setMinLaunchVelocity(const double vel);

  /*! \brief Maximum launch (optimization start) absolute s/c velocity.
   *
   * \return Setting of maximum s/c velocity at the beginning of flight.
   * \sa setMaxLaunchVelocity, getMinLaunchVelocity */
  double getMaxLaunchVelocity(void) const;

  /*! \brief Set new value of maximum launch (optimization start) velocity.
   *
   * \param vel new setting of maximum s/c launch velocity
   * \sa getMaxLaunchVelocity, setMinLaunchVelocity*/
  void setMaxLaunchVelocity(const double vel);

  /*! \brief Minimum launch (optimization start) azimuth angle of s/c.
   *
   * \return Current setting of minimum launch azimuth angle.
   * \sa setMinLaunchAzimuthAngle, getMaxLaunchAzimuthAngle*/
  double getMinLaunchAzimuthAngle(void) const;

  /*! \brief Set new value of minimum launch (optimization start) azimuth angle.
   *
   * \param angle new setting of minimum s/c launch azimuth angle
   * \sa getMinLaunchAzimuthAngle, setMaxLaunchAzimuthAngle*/
  void setMinLaunchAzimuthAngle(const double angle);

  /*! \brief Maximum launch (optimization start) azimuth angle of s/c.
   *
   * \return Current setting of maximum launch azimuth angle.
   * \sa setMaxLaunchAzimuthAngle, getMinLaunchAzimuthAngle*/
  double getMaxLaunchAzimuthAngle(void) const;

  /*! \brief Set new value of maximum launch (optimization start) azimuth angle.
   *
   * \param angle new setting of maximum s/c launch azimuth angle
   * \sa getMaxLaunchAzimuthAngle, setMinLaunchAzimuthAngle*/
  void setMaxLaunchAzimuthAngle(const double angle);

  /*! \brief Minimum launch (optimization start) elevation angle of s/c.
   *
   * \return Current setting of minimum launch elevation angle.
   * \sa setMinLaunchElevationAngle, getMaxLaunchElevationngle*/
  double getMinLaunchElevationAngle(void) const;

  /*! \brief Set new value of minimum launch (optimization start) elevation angle.
   *
   * \param angle new setting of minimum s/c launch elevation angle
   * \sa getMinLaunchElevationAngle, setMaxLaunchElevationAngle*/
  void setMinLaunchElevationAngle(const double angle);

  /*! \brief Maximum launch (optimization start) elevation angle of s/c.
   *
   * \return Current setting of maximum launch elevation angle.
   * \sa setMaxLaunchElevationAngle, getMinLaunchElevationngle*/
  double getMaxLaunchElevationAngle(void) const;

  /*! \brief Set new value of maximum launch (optimization start) elevation angle.
   *
   * \param angle new setting of maximum s/c launch elevation angle
   * \sa getMaxLaunchElevationAngle, setMinLaunchElevationAngle*/
  void setMaxLaunchElevationAngle(const double angle);

}; // end class CInputSimulation

//------------------------------------------------------------------------------
inline double CInputSimulation::getACTDistanceFactor(void) const
{
  return act_distance_factor;
} // double getACTDistanceFactor(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setACTDistanceFactor(const double factor)
{
  act_distance_factor = factor;
} // void setACTDistanceFactor(const double factor)


//------------------------------------------------------------------------------
inline double CInputSimulation::getACTExponent(void) const
{
  return act_exponent;
} // double getACTExponent(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setACTExponent(const double exponent)
{
  act_exponent = exponent;
} // void setACTExponent(const double exponent)

//------------------------------------------------------------------------------
inline bool CInputSimulation::getDynamicIntegrationIntervalReduction(void) const
{
  return dynamic_integration_interval_reduction;
} // bool CInputSiumaltion::getDynamicIntegrationIntervalReduction(void) const


//------------------------------------------------------------------------------
inline 
void CInputSimulation::setDynamicIntegrationIntervalReduction(const bool value)
{
  dynamic_integration_interval_reduction = value;
} // void setDynamicIntegrationIntervalReduction(const bool value)


//------------------------------------------------------------------------------
inline double CInputSimulation::getMinLaunchVelocity(void) const
{
  return launch_v_inf[VINF_MIN];
} // double CInputSimulation::getMinStartVelocity(void) const


//------------------------------------------------------------------------------
inline void CInputSimulation::setMinLaunchVelocity(const double vel)
{
  launch_v_inf[VINF_MIN] = vel;
} // void CInputSimulation::setMinLaunchVelocity(const double vel)


//------------------------------------------------------------------------------
inline double CInputSimulation::getMaxLaunchVelocity(void) const
{
  return launch_v_inf[VINF_MAX];
} // double CInputSimulation::getMaxStartVelocity(void) const


//------------------------------------------------------------------------------
inline void CInputSimulation::setMaxLaunchVelocity(const double vel)
{
  launch_v_inf[VINF_MAX] = vel;
} // void CInputSimulation::setMaxLaunchVelocity(const double vel)


//------------------------------------------------------------------------------
inline double CInputSimulation::getMinLaunchAzimuthAngle(void) const
{
  return launch_v_inf[VINF_AZI_MIN];
} // double CInputSimulation::getMinLaunchAzimuthAngle(void) const


//------------------------------------------------------------------------------
inline void CInputSimulation::setMinLaunchAzimuthAngle(const double angle)
{
  launch_v_inf[VINF_AZI_MIN] = angle;
} // void CInputSimulation::setMinLaunchAzimuthAngle(const double angle)


//------------------------------------------------------------------------------
inline double CInputSimulation::getMaxLaunchAzimuthAngle(void) const
{
  return launch_v_inf[VINF_AZI_MAX];
} // double CInputSimulation::getMaxLaunchAzimuthAngle(void) const


//------------------------------------------------------------------------------
inline void CInputSimulation::setMaxLaunchAzimuthAngle(const double angle)
{
  launch_v_inf[VINF_AZI_MAX] = angle;
} // void CInputSimulation::setMaxLaunchAzimuthAngle(const double angle)


//------------------------------------------------------------------------------
inline double CInputSimulation::getMinLaunchElevationAngle(void) const
{
  return launch_v_inf[VINF_ELE_MIN];
} // double CInputSimulation::getMinLaunchElevationAngle(void) const


//------------------------------------------------------------------------------
inline void CInputSimulation::setMinLaunchElevationAngle(const double angle)
{
  launch_v_inf[VINF_ELE_MIN] = angle;
} // void CInputSimulation::setMinLaunchElevationAngle(const double angle)


//------------------------------------------------------------------------------
inline double CInputSimulation::getMaxLaunchElevationAngle(void) const
{
  return launch_v_inf[VINF_ELE_MAX];
} // double CInputSimulation::getMaxLaunchElevationAngle(void) const


//------------------------------------------------------------------------------
inline void CInputSimulation::setMaxLaunchElevationAngle(const double angle)
{
  launch_v_inf[VINF_ELE_MAX] = angle;
} // void CInputSimulation::setMaxLaunchElevationAngle(const double angle)


//------------------------------------------------------------------------------
inline int CInputSimulation::getMaxStepsCount(void) const
{
  return max_steps;
} // int CInputSimulation::getMaxStepsCount(void) const


//------------------------------------------------------------------------------
inline void CInputSimulation::setMaxStepsCount(const int count)
{
  max_steps = count;
} // void CInputSimulation::setMaxStepsCount(const int count)


//------------------------------------------------------------------------------
inline int CInputSimulation::getMinOutputPoints(void) const
{
  return min_output_points;
} // int CInputSimulation::getMinOutputPoints(void) const


//------------------------------------------------------------------------------

inline void CInputSimulation::setMinOutputPoints(const int n_min_points)
{
  min_output_points = n_min_points;
} // void CInputSimulation::setMinOutputPoints(const int n_min_points)


//------------------------------------------------------------------------------

inline bool CInputSimulation::useDynamicAccuracyControl(void) const
{
  return dac_usage;
} // bool CInputSimulation::useDynamicAccuracyControl(void) const


//------------------------------------------------------------------------------
inline void CInputSimulation::setDynamicAccuracyControl(const bool use_or_not)
{
  dac_usage = use_or_not;
} // void setDynamicAccuracyControl(const bool use_or_not)


//------------------------------------------------------------------------------
inline double CInputSimulation::getDACMaxHypercubeSize(void) const
{
  return dac_hc_start_size;
} // double getDACMaxHypercubeSize(void) const


//------------------------------------------------------------------------------
inline void CInputSimulation::setDACMaxHypercubeSize(const double val)
{
  dac_hc_start_size = val;
} // void setDACMaxHypercubeSize(const double val)


//------------------------------------------------------------------------------
inline double CInputSimulation::getDACMinHypercubeSize(void) const
{
  return dac_hc_final_size;
} // double getDACMinHypercubeSize(void) const


//------------------------------------------------------------------------------
inline void CInputSimulation::setDACMinHypercubeSize(const double val)
{
  dac_hc_final_size = val;
} // void setDACMinHypercubeSize(const double val)


//------------------------------------------------------------------------------
inline double CInputSimulation::getDACMaxFactor(void) const
{
  return dac_factor_max;
} // double getDACMaxFactor(void) const


//------------------------------------------------------------------------------
inline void CInputSimulation::setDACMaxFactor(const double val)
{
  dac_factor_max = val;
} // void setDACMaxFactor(const double val)


//------------------------------------------------------------------------------
inline bool CInputSimulation::getStatusDSSC(void) const
{
  return dssc_usage;
}

//------------------------------------------------------------------------------
inline void CInputSimulation::setStatusDSSC(const bool status)
{
  dssc_usage = status;
}

//------------------------------------------------------------------------------
inline void CInputSimulation::enableDSSC(void)
{
  setStatusDSSC(true);
}

//------------------------------------------------------------------------------
inline void CInputSimulation::disableDSSC(void)
{
  setStatusDSSC(false);
}

//------------------------------------------------------------------------------
inline bool CInputSimulation::getStatusDSSCStepAngleControl(void) const
{
  return dssc_use_step_angle_control;
}

//------------------------------------------------------------------------------
inline void CInputSimulation::setStatusDSSCStepAngleControl(const bool status)
{
  dssc_use_step_angle_control = status;
} 

//------------------------------------------------------------------------------
inline void CInputSimulation::enableDSSCStepAngleControl(void)
{
  setStatusDSSCStepAngleControl(true);
}

//------------------------------------------------------------------------------
inline void CInputSimulation::disableDSSCStepAngleControl(void)
{
  setStatusDSSCStepAngleControl(false);
}

//------------------------------------------------------------------------------
inline void CInputSimulation::setDSSCMaxStepAngle(const double angle)
{
  dssc_max_step_angle = angle;
}

//------------------------------------------------------------------------------
inline double CInputSimulation::getDSSCMaxStepAngle(void) const
{
  return dssc_max_step_angle;
}

//------------------------------------------------------------------------------
inline bool CInputSimulation::getStatusDSSCStepDistanceControl(void) const
{
   return dssc_use_step_distance_control;
}

//------------------------------------------------------------------------------
inline void CInputSimulation::setStatusDSSCStepDistanceControl(const bool status)
{
  dssc_use_step_distance_control = status;
}

//------------------------------------------------------------------------------
inline void CInputSimulation::enableDSSCStepDistanceControl(void)
{
  setStatusDSSCStepDistanceControl(true);
}

//------------------------------------------------------------------------------
inline void CInputSimulation::disableDSSCStepDistanceControl(void)
{
  setStatusDSSCStepDistanceControl(false);
}

//------------------------------------------------------------------------------
inline void CInputSimulation::setDSSCMaxStepDistance(const double distance)
{
  dssc_max_step_distance = distance;
}

//------------------------------------------------------------------------------
inline double CInputSimulation::getDSSCMaxStepDistance(void) const
{
  return dssc_max_step_distance;
}

//------------------------------------------------------------------------------
inline bool CInputSimulation::getStatusDSSCApproachControl(void) const
{
  return dssc_use_approach_control;
}

inline void CInputSimulation::setStatusDSSCApproachControl(const bool status)
{
  dssc_use_approach_control = status;
}

//------------------------------------------------------------------------------
inline void CInputSimulation::enableDSSCApproachControl(void)
{
  setStatusDSSCApproachControl(true);
}

//------------------------------------------------------------------------------
inline void CInputSimulation::disableDSSCApproachControl(void)
{
  setStatusDSSCApproachControl(false);
}

//------------------------------------------------------------------------------
inline void CInputSimulation::setDSSCApproachFactor(const double factor)
{
  dssc_approach_step_size_factor = factor;
} // void CInputSimulation::setDSSCApproachFactor(const double factor)

//------------------------------------------------------------------------------
inline double CInputSimulation::getDSSCApproachFactor(void) const
{
  return dssc_approach_step_size_factor;
} // double CInputSimulation::getDSSCApproachFactor(void) const

//------------------------------------------------------------------------------
inline double CInputSimulation::getDSSCMinStepSize(void) const
{
  return dssc_min_step_size;
} // double CInputSimulation::getDSSCMinStepSize(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setDSSCMinStepSize(const double step_size)
{
  dssc_min_step_size = step_size;
} // void CInputSimulation::setDSSCMinStepSize(const double step_size)

//------------------------------------------------------------------------------
inline double CInputSimulation::getDSSCMaxStepSize(void) const
{
  return dssc_max_step_size;
} // double CInputSimulation::getDSSCMaxStepSize(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setDSSCMaxStepSize(const double step_size)
{
  dssc_max_step_size = step_size;
} // void CInputSimulation::setDSSCMaxStepSize(const double step_size)

//------------------------------------------------------------------------------
inline StateType CInputSimulation::getInitialStateType(void) const
{
  return initial_state;
} // CInputSimulation::StateType getInitialStateType(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setInitialStateType(
  const StateType state_type)
{
  initial_state = state_type;
} // void setInitialStateType(const StateType state_type)

//------------------------------------------------------------------------------
inline BodyID CInputSimulation::getInitialBody(void) const
{
  return initial_body;
} // BodyID CInputSimulation::getInitialBody(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setInitialBody(const BodyID body_id)
{
  initial_body = body_id;
} // void CInputSimulation::setInitialBody(const BodyID body_id)

//------------------------------------------------------------------------------
inline double CInputSimulation::getInitialOrbitElement(
  const SpaceLib::KEPLERIAN_ELEMENT element) const
{
    return elements_of_initial_orbit[element];
} // inline double CInputSimulation::getInitialOrbitElement(...) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setInitialOrbitElement(
  const SpaceLib::KEPLERIAN_ELEMENT element, const double value)
{
    elements_of_initial_orbit[element] = value;
} // CInputSimulation::setInitialOrbitElement(void)

//------------------------------------------------------------------------------
inline bool CInputSimulation::useTargetOrbitalElement(
  const KEPLERIAN_ELEMENT e) const
{
#ifdef development_version
  assert( size_t(e) < use_target_orbital_elements.size() );
#endif
  return use_target_orbital_elements[e];
} // bool CInputSimulation::useTargetOrbitalElement(void) const

//------------------------------------------------------------------------------
inline double CInputSimulation::getTargetOrbitalElement(
  const KEPLERIAN_ELEMENT e) const
{
#ifdef development_version
  assert( size_t(e) < target_orbital_elements.size() );
#endif
  return target_orbital_elements[e];
} // double CInputSimulation::getTargetOrbitalElement(void) const

//------------------------------------------------------------------------------
inline double CInputSimulation::getTargetOrbitalElementMaxDeviation(
    const KEPLERIAN_ELEMENT e) const
{
#if defined development_version
  assert( size_t(e) < max_target_orbital_element_deviation.size() );
#endif
  return max_target_orbital_element_deviation[e];
} // double CInputSimulation::getTargetOrbitalElementMaxDeviation(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setTargetOrbitalElementMaxDeviation(
  const KEPLERIAN_ELEMENT e, const double val)
{
#if defined development_version
  assert( size_t(e) < max_target_orbital_element_deviation.size() );
#endif
  max_target_orbital_element_deviation[e] = val;
} // void CInputSimulation::setTargetOrbitalElementMaxDeviation(void)

//------------------------------------------------------------------------------
inline BodyID CInputSimulation::getTargetBodyId(void) const
{
  return target_body;
} // BodyID CInputSimulation::getTargetBodyId(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setTargetBodyId(const BodyID target_id)
{
  target_body = target_id;
} // void CInputSimulation::setTargetBodyId(const BodyID target_id)

//------------------------------------------------------------------------------
inline double CInputSimulation::getTargetOrbitalElementProxWeight(
  const KEPLERIAN_ELEMENT e) const
{
#if defined development_version
  assert( size_t(e) < prox_weights.size() );
#endif
  return prox_weights[e];
} // double CInputSimulation::getTargetOrbitalElementProxWeight(...) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setTargetOrbitalElementProxWeight(
                                       const KEPLERIAN_ELEMENT e
                                       ,const double val)
{
#if defined development_version
  assert( size_t(e) < prox_weights.size() );
#endif
  prox_weights[e] = val;
} // CInputSimulation::setTargetOrbitalElementProxWeight(...)

//------------------------------------------------------------------------------
inline TrajectoryType CInputSimulation::getTrajectoryType(void) const
{
  return trajectory;
} // TrajectoryType CInputSimulation::getTrajectoryType(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setTrajectoryType(const TrajectoryType traj_type)
{
  trajectory = traj_type;
} // void setTrajectoryType(const TrajectoryType traj_type)

//------------------------------------------------------------------------------
inline OptimizationGoalType CInputSimulation::getOptGoalType(void) const
{
  return opt_goal;
} // OptimizationGoalType CInputSimulation::getOptGoalType(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setOptGoalType(const OptimizationGoalType nog)
{
  opt_goal = nog;
} // void setOptGoalType(const OptimizationGoalType opt_goal)

//------------------------------------------------------------------------------
inline StateType CInputSimulation::getTargetStateType(void) const
{
  return target_state;
} // StateType CInputSimulation::getTargetStateType(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setTargetStateType(const StateType state_type)
{
  target_state = state_type;
} // void setTargetStateType(const StateType)

//------------------------------------------------------------------------------
inline const std::valarray<double>& 
CInputSimulation::getTargetOrbitalElementProxWeights(void) const
{
  return prox_weights;
} // std::valarray<double>& InputSimulation::getTargetOrbitalElementProxWeights(void) const

//------------------------------------------------------------------------------
inline std::valarray<double>& 
CInputSimulation::getTargetOrbitalElementProxWeights(void)
{
  return prox_weights;
} // std::valarray<double>& getTargetOrbitalElementProxWeights(void)

//------------------------------------------------------------------------------
inline DESSolverType CInputSimulation::getIntegratorType(void) const
{
  return des_solver;
} // 

//------------------------------------------------------------------------------
inline void CInputSimulation::setIntegratorType(const DESSolverType integrator_type)
{
  des_solver = integrator_type;
} // void CInputSimulation::setIntegratorType(const DESSolverType integrator_type)

//------------------------------------------------------------------------------
inline double CInputSimulation::getMaxRelIntegrationError(void) const
{
  return max_error_rel;
} // double CInputSimulation::getMaxRelIntegrationError(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setMaxRelIntegrationError(const double rel_error)
{
  max_error_rel = rel_error;
} // void CInputSimulation::setMaxRelIntegrationError(const double rel_error)

//------------------------------------------------------------------------------
inline double CInputSimulation::getMaxAbsIntegrationError(void) const
{
  return max_error_abs;
} // double CInputSimulation::getMaxAbsIntegrationError(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setMaxAbsIntegrationError(const double abs_error)
{
  max_error_abs = abs_error;
} // void CInputSimulation::setMaxAbsIntegrationError(const double rel_error)


//------------------------------------------------------------------------------
inline double CInputSimulation::getMinSolarDistance(void) const
{
  return min_solar_distance;
} // double CInputSimulation::getMinSolarDistance(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setMinSolarDistance(const double distance)
{
  min_solar_distance = distance;
} // void CInputSimulation::setMinSolarDistance(const double)

//------------------------------------------------------------------------------
inline double CInputSimulation::getMaxFinalRelDistance(void) const
{
  return max_dr;
} // double CInputSimulation::getMaxFinalRelDistance(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setMaxFinalRelDistance(const double value)
{
  max_dr = value;
} // void CInputSimulation::setMaxFinalRelDistance(const double)

//------------------------------------------------------------------------------
inline double CInputSimulation::getMaxFinalRelVelocity(void) const
{
  return max_dv;
} // double CInputSimulation::getMaxFinalRelVelocity(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setMaxFinalRelVelocity(const double value)
{
  max_dv = value;
} // void CInputSimulation::setMaxFinalRelVelocity(const double)

//------------------------------------------------------------------------------
inline const MathLib::Vector& CInputSimulation::getInitPosVec(void) const
{
  return initial_position;
} // const MathLib::Vector& getInitPosVec(void)

//------------------------------------------------------------------------------
inline void CInputSimulation::setInitPosVec(const MathLib::Vector& pos)
{
#if defined development_version
  assert(pos.size() == 3);
#endif
  initial_position = pos;
} // void setInitPosVec(const MathLib::Vector& pos)

//------------------------------------------------------------------------------
inline void CInputSimulation::setInitPosVec(double p1, double p2, double p3)
{
  if (initial_position.size()!=3) initial_position.Set(MathLib::Vector(3));
  initial_position[0] = p1;
  initial_position[1] = p2;
  initial_position[2] = p3;
} // void setInitPosVec(double p1, double p2, double p3)

//------------------------------------------------------------------------------
inline const MathLib::Vector& CInputSimulation::getInitVelVec(void) const
{
  return initial_velocity;
} // const MathLib::Vector& getInitVelVec(void)

//------------------------------------------------------------------------------
inline void CInputSimulation::setInitVelVec(const MathLib::Vector& vel)
{
#if defined development_version
  assert(vel.size() == 3);
#endif
  initial_velocity = vel;
} // void setInitVelVec(const MathLib::Vector& vel)

//------------------------------------------------------------------------------
inline void CInputSimulation::setInitVelVec(double v1, double v2, double v3)
{
  if (initial_velocity.size()!=3) initial_velocity.Set(MathLib::Vector(3));
  initial_velocity[0] = v1;
  initial_velocity[1] = v2;
  initial_velocity[2] = v3;
} // void setInitVelVec(double v1, double v2, double v3)

//------------------------------------------------------------------------------
inline const std::string& CInputSimulation::getInitCentralBodyName(void) const
{
  return initial_central_body_name;
} // const std::string& getCentralBodyName(void) const

//------------------------------------------------------------------------------
inline std::string& CInputSimulation::getInitCentralBodyName(void)
{
  return initial_central_body_name;
} // const std::string& getCentralBodyName(void) const

//------------------------------------------------------------------------------
inline BodyID CInputSimulation::getInitCentralBody(void) const
{
  return initial_central_body;
} // void setInitCentralBodyName(const std::string& name)

//------------------------------------------------------------------------------
inline void CInputSimulation::setInitCentralBody(const BodyID central_body)
{
  initial_central_body = central_body;
} // void setInitCentralBodyName(const std::string& name)

//------------------------------------------------------------------------------
inline const SpaceLib::CSolarSystem* CInputSimulation::getSolarSystem(void)const
{
  return m_solsys;
} // CSolarSystem* getSolarSystem(void)

//------------------------------------------------------------------------------
inline SpaceLib::CSolarSystem* CInputSimulation::getSolarSystem(void)
{
  return m_solsys;
} // CSolarSystem* getSolarSystem(void)

//------------------------------------------------------------------------------
inline void CInputSimulation::setSolarSystem(CSolarSystem* sol_sys)
{
  m_solsys = sol_sys;
} // void setSolarSystem(const CSolarSystem* sol_sys)

//------------------------------------------------------------------------------
inline bool CInputSimulation::isHeliocentric(void) const
{
  return !body_centric;
} // bool isHeliocentric(void) const

//------------------------------------------------------------------------------
inline bool CInputSimulation::isBodycentric(void) const
{
  return !isHeliocentric();
} // bool isPlanetocentric(void) const

//------------------------------------------------------------------------------
inline void CInputSimulation::setHeliocentricSimulation(bool helio_centric)
{
  body_centric = !helio_centric;
} // void setHeliocentricSimulation(bool helio_centic)

//------------------------------------------------------------------------------
inline CBody* CInputSimulation::getGTOC2Body(const int group,const unsigned int no) const
{
#if defined development_version
  assert(uint(group) < gtoc2_bodies.size());  
#endif // #if defined development_version
  if (no < (gtoc2_bodies[group]).size()) return (gtoc2_bodies[group])[no];
  else return 0;
}

#endif // CINPUTSIMULATION_H_INCLUDED
