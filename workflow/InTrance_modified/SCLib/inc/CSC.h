/*! \file CSC.h \brief A class that provides the functionality
for a spacecraft.

\version

1.7

\author

Bernd Dachwald\n German Aerospace Center (DLR)\n
Institute of Space Simulation\n Linder Hoehe\n 51147 Cologne\n Germany\n
bernd.dachwald@dlr.de

\author
Andreas Ohndorf\n German Aerospace Center (DLR)\n
Institute of Space Simulation\n Linder Hoehe\n 51147 Cologne\n Germany\n
andreas.ohndorf@gmx.de

\date

22.10.2006 (AO)

\b Compiler:

Microsoft Visual C++ .NET

\b Copyright:

(c) 2001-2006 by Bernd Dachwald, Andreas Ohndorf

This software is protected by national and international copyright.
Any unauthorized use, reproduction or modificaton is unlawful and
will be prosecuted. Commercial and non-private application of the
software in any form is strictly prohibited unless otherwise granted
by the authors.

The code is provided without any warranty; without even the implied
warranty of merchantibility or fitness for a particular purpose.

*/

#ifndef CSC_H_INCLUDED
#define CSC_H_INCLUDED

#include <assert.h>
#include <SCConst.h>
#include <CInputSC.h>
#include <CSCAttitude.h>
#include <CSCSteering.h>
#include <CNEP.h>
#include <CSEP.h>
#include <CNSTAR.h>
#include <CT6.h>
#include <CRIT.h>
#include <CSolarSail.h>
#include <CPropSys.h>

const double INITIAL_SRD = 0.0; // start value of solar radiation dose

class CNC; // forward declaration

/// A class that provides the functionality for a low-thrust spacecraft.
class CSC : public CBody
{
  /*!
  \todo AO:remove moment of inertia and steering moments members and related methods
  */

  friend void CSCSteering::setNCInputCartesian
    (int, int, const CState&, const CState&, CSC*, const CState*);

  friend void CSCSteering::setNCInputPolar
    (int, int, const CState&, const CState&, CSC*, const CState*);

  friend void CSCSteering::setNCInputKepler
    (int, int, const CState&, const CState&, CSC*, const CState*);

  friend void CSCSteering::setNCInputPropellant
    (int, CSC*);

  friend void CSCSteering::performSteering
    (CSC*, CNC*nc, const CState&, const CState&, const CState*);

  //---------------------------------------------------------------------------

private:

  double m_scmass;                                   // spacecraft (dry) mass
  static unsigned long int m_count;

public:

  double getSCMass(void) const        { return m_scmass; };
  void   setSCMass(const double mass) { m_scmass=mass; };

private:

  MathLib::Vector m_moi;      ///< moments of inertia

public:

  /// get moments of inertia
  const MathLib::Vector& getMOI(void) const {return m_moi;};
  /// set moments of inertia
  void setMOI(const MathLib::Vector& moi) {(moi.size()==3)?m_moi=moi:0;};
  /// set moments of inertia
  void setMOI(const std::valarray<double>& moi)
  {
    if (moi.size()==3)
    {
      m_moi[0]=moi[0];m_moi[1]=moi[1];m_moi[2]=moi[2];
    }
  } // inline void setMOI(const std::valarray<double>& moi)

private:
  MathLib::Vector m_msm;    ///< max. steering moments about principle axes

public:
  /// get moments of inertia
  const MathLib::Vector& getMaxSteerMoments(void) const {return m_msm;};
  /// set moments of inertia
  void setMaxSteerMoments(const MathLib::Vector& msm) {(m_msm.size()==3)?m_msm=msm:0;};
  /// set moments of inertia
  inline void setMaxSteerMoments(const std::valarray<double>& msm)
  {
    if (msm.size()==3){m_msm[0]=msm[0];m_msm[1]=msm[1];m_msm[2]=msm[2];}
  } // inline void setMOI(const std::valarray<double>& moi)
private:

  Vector m_thrustacc;          // thrust acceleration vector (3D) [m/s^2]
  Vector m_delta_v;            // Delta V vector (3D) [m/s]
  Vector m_dkepler;            // time derivatives of Keplerian elements (6D)

  int    m_n_ctrl_input_dim; // number of NC input neurons
  Vector m_input_ncp;           // NC input

public:

  /// controller input dimension (# of values)
  int getCtrlInputDim(void) const;

  /// sets controller input dimension (# of controller input values)
  void setCtrlInputDim(const int input_dim);

  /// set val as new input vector value at position n
  double setCtrlInput(const int n, const double val);

  /// set complete control input vector
  void setCtrlInputVec(const MathLib::Vector& vec);

  /// input vector value at index n
  double getCtrlInput(const int n) const;

  /// controller input vector (const version)
  const MathLib::Vector& getCtrlInputVec(void) const;

  /// controller input vector
  MathLib::Vector& getCtrlInputVec(void);

private:

  int m_n_ctrl_output_dim;   // number of NC output neurons
  Vector m_nc_output;          // NC output

public:
  /// controller output dimension (# of output values)
  int getCtrlOutputDim(void) const;
  /// set controller output dimension (# of output values)
  void setCtrlOutputDim(const int output_dim);
  /// controller output vector value at position n
  double getCtrlOutput(const int n) const;
  /// get controller output vector (const)
  const MathLib::Vector& getCtrlOutputVec(void) const;
  /// get controller output vector
  MathLib::Vector& getCtrlOutputVec(void);
  /// set controller output vector value at position n
  void setCtrlOutput(const int n, const double val);

private:

  int m_DEcalls;               // DE function call counter

  CSCAttitude m_attitude;      // spacecraft attitude object

  CSCSteering m_steering;      // spacecraft steering object

public:

  CSCSteering getSteering(void)const { return m_steering; };

private:

  CPropSys* m_propsys;         // propulsion system object

public:

  /*!
  \brief Gets propulsion system object.

  \return  Pointer to object of class CPropulsionSystem (or a class
  derived from it).
  */
  CPropSys * getPropSys(void) const { return m_propsys; };


  /*!
  \brief Sets CSC object to have new propulsion system object of given type.

  \param prop_sys_type New propulsion system type.
  */
  void setPropSys(PropulsionSystemType prop_sys_type);


  /*!
  \brief Sets new propulsion system object.

  Attention:\n The object pointer will be set as a member field of CSC 
  object.
  So if you want to use a copy of a certain CPropSys object as your new
  propulsion system, you may use this method to set a copy of it, which
  you can get by calling the clone()-method on your CPropSys - object.\n
  \n Example:\n\n
  myCSC.setPropSys(anotherCSC.getPropSys()->clone()); // right\n
  myCSC.setPropSys(anotherCSC.getpropSys());          // !WRONG!\n

  \param new_prop_sys Pointer to new object of class CPropSys
  (or a class derived from it)
  */
  void setPropSys(CPropSys * new_prop_sys);

public:

  // copy operator
  CSC& operator=(const CSC& sc) 
  {
    CBody::operator=(sc); // copy operator of parent class

    // copy local data fields
    m_scmass              = sc.m_scmass;
    m_moi                 = sc.m_moi;
    m_thrustacc           = sc.m_thrustacc;
    m_delta_v             = sc.m_delta_v;
    m_dkepler             = sc.m_dkepler;
    m_n_ctrl_input_dim    = sc.m_n_ctrl_input_dim;
    m_input_ncp.Set(sc.m_input_ncp);
    m_n_ctrl_output_dim   = sc.m_n_ctrl_output_dim;
    m_nc_output.Set(sc.m_nc_output);
    m_DEcalls             = sc.m_DEcalls;
    m_attitude            = sc.m_attitude;
    m_steering            = sc.m_steering;

    // copy propulsion system
    if (m_propsys != 0 && (m_propsys->m_type==sc.m_propsys->m_type))
    {
      switch (sc.m_propsys->m_type)
      {
      case NEP:
        (*(dynamic_cast<CNEP*>(m_propsys))) = 
          (*(dynamic_cast<CNEP*>(sc.m_propsys)));
        break;

      case NSTAR:
        (*(dynamic_cast<CNSTAR*>(m_propsys))) = 
          (*(dynamic_cast<CNSTAR*>(sc.m_propsys)));
        break;

      case RIT:
        (*(dynamic_cast<CRIT*>(m_propsys))) = 
          (*(dynamic_cast<CRIT*>(sc.m_propsys)));
        break;

      case SOLARSAIL:
        (*(dynamic_cast<CSolarSail*>(m_propsys))) = 
          (*(dynamic_cast<CSolarSail*>(sc.m_propsys)));
        break;

      case T6:
        (*(dynamic_cast<CT6*>(m_propsys))) = 
          (*(dynamic_cast<CT6*>(sc.m_propsys)));
        break;

      case UNDEFINED_PST:
      default:
        assert(0);
        break;
      }
    }
    else if (sc.m_propsys != 0)
    {        
      setPropSys(sc.m_propsys->clone(this));        
    }
    m_propsys->m_state=this;

    return *this;
  };


  /// equal operator
  bool operator==(const CSC& rhs)
  {
    return (
      (m_scmass              == rhs.m_scmass) &&
      (m_moi                 == rhs.m_moi) &&
      (m_thrustacc           == rhs.m_thrustacc) &&
      (m_delta_v             == rhs.m_delta_v) &&
      (m_dkepler             == rhs.m_dkepler) &&
      (m_n_ctrl_input_dim    == rhs.m_n_ctrl_input_dim) &&
      (m_input_ncp           == rhs.m_input_ncp) &&
      (m_n_ctrl_output_dim   == rhs.m_n_ctrl_output_dim) &&
      (m_nc_output           == rhs.m_nc_output) &&
      (m_DEcalls             == rhs.m_DEcalls) &&
      (m_attitude            == rhs.m_attitude) &&
      (m_steering            == rhs.m_steering) && 
      ((*m_propsys)          == (*(rhs.m_propsys))) &&
      (CBody::operator==(rhs))
      );
  }


private:

  //-----------------------------------------------------------------





  // ~~~~~ begin inline



  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /// radial thrust acceleration component
  double ThrustAcc_r(void)const { return m_thrustacc(0); };


  double ThrustAcc_t(void)const { return m_thrustacc(1); };
  // purpose: get transversal thrust acceleration component

  double ThrustAcc_h(void)const { return m_thrustacc(2); };
  // purpose: get orbit normal thrust acceleration component

  double ThrustAcc_th(void)const {
    return sqrt(m_thrustacc(1)*m_thrustacc(1)+m_thrustacc(2)*m_thrustacc(2));
  };
  // purpose: get thrust acceleration vector component perpendicular
  // to radial direction

  double DeltaV_r(void)const { return m_delta_v(0); };
  // purpose: get radial Delta V

  double DeltaV_t(void)const { return m_delta_v(1); };
  // purpose: get transversal Delta V

  double DeltaV_h(void)const { return m_delta_v(2); };
  // purpose: get orbit normal Delta V

  double DeltaV_th(void)const {
    return sqrt(m_delta_v(1)*m_delta_v(1)+m_delta_v(2)*m_delta_v(2));
  };
  // purpose: get Delta V perpendicular to radial direction

  double DeltaV(void)const { return Norm(m_delta_v); };
  // purpose: get total Delta V

  // ~~~~~ end inline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // __update__

  void UpdateThrustAcc(void);
  // purpose: update thrust acceleration vector

  void UpdatedElements(void);
  // purpose: update time derivatives of the keplerian elements

  void UpdateDeltaV(const double delta_t);
  // time interval
  // purpose: update delta V

  //=============================================================================
  // CSC CLASS INTERFACE
  //=============================================================================
public:

  // ~~~~~ begin inline



  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





  // __accessors__


  double ThrustAcc(void)const { return Norm(m_thrustacc); };
  // purpose: get magnitude of thrust acceleration vector

  double getMass(void)const { return m_scmass+getPropellantMass(); };
  // purpose: get actual (wet) spacecraft mass

  double getClockAngle(void)const { return m_attitude.getClockAngle(); };
  // purpose: get actual clock angle

  double getConeAngle(void)const { return m_attitude.getConeAngle(); };
  // purpose: get actual cone angle

  Vector dElements(void)const { return m_dkepler; };
  // purpose: get 6D-Vector with the time derivatives of the keplerian elements

  int getNumberOfLSL(void)const { return m_steering.getNumberOfLSL(); };
  // purpose: get number of local steering laws

  int getNumberOfLStSt(void)const { return m_steering.getNumberOfLStSt(); };
  // purpose: get number of steering strategies

  StStIndex getSteeringStrategy(void)const {
    return m_steering.getSteeringStrategy();
  };
  // purpose: get index of actual steering strategy

  // __mutators__

  void IncrementDECounter(void) { m_DEcalls++; };
  // purpose: increment DE function call counter

  // __spacecraft steering__

  void setAttitude(const double alpha, const double beta);

  /*! \brief set clock and cone angle according to the optimal thrust vector
  * \param xi_r radial component of the optimal thrust vector
  * \param xi_t transversal component of the optimal thrust vector
  * \param xi_h orbit normal component of the optimal thrust vector
  */
  void setAttitude(const double xi_r,const double xi_t,const double xi_h);

  void setSteeringStrategy(const StStIndex stst)
    // steering strategy index

  { m_steering.setSteeringStrategy(stst); };
  // purpose: set actual steering strategy

  /*! \brief determine new steering commando
  * \param nc pointer to controller to use
  * \param target_state actual state of target body (neurocontroller input)
  */
  void performSteering(CNC* nc, const CState& target_state);

  // purpose: do spacecraft steering with neurocontroller
  void performSteering(const Vector&lsl_vector,
    // local steering law vector
    const CState&target_state)
    // actual state of target body
  { m_steering.performSteering(this, lsl_vector, *this, target_state); };
  // purpose: do spacecraft steering with local steering law vector

  void performSteering(const double clock_angle,
    // clock angle
    const double cone_angle)
    // cone angle
  { m_steering.performSteering(this, clock_angle, cone_angle); };
  // purpose: do spacecraft steering with direct input of attitude angles

  // __other functionality__

  void ResetDeltaV(void) { m_delta_v=0.0; };   // purpose: reset delta V

  // ~~~~~ end inline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // __constructors__


  /// CSolarSystem class constructor
  CSC(SpaceLib::CSolarSystem* solSys = 0);
  // purpose: creates spacecraft with undefined propulsion system


  /// Copy constructor
  CSC(const CSC& sc);


  CSC(const PropulsionSystemType psystype, SpaceLib::CSolarSystem* solSys);
  // spacecraft propulsion system
  // purpose: creates spacecraft with defined propulsion system


  CSC(const CInputSC&isc, SpaceLib::CSolarSystem* solSys);
  // CInputSC object
  // purpose: creates spacecraft from CInputSC object

  //  CSC(const char *scp_filename);
  // spacecraft input file
  // purpose: creates spacecraft from spacecraft input file

  /// destructor
  ~CSC(void) 
  {
    SafeDelete(m_propsys);
    --m_count;
  }


  // __update__

  void UpdateAccDeltaVec(const double delta_t);
  // purpose: update thrust acceleration vector, time derivatives of the
  //          keplerian elements and delta V

  // __accessors__

  Vector PolarStatePropellantMassSRDVector(void)const;
  // purpose get 8D-vector with the polar state and the propellant mass

  Vector ThrustAcc_xyz(void)const;
  // purpose: get thrust acceleration in cartesian coordinates



  // __mutators__

  void setPropellantMassSRDVector(const Vector& v, const double t);
  // purpose: set 8D-vector with the polar state, the propellant mass, and the solar radiation dose

  // __spacecraft steering__

  void setSteeringStrategy(const Vector&lsl_vector);
  // LSL vector
  // purpose: set actual steering strategy according to LSL vector

  void setSteeringStrategy(CNC*nc);
  // (pointer to the) neurocontroller
  // purpose: set actual steering strategy according to the used neurocontroller

  void GravityAssistSteeringOn(void) { m_steering.m_ga_steering=true; };
  // purpose: turn on gravity assist steering

  void GravityAssistSteeringOff(void) { m_steering.m_ga_steering=false; };
  // purpose: turn off gravity assist steering

  bool GravityAssistSteering(void) { return m_steering.m_ga_steering; };
  // purpose: check whether gravity assist steering is used or not


  // __other functionality__

  double MaxThrustAcceleration(void)const;
  // purpose: get maximum thrust acceleration (for use in VRML)

  bool CheckSolarDistance(const double min_solar_distance) const;
  // minimum sun-spacecraft distance

  // purpose: check, whether the spacecraft is within the minimum
  //          sun-spacecraft distance

  bool CheckTemperature(void)const;
  // purpose: check, whether the spacecraft is too hot
  //          (only necessary for real solar sails)

  void Init(const double srd_or_mp = 0.0);

  void setLaunchState(
    const double launch_time,                   // launch time
    const BodyID launch_body,              // launch body
    const double launch_velocity,               // launch velocity
    const double launch_azimuth,                // azimuth for launch
    const double launch_elevation               // elevation for launch
    );
  // purpose: set launch state for the spacecraft according to a specified
  //          launch body

  void setLaunchState(
    const double launch_time,           // launch time
    const double a,                     // orbital elements
    const double e,
    const double i,
    const double Omega,
    const double omega,
    const double M,
    const double launch_velocity,       // launch velocity
    const double launch_azimuth,        // azimuth for launch
    const double launch_elevation       // elevation for launch
    );
  // purpose: set launch state for the spacecraft according to a specified
  //          set of orbital elements

  void setDisturbingAcceleration(const double magnitude);
  // purpose: tbd

  void CheckRemainingPropellant(const double dt);
  // purpose: check if propellant tank is empty


  // __input/output functions__

  void PrintHeader(
    std::ofstream&out,       // ofstream object
    const char*comment=NULL  // optional comment
    ) const;
  // purpose: print header for spacecraft data

  void PrintDeltaHeader(
    std::ofstream&out,                  // ofstream object
    const char*comment=NULL             // optional comment
    ) const;
  // purpose: print header for "Delta data"

  void PrintNCInputHeader(std::ofstream&out) const;

  void PrintNCOutputHeader(std::ofstream&out) const;

  /*! \brief write csv-file header of delta-r and delta-v
  \param out standard output file stream object
  \param comment optional comment to add
  */
  void PrintProximityHeader(std::ofstream& out,const char*comment = 0) const;

  void Print(std::ofstream&out) const;
  // ofstream object
  // purpose: print spacecraft data

  void PrintDelta(
    std::ofstream&out,                  // ofstream object
    const int planet,                   // planet, for which the data is valid
    bool spacecraft_minus_planet=false
    // if true : spacecraft value - planet value
    // if false: planet value - spacecraft value
    ) const;
  // purpose: print "Delta data"

  void PrintNCInput(std::ofstream&out) const;

  void PrintNCOutput(std::ofstream&out) const;

  void PrintProximity(
    std::ofstream&out,                  // ofstream object
    const int planet                    // planet, for which the data is valid
    ) const;
  // purpose: print Delta r and Delta v

  /*! \brief write VRML-script code into given output stream object to
   *display propulsion system type-dependent spacecraft shape
   *
   * The appearance of every spacecraft type, more precisely the propulsion
   * system type used by the spacecraft to fly through space, shall be
   * distinguishable in the possibly produced VRML-script file. In order to
   * support this, a CSC class object is able to produce it's own appearance
   * in VRML-code and writes it into the given output stream, which may be a
   * text stream or a file stream.\
   * Basically it calls a similar method of the member object representing
   * the propulsion system and this finally produces the VRML-code.
   * Consequently, this is also the place to change an Spacecraft's object
   * VRML-shape by overloading printVrmlAppearance in the appropriate propulsion
   * system object.

   * \param os reference to given output stream object
   * \throws nothing
   * \return nothing
   * \sa CPropSys::printVrmlAppearance*/
  void printVrmlAppearance(std::ostream& os) const throw ();


  // ************************************************************************
  // ** FUNCTIONS ONLY FOR THE PROPULSION SYSTEM                           **
  // ************************************************************************

  PropulsionSystemType getPropulsionSystemType(void) const {
    return m_propsys->getPropulsionSystemType();
  };

  bool isSolarSail(void) const { return m_propsys->isSolarSail(); };
  bool isNEP(void)       const { return m_propsys->isNEP(); };
  bool isNSTAR(void)     const { return m_propsys->isNSTAR(); };
  bool isT6(void)        const { return m_propsys->isT6(); };
  bool isRIT(void)       const { return m_propsys->isRIT(); };
  bool isSEP(void)       const { return m_propsys->isSEP(); };
  bool isEP(void)        const { return m_propsys->isEP(); };


  void setInitialPropellantMassInterval(
    const double min_mass,
    const double max_mass) {

      m_propsys->setInitialPropellantMassInterval(min_mass, max_mass);
  };


  double getInitialPropellantMass(void)const {
    return m_propsys->getInitialPropellantMass();
  };


  void setInitialPropellantMass(const double mass) {
    m_propsys->setInitialPropellantMass(mass);
  };


  double getPropellantMass(void)const {
    return m_propsys->getPropellantMass();
  };


  void setPropellantMass(const double mass) {
    m_propsys->setPropellantMass(mass);
  };


  void setSolarRadiationDose(const double srd) {
    m_propsys->setSolarRadiationDose(srd);
  };


  void setEmpty(const bool yes_or_no) { m_propsys->setEmpty(yes_or_no); };

  bool isEmpty(void)const { return m_propsys->isEmpty(); };


  // ========================================
  // == FUNCTIONS ONLY FOR SOLAR SAILCRAFT ==
  // ========================================

  // __ mutators__

  void CSC::setLightness(const double lambda) {
    m_propsys->setLightness(lambda);
  };


  void CSC::setOpticalModel(
    const bool   ideal_or_real,
    const double srd = 0.0,
    const double max_temperature = HUGE_VAL,
    const double rho_0 = 0.0,
    const double s_0 = 0.0,
    const double epsilonf_0 = 0.0,
    const double epsilonb_0 = 0.0,
    const double rho_infty = 0.0,
    const double s_infty = 0.0,
    const double epsilonf_infty = 0.0,
    const double epsilonb_infty = 0.0,
    const double rho_half = 0.0,
    const double s_half = 0.0,
    const double epsilonf_half = 0.0,
    const double epsilonb_half = 0.0,
    const double Bf = 0.0,
    const double Bb = 0.0)
  {
    m_propsys->setOpticalModel(
      ideal_or_real,
      srd,
      max_temperature,
      rho_0,
      s_0,
      epsilonf_0,
      epsilonb_0,
      rho_infty,
      s_infty,
      epsilonf_infty,
      epsilonb_infty,
      rho_half,
      s_half,
      epsilonf_half,
      epsilonb_half,
      Bf,
      Bb);
  };


  // __accessors__

  double CSC::getLightness(void)const { return m_propsys->getLightness(); };
  bool   CSC::isIdealSail(void)const { return m_propsys->isIdealSail(); };
  bool   CSC::isRealSail(void)const { return !m_propsys->isIdealSail(); };

  double CSC::getSolarRadiationDose(void) const {
    return m_propsys->getSolarRadiationDose();
  };

  double CSC::getMaxTemperature(void) const {
    return m_propsys->getMaxTemperature();
  };

  double CSC::getrho(void) const
  { return m_propsys->getrho(); };

  double CSC::gets(void) const { return m_propsys->gets(); };

  double CSC::getepsilonf(void) const {
    return m_propsys->getepsilonf();
  };

  double CSC::getepsilonb(void) const {
    return m_propsys->getepsilonb();
  };

  double CSC::getBf(void)const { return m_propsys->getBf(); };

  double CSC::getBb(void)const { return m_propsys->getBb(); };

  double CSC::getG(void) const { return m_propsys->getG(); };

  double CSC::getH(void) const { return m_propsys->getH(); };

  double CSC::getK(void) const { return m_propsys->getK(); };

  double CSC::getG0(void) const { return m_propsys->getG0(); };

  double CSC::getH0(void) const { return m_propsys->getH0(); };

  double CSC::getK0(void) const { return m_propsys->getK0(); };

  double CSC::getTemperature(void) const {
    return m_propsys->getTemperature(CState::r(), getConeAngle());
  };

  double CSC::getSolarEnergyFlux(void) const {
    return m_propsys->getSolarEnergyFlux(CState::r(), getConeAngle());
  };

  // ===============================================
  // == FUNCTIONS ONLY FOR NEP AND SEP SPACECRAFT ==
  // ===============================================

  void CSC::setThrottle(const double throttle) {
    m_propsys->setThrottle(throttle);
  };


  void CSC::setThrottleType(const ThrottleType throttletype) {
    m_propsys->setThrottleType(throttletype);
  };


  // __accessors__

  double CSC::getIsp(void)const { return m_propsys->getIsp(); };

  double CSC::getThrottle(void)const { return m_propsys->getThrottle(); };

  ThrottleType CSC::getThrottleType(void)const {
    return m_propsys->getThrottleType();
  };


  // =======================================
  // == FUNCTIONS ONLY FOR NEP SPACECRAFT ==
  // =======================================

  void CSC::setMaxThrust(const double thrust_max) {
    m_propsys->setMaxThrust(thrust_max);
  };


  void CSC::setIsp(const double isp) { m_propsys->setIsp(isp); };


  // =======================================
  // == FUNCTIONS ONLY FOR SEP SPACECRAFT ==
  // =======================================

  // __mutators__

  void CSC::setCharacteristicPower(const double characteristic_power) {
    m_propsys->setCharacteristicPower(characteristic_power);
  };


  void CSC::setPowerExponent(const double _p) {
    m_propsys->setPowerExponent(_p);
  };


  void CSC::setNumberOfThrusters(const double thrusters) {
    m_propsys->setNumberOfThrusters(thrusters);
  };


  // __accessors__

  double CSC::getNumberOfThrusters(void)const {
    return m_propsys->getNumberOfThrusters();
  };

  double CSC::getNumberOfUsedThrusters(void)const {
    return m_propsys->getNumberOfUsedThrusters();
  };

  double CSC::getPowerOfUsedThruster(void)const {
    return m_propsys->getPowerOfUsedThruster();
  };

  double CSC::getPower(void)const { return m_propsys->getPower(); };

  // implemented in class CNSTAR, CT6, and RIT
  double CSC::getThrust(void)const
  { return m_propsys->getThrust(); };


  // implemented in class CNSTAR, CT6, and RIT
  double CSC::getMaxThrust(void)const
  { return m_propsys->getMaxThrust(); };


  // implemented in class CNSTAR, CT6, and RIT
  double CSC::getMassFlow(void)const
  { return m_propsys->getMassFlow(); };

}; // end class CSC

//------------------------------------------------------------------------------
inline void CSC::performSteering(CNC* nc, const CState& target_state)
{ 
  if (getCenterBodyId()!=SUN && isSolarSail()==true)
  {
    m_attitude.setSunVector(getRefSolSys()->Cart(getId(), SUN, MJD(), AUD));
  }
  m_steering.performSteering(this, nc, *this, target_state); 
} // void performSteering(CNC* nc, const CState& target_state)

//------------------------------------------------------------------------------
inline void CSC::setAttitude(const double alpha, const double beta)
{ 
  m_attitude.setAttitude(alpha, beta); 
} // void setAttitude(const double alpha, const double beta)

//------------------------------------------------------------------------------
inline void CSC::setAttitude(const double xi_r,const double xi_t
                            ,const double xi_h)
{
  m_attitude.setAttitude(xi_r, xi_t, xi_h);
} // void setAttitude(const double xi_r,const double xi_t,const double xi_h)

//------------------------------------------------------------------------------
inline int CSC::getCtrlInputDim(void) const
{ 
#if defined development_version
  assert (m_input_ncp.size()== m_n_ctrl_input_dim);
#endif
  return m_input_ncp.size();
} // int getCtrlInputDim(void) const

//------------------------------------------------------------------------------
inline void CSC::setCtrlInputDim(const int input_dim)
{
  if (getCtrlInputDim() != input_dim)
  {
    m_input_ncp.Set(MathLib::Vector(m_n_ctrl_input_dim = input_dim));
  }
} // void setCtrlInputDim(const int input_dim)

//------------------------------------------------------------------------------
inline double CSC::setCtrlInput(const int n, const double val)
{
  m_input_ncp(n) = val;
} // double setCtrlInput(const int n, const double val)

//------------------------------------------------------------------------------
inline void CSC::setCtrlInputVec(const MathLib::Vector& vec)
{
  m_input_ncp = vec;
} // void setCtrlInputVec(const MathLib::Vector& vec)

//------------------------------------------------------------------------------
inline double CSC::getCtrlInput(const int n) const
{ 
  return m_input_ncp(n);
} // double getCtrlInput(const int n) const

//------------------------------------------------------------------------------
inline MathLib::Vector& CSC::getCtrlInputVec(void)
{
  return m_input_ncp;
} // inline CSC::MathLib::Vector& getCtrlInputVec(void);

//------------------------------------------------------------------------------
inline const MathLib::Vector& CSC::getCtrlInputVec(void) const
{
  return m_input_ncp;
} // const MathLib::Vector& getCtrlInputVec(void) const

//------------------------------------------------------------------------------
inline int CSC::getCtrlOutputDim(void) const 
{
#if defined development_version
  assert (m_nc_output.size()== m_n_ctrl_output_dim);
#endif
  return m_nc_output.size(); 
} // int getCtrlOutputDim(void) const 

//------------------------------------------------------------------------------
inline void CSC::setCtrlOutputDim(const int output_dim)
{
  if (getCtrlOutputDim() != output_dim)
  {
    m_nc_output.Set(MathLib::Vector(m_n_ctrl_output_dim = output_dim));
  }
} // void setCtrlOutputDim(const int output_dim)

//------------------------------------------------------------------------------
inline double CSC::getCtrlOutput(const int n) const 
{
  return m_nc_output(n);
} // double getCtrlOutput(const int n) const 

//------------------------------------------------------------------------------
inline const MathLib::Vector& CSC::getCtrlOutputVec(void) const
{
  return m_nc_output;
} // MathLib::Vector& getCtrlOutputVector(void) const

//------------------------------------------------------------------------------
inline MathLib::Vector& CSC::getCtrlOutputVec(void)
{
  return m_nc_output;
} // MathLib::Vector& getCtrlOutputVector(void) const

//------------------------------------------------------------------------------
inline void CSC::setCtrlOutput(const int n, const double val)
{
  m_nc_output(n)=val;
} // void setCtrlOutput(const int n, const double val)

//------------------------------------------------------------------------------
inline void CSC::CheckRemainingPropellant(const double dt)
{
  m_propsys->CheckRemainingPropellant(dt);
} // void CheckRemainingPropellant(const double dt)

//------------------------------------------------------------------------------
inline void CSC::printVrmlAppearance(std::ostream& os) const throw ()
{
  m_propsys->printVrmlAppearance(os, getName());
} // void printVrmlAppearance(std::ostream& os) const throw ();
#endif // CSC_H_INCLUDED
