/*! \file CSC.cpp \brief A class that provides the functionality
for spacecraft.

\version

1.7

\author

Bernd Dachwald\n German Aerospace Center (DLR)\n
Institute of Space Simulation\n Linder Hoehe\n 51147 Cologne\n Germany\n
bernd.dachwald@dlr.de

\author

  Andreas Ohndorf\n
  German Aerospace Center (DLR)\n
  Institute of Space Simulation\n
  Linder Hoehe\n
  51147 Cologne\n
  Germany\n
  andreas.ohndorf@gmx.de

\date

24.10.2006 (AO)

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
warranty of merchantibility or fitness for a particular purpose. */

// STL
#include <sstream>
#include <new>

// SCLib
#include <configSC.h>
#include <CSC.h>

// CNCLib
#include <CNC.h>


using namespace std;

unsigned long int CSC::m_count = 0;


//------------------------------------------------------------------------------
void CSC::setPropSys(CPropSys* new_prop_sys)
{
  SafeDelete(m_propsys);
  m_propsys = new_prop_sys;

  if (new_prop_sys == 0)
  {
    ostringstream os;
    os << __FILE__ << ", L" << __LINE__ << ": CSC::setPropulsionSystem(";
    os << "CPropSys*): pointer of new value is zero!\n";
    Warning(os.str());
  } // if (new_prop_sys == 0)
} // void setPropSys(CPropSys* new_prop_sys)


//------------------------------------------------------------------------------
void CSC::setPropSys(const PropulsionSystemType prop_sys_type) 
{
  SafeDelete(m_propsys);
  try 
  {
    switch (prop_sys_type)
    {
    case SOLARSAIL:m_propsys = new CSolarSail(this);break;
    case NEP:      m_propsys = new CNEP(this);break;
    case NSTAR:    m_propsys = new CNSTAR(this);break;
    case T6:       m_propsys = new CT6(this);break;
    case RIT:      m_propsys = new CRIT(this);break;
    default:       m_propsys = 0;break;
    } // switch (prop_sys_type)
  } // try 
  catch (std::bad_alloc& ) // exception handler if new throws exception
  {
    ostringstream os;
    os << __FILE__ << ", L" << __LINE__ << ": CSC::setPropSys(";
    os << "const PropulsionSystemType): bad_alloc excp. for propulsion\n";
    os << " system type:" << prop_sys_type << endl;
    Warning(os.str());
  } // catch (std::bad_alloc& )
} // void setPropSys(const PropulsionSystemType prop_sys_type) 


//------------------------------------------------------------------------------
void CSC::UpdateThrustAcc(void)
{
  double acc;

  if (isSolarSail())
    acc=(getLightness()*GM_SUN_SI/(CState::r()*AU*CState::r()*AU))*pow(cos(getConeAngle()), 2);
  else if (isEP())
    acc=(!isEmpty())?getThrust()/getMass():0.0;
  else
    FatalError("CSC::UpdateThrustAcc(): execution not possible for undefined propulsion system");

  m_thrustacc(0)=acc*cos(getConeAngle()); // a_r
  m_thrustacc(1)=acc*cos(getClockAngle())*sin(getConeAngle()); // a_t
  m_thrustacc(2)=acc*sin(getClockAngle())*sin(getConeAngle()); // a_n
} // void UpdateThrustAcc(void)


//------------------------------------------------------------------------------
void CSC::UpdatedElements(void)
{
  const double n=CState::MeanMotion();
  const double b=CState::SemiMinorAxis();

  const double a_r    =m_thrustacc(0);
  const double a_phi  =m_thrustacc(1);
  const double a_theta=m_thrustacc(2);

  const double sinf  =sin(CState::f());
  const double cosf  =cos(CState::f());
  const double sinphi=sin(CState::phi());
  const double cosphi=cos(CState::phi());
  const double sini  =sin(CState::i());
  const double cosi  =cos(CState::i());

  const double h=CState::Momentum();
  const double p=CState::SemilatusRectum();

  const double pr=p+CState::r();

  const double da=2.0*CState::a()*CState::a()*((CState::e()*sinf*a_r)+(p/CState::r()*a_phi))/h;

  const double de=((p*sinf*a_r)+((pr*cosf+CState::r()*CState::e())*a_phi))/h;

  const double di=CState::r()*cosphi*a_theta/h;

  double dOmega;
  if (sini!=0.0)
    dOmega=CState::r()*sinphi*a_theta/(h*sini); else
    dOmega=0.0;

  double domega;
  if (CState::e()!=0.0&&sini!=0.0)
    domega=((-p*cosf*a_r)+(pr*sinf*a_phi))/(h*CState::e())-CState::r()*sinphi*cosi*
     a_theta/(h*sini); 
  else
  {
    if (CState::e()!=0.0)
      domega=((-p*cosf*a_r)+(pr*sinf*a_phi))/(h*CState::e()); else
      domega=0.0;
  }

  double dM;
  if (CState::e()!=0.0)
    dM=n+b*((p*cosf-2.0*CState::r()*CState::e())*a_r-pr*sinf*a_phi)/(CState::a()*h*CState::e()); else
    dM=n;

  m_dkepler(0) = da;
  m_dkepler(1) = de;
  m_dkepler(2) = di;
  m_dkepler(3) = dOmega;
  m_dkepler(4) = domega;
  m_dkepler(5) = dM;
} // void UpdatedElements(void)


//------------------------------------------------------------------------------
void CSC::UpdateDeltaV(const double delta_t)
{
  m_delta_v(0)+=delta_t*DAY2SEC*fabs(m_thrustacc(0));
  m_delta_v(1)+=delta_t*DAY2SEC*fabs(m_thrustacc(1));
  m_delta_v(2)+=delta_t*DAY2SEC*fabs(m_thrustacc(2));
} // void UpdateDeltaV(const double delta_t)


//------------------------------------------------------------------------------
CSC::CSC(SpaceLib::CSolarSystem* sol_sys) :
CBody(),
m_scmass(0.0),
m_moi(3),
m_msm(3),
m_thrustacc(3),
m_delta_v(3),
m_dkepler(6),
m_n_ctrl_input_dim(55),
m_input_ncp(55),
m_n_ctrl_output_dim(13),
m_nc_output(13),
m_DEcalls(0),
m_propsys(NULL) 
{
  // store pointer of reference CSolarSystem object and center state
  if (sol_sys != 0)
  {
    setRefSolSys(sol_sys);           // our s/c state is sun centered by default
    setCenterBodyId(SUN);
    // after storing of sun's ID we also have to store the pointer to the
    // particular object
    setCenterState(&(sol_sys->State(getCenterBodyId(), MJD_J2000)));
  } // if (sol_sys != 0)
  m_attitude.m_sc=this;
  m_count++;
} // CSC(SpaceLib::CSolarSystem* sol_sys)

//------------------------------------------------------------------------------
CSC::CSC(
 const CInputSC&isc,
 SpaceLib::CSolarSystem* solSys) :
  CBody(),
  m_scmass(isc.sc_mass),
  m_moi(isc.moments_of_inertia),
  m_msm(isc.max_steering_moments),
  m_thrustacc(3),
  m_delta_v(3),
  m_dkepler(6),
  m_n_ctrl_input_dim(55),
  m_input_ncp(55),
  m_n_ctrl_output_dim(13),
  m_nc_output(13),
  m_DEcalls(0),
  m_propsys(0)
{

  CBody::setName(isc.name);

   if (solSys!=0)   // store pointer of reference CSolarSystem object
   {
     setRefSolSys(solSys);     
     setCenterBodyId(SUN);  // the state is sun centered by default
     setCenterState(&(solSys->State(SUN, MJD_J2000)));
   } // if (solSys!=0)

   setPropSys(isc.psystype);
   m_attitude.m_sc=this;

   if (isSolarSail())
   {
     setLightness(isc.sailcraft_performance);
     setOpticalModel(
      isc.ideal_or_real_sail,
      0.0, // srd
      isc.max_temperature,
      isc.reflection_coefficient_new,
      isc.specular_reflection_coefficient_new,
      isc.emission_coefficient_front_new,
      isc.emission_coefficient_back_new,
      isc.reflection_coefficient_aged,
      isc.specular_reflection_coefficient_aged,
      isc.emission_coefficient_front_aged,
      isc.emission_coefficient_back_aged,
      isc.reflection_coefficient_halflife,
      isc.specular_reflection_coefficient_halflife,
      isc.emission_coefficient_front_halflife,
      isc.emission_coefficient_back_halflife,
      isc.non_lambertian_coefficient_front,
      isc.non_lambertian_coefficient_back);
   } else if (isNEP()) {
     setSCMass(isc.sc_mass);
     setInitialPropellantMassInterval(isc.min_mass_initial_propellant,
      isc.max_mass_initial_propellant);
     setMaxThrust(isc.thrust_max);
     setIsp(isc.isp);
     setThrottleType(isc.throttletype);
   } else if (isSEP()) 
   {
     setSCMass(isc.sc_mass);
     setInitialPropellantMassInterval(isc.min_mass_initial_propellant,
      isc.max_mass_initial_propellant);
     setCharacteristicPower(isc.char_power);
     setPowerExponent(isc.p_exponent);
     setNumberOfThrusters(isc.thrusters);
     setThrottleType(isc.throttletype);
   } 
   else
     FatalError("CSC::CSC(): unknown spacecraft type");

     m_count++;
} // CSC(const CInputSC&isc,SpaceLib::CSolarSystem* solSys)


//------------------------------------------------------------------------------
//! Attention: Never think of removing CState- and CBody copy constructor
//             calls from initializer list, because otherwise only the
//             standard constructor of these base classes will be used,
//             which obviously is wrong.
 CSC::CSC(const CSC& sc)
   :CState(sc)
   ,CBody(sc)   
   ,m_scmass(sc.m_scmass)
   ,m_moi(sc.m_moi)
   ,m_msm(sc.m_msm)
   ,m_thrustacc(sc.m_thrustacc)
   ,m_delta_v(sc.m_delta_v)
   ,m_dkepler(sc.m_dkepler)
   ,m_n_ctrl_input_dim(sc.m_n_ctrl_input_dim)
   ,m_input_ncp(sc.m_input_ncp)
   ,m_n_ctrl_output_dim(sc.m_n_ctrl_output_dim)
   ,m_nc_output(sc.m_nc_output)
   ,m_DEcalls(sc.m_DEcalls)
   ,m_attitude(sc.m_attitude)
   ,m_steering(sc.m_steering)   
 {
   m_propsys = (sc.m_propsys!=0?sc.m_propsys->clone(this):0);
   CSolarSystem* solar_system = sc.getRefSolSys();
   if (solar_system != 0) 
   {
     setRefSolSys(solar_system);        
     setCenterBodyId(SUN); // the state is sun centered by default
     setCenterState(&(solar_system->State(getCenterBodyId(), MJD_J2000)));
   }
   m_count++;
 } // CSC(const CSC& sc)

//------------------------------------------------------------------------------
CSC::CSC(
 const PropulsionSystemType psystype,
 SpaceLib::CSolarSystem* solSys) 
 :
  CBody(),
  m_scmass(0.0),
  m_moi(3),
  m_msm(3),
  m_thrustacc(3),
  m_delta_v(3),
  m_dkepler(6),
  m_n_ctrl_input_dim(55),
  m_input_ncp(55),
  m_n_ctrl_output_dim(13),
  m_nc_output(13),
  m_DEcalls(0),
  m_propsys(0)
  {   
   if (solSys!=0) // store pointer of reference CSolarSystem object
   {
     setRefSolSys(solSys);    
     setCenterBodyId(SUN);// the state is sun centered by default
     setCenterState(&(solSys->State(SUN, MJD_J2000)));
   }
   setPropSys(psystype);
   m_attitude.m_sc=this;
   m_count++;
} // CSC::CSC(const PropulsionSystemType psystype,SpaceLib::CSolarSystem* solSys)

//------------------------------------------------------------------------------
void CSC::UpdateAccDeltaVec(const double delta_t)
{
  UpdateThrustAcc();
  UpdatedElements();
  UpdateDeltaV(delta_t);
} // void UpdateAccDeltaVec(const double delta_t)


//------------------------------------------------------------------------------
Vector CSC::PolarStatePropellantMassSRDVector(void) const
throw (std::logic_error)
{
  const std::string method_name =
    "Vector CSC::PolarStatePropellantMassSRDVector(void)const";
  const std ::string err_msg_01 = 
    ": not implemented state vector initialization for given coordinate type";

  switch (getPropSys()->getEomCoordinateType())
  {
  case POLAR:
    return Vector(CState::r(),CState::phi(),CState::theta()
      ,CState::dr(),CState::dphi(),CState::dtheta()
      ,getPropellantMass(),getSolarRadiationDose());
    break;
  case CARTESIAN:
    return Vector(Stack(CState::Cart()
      ,Vector(getPropellantMass(), getSolarRadiationDose())));
    break;
  default:
    throw std::logic_error(method_name + err_msg_01);
    break;
  } // switch (getPropSys()->getEomCoordinateType())
} // Vector CSC::PolarStatePropellantMassSRDVector(void) const


//------------------------------------------------------------------------------
Vector CSC::ThrustAcc_xyz(void) const
{
  const double sp=sin(CState::phi());
  const double cp=cos(CState::phi());
  const double st=sin(CState::theta());
  const double ct=cos(CState::theta());

  const double zeta=atan2(CState::dtheta(), CState::dphi()*ct);
  const double sz  =sin(zeta);
  const double cz  =cos(zeta);

  Matrix A(3, 3);
  A.SetRow(0, Vector(cp*ct, sp*ct, st));
  A.SetRow(1, Vector(-sp*cz-cp*st*sz, cp*cz-sp*st*sz, ct*sz));
  A.SetRow(2, Vector(sp*sz-cp*st*cz, -cp*sz-sp*st*cz, ct*cz));
  Matrix InvA=Inv(A);
  Vector axyz=InvA*m_thrustacc;
  return axyz;
} // Vector ThrustAcc_xyz(void) const


//------------------------------------------------------------------------------
void CSC::setPropellantMassSRDVector(const Vector& v, const double t)
{
#ifdef development_version
  if (v.getDimension()!=8)
    FatalError("CSC::setPolarStatePropellantMassSRDVector(): vector has wrong dimension");
#endif

  switch (getPropSys()->getEomCoordinateType())
  {
  case POLAR:
    CState::SetState(v(0),v(1),v(2),v(3),v(4),v(5)
      , POLAR, t, getCenterBodyId());
    setPropellantMass(v(6));
    setSolarRadiationDose(v(7));
    break;
  case CARTESIAN:
    CState::SetState(v(0),v(1),v(2),v(3),v(4),v(5)
      , CARTESIAN, t, getCenterBodyId());
    setPropellantMass(v(6));
    setSolarRadiationDose(v(7));
    break;
  default:
    assert(1==0); // intentionally bug to signal error (debugging purpose)
    break;
  }
} // void setPropellantMassSRDVector(const Vector& v, const double t)


//------------------------------------------------------------------------------
void CSC::setSteeringStrategy(const Vector&lsl_vector)
{
  switch (lsl_vector.getDimension())
  {
    case 2:
      setSteeringStrategy(ALPHABETA);
      break;
    case 5:
      setSteeringStrategy(DEOPT);
      break;
    case 10:
      setSteeringStrategy(DEMINMAX);
      break;
    default:
      FatalError("CSC::setSteeringStrategy(): local steering vector not suitable for spacecraft steering");
  } // switch (lsl_vector.getDimension())
} // void CSC::setSteeringStrategy(const Vector&lsl_vector)


//------------------------------------------------------------------------------
void CSC::setSteeringStrategy(CNC*nc)
{
  int nc_output_neurons=nc->getNumberOfOutputNeurons();

  if (isSolarSail()) {
    if (GravityAssistSteering()) {
      switch (nc_output_neurons) {
        case 5:
          setSteeringStrategy(DIRECTION);
          break;
        case 7:
          setSteeringStrategy(DEOPT);
          break;
        case 12:
          setSteeringStrategy(DEMINMAX);
          break;
        default:
          FatalError("CSC::setSteeringStrategy(): neurocontroller topology not suitable for sailcraft steering");
      }
    } else // if (!GravityAssistSteering())
    {
      switch (nc_output_neurons) {
        case 3:
          setSteeringStrategy(DIRECTION);
          break;
        case 5:
          setSteeringStrategy(DEOPT);
          break;
        case 10:
          setSteeringStrategy(DEMINMAX);
          break;
        default:
          FatalError("CSC::setSteeringStrategy(): neurocontroller topology not suitable for sailcraft steering");
      }
    }
  } else {
    if (GravityAssistSteering()) {
      switch (nc_output_neurons) {
        case 6:
          setSteeringStrategy(DIRECTION);
          break;
        case 8:
          setSteeringStrategy(DEOPT);
          break;
        case 13:
          setSteeringStrategy(DEMINMAX);
          break;
        default:
          FatalError("CSC::setSteeringStrategy(): neurocontroller topology not suitable for spacecraft steering");
      }
    } else // if (!GravityAssistSteering())
    {
      switch (nc_output_neurons) {
        case 5:
          setSteeringStrategy(DIRECTION);
          break;
        case 6:
          setSteeringStrategy(DEOPT);
          break;
        case 11:
          setSteeringStrategy(DEMINMAX);
          break;
        default:
          FatalError("CSC::setSteeringStrategy(): neurocontroller topology not suitable for spacecraft steering");
      }
    }
  }
} // void setSteeringStrategy(CNC*nc)


//------------------------------------------------------------------------------
double CSC::MaxThrustAcceleration(void) const
{
  double maxacc;

  if (isSolarSail())
  {
    maxacc=m_propsys->getLightness()*GM_SUN_AUD*AU/(DAY2SEC*DAY2SEC);
  } // if (isSolarSail())
  else if (isSEP()||isNEP()) 
  { 
    maxacc=getMaxThrust()/getSCMass();
  }  // else if (isSEP()||isNEP())
  else
    FatalError("CSC::MaxThrustAcceleration(): unknown spacecraft type");

  return maxacc;
} // double MaxThrustAcceleration(void) const


//------------------------------------------------------------------------------
bool CSC::CheckSolarDistance(const double min_solar_distance) const
{
  const bool r_too_small=(CState::r()<min_solar_distance)?true:false;
  const bool a_too_small=((CState::a()<min_solar_distance)&&(CState::a()>0.0))?true:false;
  return (r_too_small || a_too_small)?false:true;
} // bool CheckSolarDistance(const double min_solar_distance) const


//------------------------------------------------------------------------------
bool CSC::CheckTemperature(void) const
{
  bool result=true;
  if (isSolarSail()) 
  {
    if (getTemperature() > m_propsys->getMaxTemperature())
      result=false;
  }
  return result;
} // bool CheckTemperature(void) const


//------------------------------------------------------------------------------
void CSC::Init(const double srd_or_mp /* = 0.0 */)
{
  if (isSolarSail())
  {
    m_propsys->Init(srd_or_mp);
  }
  else if (isEP())
  {
    m_propsys->setInitialPropellantMass(srd_or_mp);
    m_propsys->Refuel();
  }
  else
    FatalError("CSC::Init(): unknown spacecraft type");
} // void Init(const double srd_or_mp /* = 0.0 */)


//------------------------------------------------------------------------------
void CSC::setLaunchState(
 const double launch_time,
 const BodyID launch_body,
 const double launch_velocity,
 const double launch_azimuth,
 const double launch_elevation)
{
   CSolarSystem* sol_sys(getRefSolSys());
   if (sol_sys != 0)                // check for existing reference solar system
   {
     const SpaceLib::UnitType ut = 
       sol_sys->State(launch_body,launch_time).Unit();

     const double dcf(ut == SI ? M2AU : 1.0);    // distance conversion factor
     const double vcf(ut == SI ? MS2AUD : 1.0);  // velocity conversion factor
     const double r(sol_sys->State(launch_body, launch_time).r());
     const double phi(sol_sys->State(launch_body, launch_time).phi());
     const double theta(sol_sys->State(launch_body, launch_time).theta());
     const double vr(sol_sys->State(launch_body, launch_time).vr());
     const double vphi(sol_sys->State(launch_body, launch_time).vphi());
     const double vtheta(sol_sys->State(launch_body, launch_time).vtheta());

     // radial and tangential velocities
     const double vrg = 
       vcf*vr+launch_velocity*cos(launch_elevation)*sin(launch_azimuth);
     const double vphig = 
       vcf*vphi+launch_velocity*cos(launch_elevation)*cos(launch_azimuth);
     const double vthetag = 
       vcf*vtheta +launch_velocity*sin(launch_elevation);

     const double drg     = vrg;
     const double dphig   = vphig/(r*dcf)/cos(theta);
     const double dthetag = vthetag/(r*dcf);

     CState::SetState(r*dcf, phi, theta, drg, dphig, dthetag, POLAR
       , launch_time, getCenterBodyId());
   }
   else Warning("CSC::setLaunchState(.....): CSolarSystem pointer is null!\n");
} // void setLaunchState(...)


//------------------------------------------------------------------------------
void CSC::setLaunchState(
  const double launch_time
 ,const double a
 ,const double e
 ,const double i
 ,const double Omega
 ,const double omega
 ,const double M
 ,const double launch_velocity
 ,const double launch_azimuth
 ,const double launch_elevation) 
{
   CState::SetState(a, e, i, Omega, omega, M, KEPLERM
     , launch_time, getCenterBodyId());

   const double vrg=vr()+launch_velocity*cos(launch_elevation)*sin(launch_azimuth);
   const double vphig=vphi()+launch_velocity*cos(launch_elevation)*cos(launch_azimuth);
   const double vthetag=vtheta()+launch_velocity*sin(launch_elevation);

   const double drg    =vrg;
   const double dphig  =vphig/CState::r()/cos(CState::theta());
   const double dthetag=vthetag/CState::r();

   CState::SetState(CState::r(), CState::phi(), CState::theta(),
    drg, dphig, dthetag, POLAR, launch_time, getCenterBodyId());
} // void setLaunchState(...)


//------------------------------------------------------------------------------
void CSC::PrintHeader(ofstream&out, const char*comment /* = NULL */) const
{
  CState::PrintHeader(out, comment);

  out<<", ";
  if (!comment)
  {
    if (isNEP()||isSEP())
    {
      out<<"throttle, ";
      if (isSEP())
      {
        out<<"P [kW], ";
        if (getNumberOfThrusters()>1.0)
        {
          out<<"ni, ";
          out<<"Pi [kW], ";
        } // if (getNumberOfThrusters()>1.0)
      } // if (isSEP())
      out<<"F [N], ";
      out<<"m [kg], ";
      out<<"m_p [kg], ";
      out<<"dm_p [kg/d], ";
      if (isSEP())
        out<<"Isp [s], ";
    } // if (isNEP()||isSEP())
    out<<"alpha-+180 [rad], ";
    out<<"alpha0-360 [rad], ";
    out<<"beta [rad], ";
    if (isSolarSail() && isRealSail())
    {
      out<<"T [K], ";
      out<<"Sigma [S0*d], ";
      out<<"rho, ";
      out<<"s, ";
      out<<"epsilonf, ";
      out<<"epsilonb, ";
      out<<"Bf, ";
      out<<"Bb, ";
    } // if (isSolarSail() && isRealSail())
    out<<"a [mm/s^2], ";
    out<<"a_r [mm/s^2], ";
    out<<"a_t [mm/s^2], ";
    out<<"a_h [mm/s^2], ";
    out<<"a_th [mm/s^2], ";
    out<<"DeltaV [m/s], ";
    out<<"DeltaV_r [m/s], ";
    out<<"DeltaV_t [m/s], ";
    out<<"DeltaV_h [m/s], ";
    out<<"DeltaV_th [m/s], ";
    out<<"da [AU/d], ";
    out<<"de [1/d], ";
    out<<"di [rad/d], ";
    out<<"dOmega [rad/d], ";
    out<<"domega [rad/d], ";
    out<<"dM [rad/d]";

    if (getSteeringStrategy()==DEMINMAX)
    {
      out<<", ";
      out<<"a_max, ";
      out<<"a_min, ";
      out<<"e_max, ";
      out<<"e_min, ";
      out<<"i_max, ";
      out<<"i_min, ";
      out<<"Omega_max, ";
      out<<"Omega_min, ";
      out<<"omega_max, ";
      out<<"omega_min";
    } // if (getSteeringStrategy()==DEMINMAX)
    else if (getSteeringStrategy()==DEOPT)
    {
      out<<", ";
      out<<"a_opt, ";
      out<<"e_opt, ";
      out<<"i_opt, ";
      out<<"Omega_opt, ";
      out<<"omega_opt";
    } // else if (getSteeringStrategy()==DEOPT)
  } //if (!comment)
  else 
  {
    if (isNEP()||isSEP())
    {
      out<<"throttle (" <<comment<<"), ";
      if (isSEP())
      {
        out<<"P (" <<comment<<") [kW], ";
        if (getNumberOfThrusters()>1.0)
        {
          out<<"ni (" <<comment<<"), ";
          out<<"Pi (" <<comment<<") [kW], ";
        } // if (getNumberOfThrusters()>1.0)
      } // if (isSEP())
      out<<"F (" <<comment<<") [N], ";
      out<<"m (" <<comment<<") [kg], ";
      out<<"m_p (" <<comment<<") [kg], ";
      out<<"dm_p (" <<comment<<") [kg/d], ";
      if (isSEP())
        out<<"Isp (" <<comment<<") [s], ";
    } // if (isNEP()||isSEP())
    out<<"alpha-+180 (" <<comment<<") [rad], ";
    out<<"alpha0-360 (" <<comment<<") [rad], ";
    out<<"beta (" <<comment<<") [rad], ";
    if (isSolarSail() && isRealSail())
    {
      out<<"T (" <<comment<<") [K], ";
      out<<"Sigma (" <<comment<<") [S0*d], ";
      out<<"rho (" <<comment<<"), ";
      out<<"s (" <<comment<<"), ";
      out<<"epsilonf (" <<comment<<"), ";
      out<<"epsilonb (" <<comment<<"), ";
      out<<"Bf (" <<comment<<"), ";
      out<<"Bb (" <<comment<<"), ";
    } // if (isSolarSail() && isRealSail())
    out<<"a (" <<comment<<") [mm/s^2], ";
    out<<"a_r (" <<comment<<") [mm/s^2], ";
    out<<"a_t (" <<comment<<") [mm/s^2], ";
    out<<"a_h (" <<comment<<") [mm/s^2], ";
    out<<"a_th (" <<comment<<") [mm/s^2], ";
    out<<"DeltaV (" <<comment<<") [m/s], ";
    out<<"DeltaV_r (" <<comment<<") [m/s], ";
    out<<"DeltaV_t (" <<comment<<") [m/s], ";
    out<<"DeltaV_h (" <<comment<<") [m/s], ";
    out<<"DeltaV_th (" <<comment<<") [m/s], ";
    out<<"da (" <<comment<<") [AU/d], ";
    out<<"de (" <<comment<<") [1/d], ";
    out<<"di (" <<comment<<") [rad/d], ";
    out<<"dOmega (" <<comment<<") [rad/d], ";
    out<<"domega (" <<comment<<") [rad/d], ";
    out<<"dM (" <<comment<<") [rad/d]";

    if (getSteeringStrategy()==DEMINMAX)
    {
      out<<", ";
      out<<"a_max (" <<comment<<"), ";
      out<<"a_min (" <<comment<<"), ";
      out<<"e_max (" <<comment<<"), ";
      out<<"e_min (" <<comment<<"), ";
      out<<"i_max (" <<comment<<"), ";
      out<<"i_min (" <<comment<<"), ";
      out<<"Omega_max (" <<comment<<"), ";
      out<<"Omega_min (" <<comment<<"), ";
      out<<"omega_max (" <<comment<<"), ";
      out<<"omega_min (" <<comment<<")";
    } // if (getSteeringStrategy()==DEMINMAX)
    else if (getSteeringStrategy()==DEOPT)
    {
      out<<", ";
      out<<"a_opt (" <<comment<<"), ";
      out<<"e_opt (" <<comment<<"), ";
      out<<"i_opt (" <<comment<<"), ";
      out<<"Omega_opt (" <<comment<<"), ";
      out<<"omega_opt (" <<comment<<")";
    } // else if (getSteeringStrategy()==DEOPT)
  } // if(!comment) else branch
} // void PrintHeader(ofstream&out, const char*comment /* = NULL */) const


//------------------------------------------------------------------------------
void CSC::PrintDeltaHeader(ofstream&out, const char*comment /* = NULL */) const
{
  if (!comment) 
  {
    out<<"D x " <<"[AU], ";
    out<<"D y " <<"[AU], ";
    out<<"D z " <<"[AU], ";
    out<<"D dx " <<"[AU/d], ";
    out<<"D dy " <<"[AU/d], ";
    out<<"D dz " <<"[AU/d], ";
    out<<"D r " <<"[AU], ";
    out<<"D phi " <<"[rad], ";
    out<<"D theta " <<"[rad], ";
    out<<"D dr " <<"[AU/d], ";
    out<<"D dphi " <<"[rad/d], ";
    out<<"D dtheta " <<"[rad/d], ";
    out<<"D a " <<"[AU], ";
    out<<"D e, ";
    out<<"D i " <<"[rad], ";
    out<<"D Omega " <<"[rad], ";
    out<<"D omega " <<"[rad], ";
    out<<"D curlypi " <<"[rad], ";
    out<<"D M " <<"[rad], ";
    out<<"D E " <<"[rad], ";
    out<<"D f " <<"[rad], ";
    out<<"D L " <<"[rad], ";
    out<<"D v " <<"[AU/d], ";
    out<<"D Ekin " <<"[AU^2/d^2], ";
    out<<"D Epot " <<"[AU^2/d^2], ";
    out<<"D Etot " <<"[AU^2/d^2]";
  } // if (!comment) 
  else 
  {
    out<<"D x (" <<comment<<") [AU], ";
    out<<"D y (" <<comment<<") [AU], ";
    out<<"D z (" <<comment<<") [AU], ";
    out<<"D dx (" <<comment<<") [AU/d], ";
    out<<"D dy (" <<comment<<") [AU/d], ";
    out<<"D dz (" <<comment<<") [AU/d], ";
    out<<"D r (" <<comment<<") [AU], ";
    out<<"D phi (" <<comment<<") [rad], ";
    out<<"D theta (" <<comment<<") [rad], ";
    out<<"D dr (" <<comment<<") [AU/d], ";
    out<<"D dphi (" <<comment<<") [rad/d], ";
    out<<"D dtheta (" <<comment<<") [rad/d], ";
    out<<"D a (" <<comment<<") [AU], ";
    out<<"D e (" <<comment<<"), ";
    out<<"D i (" <<comment<<") [rad], ";
    out<<"D Omega (" <<comment<<") [rad], ";
    out<<"D omega (" <<comment<<") [rad], ";
    out<<"D curlypi (" <<comment<<") [rad], ";
    out<<"D M (" <<comment<<") [rad], ";
    out<<"D E (" <<comment<<") [rad], ";
    out<<"D f (" <<comment<<") [rad], ";
    out<<"D L (" <<comment<<") [rad], ";
    out<<"D v (" <<comment<<") [AU/d], ";
    out<<"D Ekin (" <<comment<<") [AU^2/d^2], ";
    out<<"D Epot (" <<comment<<") [AU^2/d^2], ";
    out<<"D Etot (" <<comment<<") [AU^2/d^2]";
  } // if (!comment) else branch
} // void PrintDeltaHeader(ofstream&out, const char*comment /* = NULL */) const


//------------------------------------------------------------------------------
void CSC::PrintNCInputHeader(ofstream&out) const
{
  for (int i=0; i<getCtrlInputDim(); ++i)
    out<<", NC input (" <<i<<")";
} // void PrintNCInputHeader(ofstream&out) const


//------------------------------------------------------------------------------
void CSC::PrintNCOutputHeader(ofstream&out) const
{
  for (int i=0; i<getCtrlOutputDim(); ++i)
    out<<", NC output (" <<i<<")";
} // void PrintNCOutputHeader(ofstream&out) const


//------------------------------------------------------------------------------
void CSC::PrintProximityHeader(ofstream&out, const char*comment) const
{
  if (!comment)
  {
    out<<"Dist " <<"[km], ";
    out<<"Vel " <<"[km/s]";
  } // if (!comment)
  else 
  {
    out<<"Dist (" <<comment<<") [km], ";
    out<<"Vel (" <<comment<<") [km/s]";
  } // if (!comment) else branch
} // void PrintProximityHeader(ofstream&out, const char*comment) const


//------------------------------------------------------------------------------
void CSC::Print(ofstream&out) const
{
  CState::Print(out);

  out<<", ";
  if (isNEP()||isSEP())
  {
    out<<getThrottle()<<", ";
    if (isSEP())
    {
      out<<getPower()<<", ";
      if (getNumberOfThrusters()>1.0)
      {
        out<<getNumberOfUsedThrusters()<<", ";
        out<<getPowerOfUsedThruster()<<", ";
      } // if (getNumberOfThrusters()>1.0)
    } // if (isSEP())
    out<<getThrust()<<", ";
    out<<getMass()<<", ";
    out<<getPropellantMass()<<", ";
    out<<getMassFlow()*MASSFLOW_SI2AUD<<", ";
    if (isSEP())
      out<<getIsp()<<", ";
  } // if (isNEP()||isSEP())
  double alpha1=getClockAngle();
  out<<alpha1<<", ";
  double alpha2=(alpha1>0.0)?alpha1:alpha1+TWOPI;
  out<<alpha2<<", ";
  out<<getConeAngle()<<", ";
  if (isSolarSail() && isRealSail())
  {
    out<<getTemperature()<<", ";
    out<<getSolarRadiationDose()<<", ";
    out<<getrho()<<", ";
    out<<gets()<<", ";
    out<<getepsilonf()<<", ";
    out<<getepsilonb()<<", ";
    out<<getBf()<<", ";
    out<<getBb()<<", ";
  } // if (isSolarSail() && isRealSail())
  out<<ThrustAcc()*1000.0<<", ";
  out<<ThrustAcc_r()*1000.0<<", ";
  out<<ThrustAcc_t()*1000.0<<", ";
  out<<ThrustAcc_h()*1000.0<<", ";
  out<<ThrustAcc_th()*1000.0<<", ";
  out<<DeltaV()<<", ";
  out<<DeltaV_r()<<", ";
  out<<DeltaV_t()<<", ";
  out<<DeltaV_h()<<", ";
  out<<DeltaV_th()<<", ";
  out<<dElements();

  if (getSteeringStrategy()==DEMINMAX)
  {
    out<<", ";
    out<<m_steering.m_lslvec(amax)<<", ";
    out<<m_steering.m_lslvec(amin)<<", ";
    out<<m_steering.m_lslvec(emax)<<", ";
    out<<m_steering.m_lslvec(emin)<<", ";
    out<<m_steering.m_lslvec(imax)<<", ";
    out<<m_steering.m_lslvec(imin)<<", ";
    out<<m_steering.m_lslvec(Omegamax)<<", ";
    out<<m_steering.m_lslvec(Omegamin)<<", ";
    out<<m_steering.m_lslvec(omegamax)<<", ";
    out<<m_steering.m_lslvec(omegamin);
  } // if (getSteeringStrategy()==DEMINMAX)
  if (getSteeringStrategy()==DEOPT)
  {
    out<<", ";
    out<<m_steering.m_lslvec(aopt)<<", ";
    out<<m_steering.m_lslvec(eopt)<<", ";
    out<<m_steering.m_lslvec(iopt)<<", ";
    out<<m_steering.m_lslvec(Omegaopt)<<", ";
    out<<m_steering.m_lslvec(omegaopt);
  } // if (getSteeringStrategy()==DEOPT)
} // void Print(ofstream&out) const


//------------------------------------------------------------------------------
void CSC::PrintDelta(
  ofstream& out
 ,const int planet
 ,bool spacecraft_minus_planet /* = false */) const 
{
   if (!spacecraft_minus_planet) 
   {
     out<< getRefSolSys()->x(planet, MJD())-CState::x()<<", ";
     out<<getRefSolSys()->y(planet, MJD())-CState::y()<<", ";
     out<<getRefSolSys()->z(planet, MJD())-CState::z()<<", ";
     out<<getRefSolSys()->dx(planet, MJD())-CState::dx()<<", ";
     out<<getRefSolSys()->dy(planet, MJD())-CState::dy()<<", ";
     out<<getRefSolSys()->dz(planet, MJD())-CState::dz()<<", ";
     out<<getRefSolSys()->r(planet, MJD())-CState::r()<<", ";
     out<<DifferenceAngle(getRefSolSys()->phi(planet, MJD()), CState::phi())<<", ";
     out<<getRefSolSys()->theta(planet, MJD())-CState::theta()<<", ";
     out<<getRefSolSys()->dr(planet, MJD())-CState::dr()<<", ";
     out<<getRefSolSys()->dphi(planet, MJD())-CState::dphi()<<", ";
     out<<getRefSolSys()->dtheta(planet, MJD())-CState::dtheta()<<", ";
     out<<getRefSolSys()->a(planet, MJD())-CState::a()<<", ";
     out<<getRefSolSys()->e(planet, MJD())-CState::e()<<", ";
     out<<getRefSolSys()->i(planet, MJD())-CState::i()<<", ";
     out<<DifferenceAngle(getRefSolSys()->Omega(planet, MJD()), CState::Omega())<<", ";
     out<<DifferenceAngle(getRefSolSys()->omega(planet, MJD()), CState::omega())<<", ";
     out<<DifferenceAngle(getRefSolSys()->curlypi(planet, MJD()), CState::curlypi())<<", ";
     out<<DifferenceAngle(getRefSolSys()->M(planet, MJD()), CState::M())<<", ";
     out<<DifferenceAngle(getRefSolSys()->E(planet, MJD()), CState::E())<<", ";
     out<<DifferenceAngle(getRefSolSys()->f(planet, MJD()), CState::f())<<", ";
     out<<DifferenceAngle(getRefSolSys()->L(planet, MJD()), CState::L())<<", ";
     out<<getRefSolSys()->Velocity(planet, BODY_STATE, 0.0, 0.0, 0.0, 0.0, 0.0,
      (*this))<<", ";
     out<<getRefSolSys()->KineticEnergy(planet, MJD())-CState::KineticEnergy()<<", ";
     out<<getRefSolSys()->PotentialEnergy(planet, MJD())-CState::PotentialEnergy()<<", ";
     out<<getRefSolSys()->Energy(planet, MJD())-CState::Energy();
   } else {
     out<<CState::x()-getRefSolSys()->x(planet, MJD())<<", ";
     out<<CState::y()-getRefSolSys()->y(planet, MJD())<<", ";
     out<<CState::z()-getRefSolSys()->z(planet, MJD())<<", ";
     out<<CState::dx()-getRefSolSys()->dx(planet, MJD())<<", ";
     out<<CState::dy()-getRefSolSys()->dy(planet, MJD())<<", ";
     out<<CState::dz()-getRefSolSys()->dz(planet, MJD())<<", ";
     out<<CState::r()-getRefSolSys()->r(planet, MJD())<<", ";
     out<<DifferenceAngle(CState::phi(), getRefSolSys()->phi(planet, MJD()))<<", ";
     out<<CState::theta()-getRefSolSys()->theta(planet, MJD())<<", ";
     out<<CState::dr()-getRefSolSys()->dr(planet, MJD())<<", ";
     out<<CState::dphi()-getRefSolSys()->dphi(planet, MJD())<<", ";
     out<<CState::dtheta()-getRefSolSys()->dtheta(planet, MJD())<<", ";
     out<<CState::a()-getRefSolSys()->a(planet, MJD())<<", ";
     out<<CState::e()-getRefSolSys()->e(planet, MJD())<<", ";
     out<<CState::i()-getRefSolSys()->i(planet, MJD())<<", ";
     out<<DifferenceAngle(CState::Omega(), getRefSolSys()->Omega(planet, MJD()))<<", ";
     out<<DifferenceAngle(CState::omega(), getRefSolSys()->omega(planet, MJD()))<<", ";
     out<<DifferenceAngle(CState::curlypi(), getRefSolSys()->curlypi(planet, MJD()))<<", ";
     out<<DifferenceAngle(CState::M(), getRefSolSys()->M(planet, MJD()))<<", ";
     out<<DifferenceAngle(CState::E(), getRefSolSys()->E(planet, MJD()))<<", ";
     out<<DifferenceAngle(CState::f(), getRefSolSys()->f(planet, MJD()))<<", ";
     out<<DifferenceAngle(CState::L(), getRefSolSys()->L(planet, this->MJD()))<<", ";
     out<<getRefSolSys()->Velocity(planet, BODY_STATE, 0.0, 0.0, 0.0, 0.0, 0.0,
      (*this))<<", ";
     out<<CState::KineticEnergy()-(getRefSolSys()->KineticEnergy(planet, MJD()))<<", ";
     out<<CState::PotentialEnergy()-(getRefSolSys()->PotentialEnergy(planet, MJD()))<<", ";
     out<<CState::Energy()-(getRefSolSys()->Energy(planet, MJD()));
     //AO
   }
} // void PrintDelta(...)


//------------------------------------------------------------------------------
void CSC::PrintNCInput(ofstream&out) const
{
  for (int i=0; i<getCtrlInputDim(); ++i)
    out<<", " <<m_input_ncp(i);
} // void PrintNCInput(ofstream&out) const


//------------------------------------------------------------------------------
void CSC::PrintNCOutput(ofstream&out) const
{
  for (int i=0; i<getCtrlOutputDim(); ++i)
    out<<", " <<m_nc_output(i);
} // void CSC::PrintNCOutput(ofstream&out) const


//------------------------------------------------------------------------------
void CSC::PrintProximity(ofstream&out, const int planet) const
{
  out << getRefSolSys()->Distance_km(
    planet, BODY_STATE, 0.0, 0.0, 0.0, 0.0, 0.0, (*this))<<", ";
  out << getRefSolSys()->Velocity_kms(
    planet, BODY_STATE, 0.0, 0.0, 0.0, 0.0, 0.0, (*this));
} // void PrintProximity(ofstream&out, const int planet) const
