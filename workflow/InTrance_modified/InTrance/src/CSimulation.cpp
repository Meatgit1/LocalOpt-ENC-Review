/*!
\file CSimulation.cpp
\version 1.7
\brief CSimulation class implementation.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de)> (AO)
\date 04.11.2006 (AO)

\b Purpose:
   This file provides the interface for the class CSimulation.
   CSimulation is the class that implements the simulation functionality.

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

// C standard library
#include <conio.h>

// STL
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <cassert>
#include <functional>
#include <numeric>
#include <time.h>
#include <time.h>
#include <algorithm>
#include <valarray>

// SpaceLib
#include <MLQuaternion.h>
#include <SpaceLib.h>
#include <MathLib.h>

// ENCLib
#include <CIndividual.h>
#include <CNC.h>

// SCLib
#include <CSCAttitudeHistory.h>
#include <CSC.h>
#include <NativeBridgeProbe.h>

// InTrance
#include <InTranceConst.h>
//#include <CInputSimulation.h>
#include <CInTranceFitness.h>
#include <CApproxFitness.h>
#include <CSimulation.h>
#include <VRMLSkel1.h>
#include <CConstraint.h>
#if defined using_constraints
#include <CConstraint.h>
#include <CMinCentralDistConstraint.h>
#include <CMaxSailTempConstraint.h>
#include <CProximityConstraint.h>
#endif // using_constraints


// COMMONCPP library
#if defined COMMONCPP
#include <TQueue.h>
#undef max
#undef min
#endif //COMMONCPP

using namespace std;
using namespace MathLib;
using namespace SpaceLib;

// see calcConeAngleApproximation() for explanation
const double HALFPI_MODIFICATION = 1E-8;

const double SAILJETTISONDISTANCE(6.0);
///< solar sail jettison distance for missions into the outer solar system [AU]

const double SAILJETTISONINTEGRATIONSTEPWIDTH(100.0);
///< solar sail jettison integration step width [days]

extern TrajectoryType tt; // global in CInTranceFitness.cpp
extern FitnessFunctionType jt; // global in CInTranceFitness.cpp

// initialize sensitivity factors (only once)
static bool init_sensitivity_factors=true;

SpaceLib::CState temp_state;

std::auto_ptr<constraints::CConstraint> MIN_CENTRAL_BODY_CONSTRAINT;
std::auto_ptr<constraints::CConstraint> MAX_SAIL_TEMP_CONSTRAINT;
std::auto_ptr<constraints::CConstraint> TRGT_PROX_CONSTRAINT;

// special constraint (to be removed, because this is not a proper design)
namespace constraints
{

//#define CONSEP 1
#undef CONSEP
#ifdef CONSEP
  /// special constraint object for CONSEP study
  class CCONSEPConstraint : virtual public CConstraint
  {
  private:
    bool m_switched;
    double m_sc_init_mass;

  public:
    CCONSEPConstraint(void) throw();
    CCONSEPConstraint(const CCONSEPConstraint& cs) throw();
    CCONSEPConstraint(CSimulation* sim) throw();
    ~CCONSEPConstraint(void) throw();
    virtual bool isMatched(void) throw();
  };

  //----------------------------------------------------------------------------
  inline CCONSEPConstraint::CCONSEPConstraint(void) throw() 
    : CConstraint(), m_switched(false), m_sc_init_mass(-1.0) {};

  //----------------------------------------------------------------------------
  inline CCONSEPConstraint::CCONSEPConstraint(const CCONSEPConstraint& cs)
    throw() : CConstraint(cs), m_switched(false), m_sc_init_mass(-1.0) {}

  //----------------------------------------------------------------------------
  inline CCONSEPConstraint::CCONSEPConstraint(CSimulation* sim) throw()
    : CConstraint(sim), m_switched(false), m_sc_init_mass(-1.0) {}

  //----------------------------------------------------------------------------
  inline CCONSEPConstraint::~CCONSEPConstraint(void) throw() {};

  //----------------------------------------------------------------------------
  inline bool CCONSEPConstraint::isMatched(void) throw()
  {
    if (m_sim->getFlightTime() == 0.0)
    { 
      CRIT* rit = dynamic_cast<CRIT*>(m_sim->getCurrentSC().getPropSys());
      if (rit !=0)
      {
        if (m_sc_init_mass<0.0)
        {
          m_sc_init_mass = m_sim->getCurrentSC().getSCMass();
        }
        else
        {
          m_sim->getCurrentSC().setSCMass(m_sc_init_mass);
        }
        rit->setCoefficient(MAX_THRUSTER_POWER, 11.11);
        rit->setCoefficient(OPT_THRUSTER_POWER, 11.11);
        rit->setCoefficient(MIN_THRUSTER_POWER, 7.2215);
        rit->setCoefficient(THRUST_0, 0.0);
        rit->setCoefficient(THRUST_1, 21.7E-3);
        rit->setCoefficient(M_FLOW_0, 0.335E-6);
        rit->setCoefficient(M_FLOW_1, 0.0);
        rit->update();                        // to update internal parameters
      }
      m_switched = false;
    } // if (m_sim->getFlightTime() == 0.0)
    else
      if (m_switched == false && m_sim->getCurrentSC().CState::r() >= 3.026)
      {
        CRIT* rit = dynamic_cast<CRIT*>(m_sim->getCurrentSC().getPropSys());
        if (rit !=0)
        {
          rit->setCoefficient(MAX_THRUSTER_POWER, 4.02);
          rit->setCoefficient(OPT_THRUSTER_POWER, 4.02);
          rit->setCoefficient(MIN_THRUSTER_POWER, 2.613);
          rit->setCoefficient(THRUST_0, 0.0);
          rit->setCoefficient(THRUST_1, 33.73E-3);
          rit->setCoefficient(M_FLOW_0, 0.926E-6);
          rit->setCoefficient(M_FLOW_1, 0.0);
          rit->update();                      // to update internal parameters
        }
        // reduce spacecraft mass by 545.0 kg (jettison of 1st stage)
        const double mass_stage_1 = 545.0; // [kg]
        const double mass = m_sim->getCurrentSC().getSCMass();
        m_sim->getCurrentSC().setSCMass(Global::max(0.0, (mass-mass_stage_1)));
        m_switched = true;
      }
      return true;
  } // bool CCONSEPConstraint::isMatched(void) const throw()
#endif // CONSEP

  //----------------------------------------------------------------------------
  /// special constraint object for GTOC2 problem
  class CGTOC2Constraint : virtual public CConstraint
  {
  private:
 
  public:
    CGTOC2Constraint(void) throw();
    CGTOC2Constraint(const CGTOC2Constraint& cs) throw();
    CGTOC2Constraint(CSimulation* sim) throw();
    ~CGTOC2Constraint(void) throw();
    virtual bool isMatched(void) throw();
  };

  //----------------------------------------------------------------------------
  inline CGTOC2Constraint::CGTOC2Constraint(void) throw() 
    : CConstraint() {}

  //----------------------------------------------------------------------------
  inline CGTOC2Constraint::CGTOC2Constraint(const CGTOC2Constraint& cs)
    throw() : CConstraint(cs) {}

  //----------------------------------------------------------------------------
  inline CGTOC2Constraint::CGTOC2Constraint(CSimulation* sim) throw()
    : CConstraint(sim) {}

  //----------------------------------------------------------------------------
  inline CGTOC2Constraint::~CGTOC2Constraint(void) throw() {};

  //----------------------------------------------------------------------------
  inline bool CGTOC2Constraint::isMatched(void) throw()
  {
      return true;
  } // bool CGTOC2Constraint::isMatched(void) const throw()


} // namespace constraints

//------------------------------------------------------------------------------
bool IntegrationStopFunction(double t, Vector& y, void* pA)
{
  CSimulation* sim = static_cast<CSimulation*>(pA);
  if (sim != 0)
  {
    sim->updateProximity();
    if (sim->checkForNewBestStep() == true) return true;
  } // if (sim != 0)
  return false;
} // bool GTOC2StopFunction(double t, Vector y, void* pA)

std::auto_ptr<constraints::CConstraint> GTOC2_CONSTRAINT;

#if defined CONSEP
std::auto_ptr<constraints::CConstraint> CONSEP_CONSTRAINT;
#endif // CONSEP

//-----------------------------------------------------------------------------
inline double calcPropMassFitness(const double m_prop_init
                                , const double m_prop_current)
{
  return 2.0e3*(m_prop_init/(2.0*m_prop_init-m_prop_current)-0.5);
}

//-----------------------------------------------------------------------------
template<class T, class U>
inline void CopyDataStructure(const T& lhs, U& s2)
{
  if (lhs.size() != 0)
  {
    s2.resize(lhs.size());
    for (unsigned int i=0; i<lhs.size();++i) s2[i]=lhs[i];
  } // if (lhs.size() != 0)
} // void CopyDataStructure(const T& lhs, U& s2)


//-----------------------------------------------------------------------------
double calcPerihelionTime(BodyID body, double time,
                   CSolarSystem* solsys,double eps)
{
  double d_t(solsys->f(body,time)/TWOPI*solsys->OrbitPeriod(body,time)/2.0);
  double f_n(0.0);
  if ((f_n=solsys->f(body,time))>PI)
  {
    f_n=solsys->f(body,time-=solsys->OrbitPeriod(body,time)/2.0);
    d_t=f_n/TWOPI*solsys->OrbitPeriod(body,time)/2.0;
  }
  do
  {
    f_n = solsys->f(body,time-=d_t);
    if (f_n>PI)
    {
      time+=d_t;
      d_t*=0.5;
    }
  } while (d_t>eps || f_n>PI);
  return time;
}


//-----------------------------------------------------------------------------
double calcAngleTime(BodyID body, double angle, double t_start
                     ,CSolarSystem* solsys, double eps)
{
  double d_t(solsys->OrbitPeriod(body,t_start)/12.0);
  double t_l(t_start),t_u(t_start+d_t);
  double angle1(AssertAngularRange2pi(angle));
  double angle2(angle1+TWOPI);
  double a_l(0.0),a_u(0.0);
  bool between(false);
  do
  {
    do
    {
      a_l=solsys->phi(body,t_l);a_u=solsys->phi(body,t_u);
      while (a_l<0.0 || a_u<0.0) {a_l+=TWOPI;a_u+=TWOPI;}
      if ((a_u>=angle1 && angle1>=a_l) || (a_u>=angle2 && angle1>=a_l)) between=true;
      else {t_l+=d_t;t_u+=d_t;}
    } while (between==false);
    if (d_t>eps) { between=false;t_l-=d_t;t_u=(t_l+(d_t*=0.1));}
  }
  while (between==false);
  return 0.5*(t_l+t_u);
}


//------------------------------------------------------------------------------
void CSimulation::fillStatusRecRendezvous(SimStatusRec* s) 
{
  const int min_rows = 30, fw[] = {20,30,20,20};
  const int dim(sizeof(fw)/sizeof(fw[0]));
  const std::string col_names[dim] = {"variable","launch","s/c","best"};
  while(s->data.size()<min_rows) s->data.push_back("");
  std::ostringstream os,os2;os << left;
  int row = -1;
  // column header
  for (int i=0;i<dim;++i) os << " " << setw(fw[i]) << col_names[i];
  s->data[++row]= os.str();

  //// body name
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "name";
  os << " " << setw(fw[1]) << m_solsys->Name(getInitialBody());
  os << " " << setw(fw[2]) << "spacecraft";
  os << " " << setw(fw[3]) << m_solsys->Name(getTargetBodyId());
  s->data[++row]=os.str();

  double t_best = getStartTime()+m_time_best;

  // central body name
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "central body";
  os << " " << setw(fw[1]) << m_solsys->Name(m_solsys->State(
    getInitialBody(),getStartTime()).getCenterBodyId());
  os << " " << setw(fw[2]) << m_solsys->Name(
    m_spacecraft.getCenterBodyId());
  os << " " << setw(fw[3]) << m_solsys->Name(m_solsys->State(
    getTargetBodyId(),t_best).getCenterBodyId());
  s->data[++row]=os.str();

  // time [MJD]
  os.clear();os.str("");
  os2.clear();os2.str("");os2 << fixed << setprecision(1);
  os2 << " (" << getMinStartTime() << ") " << getStartTime()
      << " (" << getMaxStartTime() << ")";
  os << " " << setw(fw[0]) << "date [MJD]";
  os << ":" << setw(fw[1]) << os2.str();
  os << " " << setw(fw[2]) << fixed << setprecision(2) << t_best;
  os << " " << setw(fw[3]) << fixed << setprecision(2) << t_best;
  s->data[++row]=os.str();

  // spacecraft mass [kg]
  os.clear();os.str("");
  os2.clear();os2.str("");os2 << fixed << setprecision(1);
  os2 << " (" << m_spacecraft.getPropSys()->getMinInitialPropellantMass()<< ")";
  os2 << " " << m_spacecraft.getInitialPropellantMass();
  os2 << " (" << m_spacecraft.getPropSys()->getMaxInitialPropellantMass()<< ")";
  os << " " << setw(fw[0]) << "prop. mass [kg]";
  os << ":" << setw(fw[1]) << os2.str();
  os << " " << setw(fw[2]);
  os << fixed << setprecision(2) << m_best_sc.getPropellantMass();
  s->data[++row]=os.str();

  // flight time [days]
  os.clear();os.str("");os2.clear();os2.str("");
  os << " " << setw(fw[0]) << "flight time [days]";
  os << ":" << setw(fw[1]) << " ";
  os2 << fixed<<setprecision(2)<<getBestFlightTime()<<"d";
  os2 << ", " << getBestFlightTime()*DAY2JYR<<"y";
  os2 << " (" << fixed<<setprecision(2)<<getIntegrationInterval()<< "d";
  os2 << ", " << fixed<<setprecision(2)<<getIntegrationInterval()*DAY2JYR;
  os2 <<")";
  os  << " " << setw(fw[3]) << os2.str();
  s->data[++row]=os.str();

  // semimajor axis [AU]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "a [AU]";
  os << ":" << setw(fw[1]) << showpos << scientific
    << setprecision(4)
    << m_solsys->a(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::a();
  os << " " << setw(fw[3]) << showpos << scientific
    << setprecision(4)
    << m_solsys->a(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // eccentricity
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "e";
  os << ":" << setw(fw[1]) << noshowpos << fixed
    << setprecision(4)
    << m_solsys->e(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::e();;
  os << " " << setw(fw[3]) << m_solsys->e(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // inclination [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "i [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific
    << setprecision(3) << m_solsys->i(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::i();
  os << " " << setw(fw[3]) << m_solsys->i(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // Omega [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "Omega [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->Omega(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::Omega();
  os << " " << setw(fw[3]) << m_solsys->Omega(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // omega [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "omega [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->omega(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::omega();
  os << " " << setw(fw[3]) << m_solsys->omega(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // anomaly M [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "anomaly M [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->M(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::M();
  os << " " << setw(fw[3]) << m_solsys->M(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // distance from cental body r [AU]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "distance r [AU]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->r(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::r();
  os << " " << setw(fw[3]) << m_solsys->r(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // angle phi [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "phi [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->phi(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::phi();
  os << " " << setw(fw[3]) << m_solsys->phi(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // angle theta [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "theta [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->theta(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::theta();
  os << " " << setw(fw[3]) << m_solsys->theta(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // dr/dt [AU/d]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "dr/dt [AU/d]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->dr(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::dr()*m_max_dr_AU;
  os << " " << setw(fw[3]) << m_solsys->dr(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // dphi/dt [rad/d]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "dphi/dt [rad/d]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->dphi(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::dphi();
  os << " " << setw(fw[3]) << m_solsys->dphi(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // dtheta/dt [rad/d]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "dtheta/dt [rad/d]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->dtheta(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::dtheta();
  os << " " << setw(fw[3]) << m_solsys->dtheta(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

//> may be deleted
  //if (getTargetStateType()==ORBITAL_ELEMENTS_STATE)
  //{
  //  const std::string el_names[] =
  //  {
  //    "rem. time a [d,y]"
  //   ,"rem. time e [d,y]"
  //   ,"rem. time i [d,y]"
  //   ,"rem. time Om[d,y]"
  //   ,"rem. time om[d,y]"
  //   ,"rem. time M [d,y]"
  //  };

  //  for (int i=0;i<N_KEPLERIAN_ELEMENTS;++i)
  //  {
  //    const KEPLERIAN_ELEMENT e = KEPLERIAN_ELEMENT(i);
  //    if (useTargetOrbitalElement(e) == true)
  //    {
  //      // remaining time to change semimajor axis
  //      os.clear();os.str("");
  //      os << " " << setw(fw[0]) << el_names[e];
  //      os << ":" << setw(fw[1]) << " ";
  //      os << " " << setw(fw[2]/2) << noshowpos << fixed << right << setprecision(3) 
  //        << m_best_element_change_rest_times[e];
  //      os << "," << setw(fw[2]/2) << noshowpos << fixed << left << setprecision(3)
  //        << m_best_element_change_rest_times[e]*DAY2JYR;
  //      s->data[++row]=os.str();
  //    }
  //  }

  //  const double change_time1 = 
  //    sqrt(pow(m_best_element_change_rest_times,2.0).sum());
  //  os.clear();os.str("");
  //  os << " " << setw(fw[0]) << "rem. time  [d,y]";
  //  os << ":" << setw(fw[1]) << " ";
  //  os << " " << setw(fw[2]/2) << noshowpos << fixed << right << setprecision(3)
  //    << change_time1;
  //  os << " " << setw(fw[2]/2) << noshowpos << fixed << left << setprecision(3)
  //    << change_time1*DAY2JYR;
  //  s->data[++row] = os.str();

  //  const double change_time2 = m_best_element_change_rest_times.sum();
  //  os.clear();os.str("");
  //  os << " " << setw(fw[0]) << "r. time sum[d,y]";
  //  os << ":" << setw(fw[1]) << " ";
  //  os << " " << setw(fw[2]/2) << noshowpos << fixed << right << setprecision(3)
  //    << change_time2;
  //  os << " " << setw(fw[2]/2) << noshowpos << fixed << left << setprecision(3)
  //    << change_time2*DAY2JYR;
  //  s->data[++row] = os.str();
  //}
//< may be deleted

  // initial excess velocity
  os.clear();os.str("");
  os << "" << setw(31) << "C3          [AU/d](min,cur,max)" << ":";
  os << showpos << scientific << setprecision(2) <<"("<< getMinLaunchVelocity()
    <<","<< m_initial_velocity  << "," << getMaxLaunchVelocity() << ")";
  s->data[++row]=os.str();

  // initial excess velocity vector azimuth angle
  os.clear();os.str("");
  os << "" << setw(31) << "C3 azimuth   [rad](min,cur,max)" << ":";
  os << showpos << scientific << setprecision(2) <<"("<< getMinLaunchAzimuthAngle()
    <<","<< m_initial_azimuth  << "," << getMaxLaunchAzimuthAngle() << ")";
  s->data[++row]=os.str();

  // initial excess velocity vector elevation angle
  os.clear();os.str("");
  os << "" << setw(31) << "C3 elevation [rad](min,cur,max)" << ":";
  os << showpos << scientific << setprecision(2) <<"("<< getMinLaunchElevationAngle()
    <<","<< m_initial_elevation  << "," << getMaxLaunchElevationAngle() << ")";
  s->data[++row]=os.str();

}

//-----------------------------------------------------------------------------
void CSimulation::fillStatusRecFlyby(SimStatusRec* s)
{
  fillStatusRecRendezvous(s);
} // void fillStatusRecFlyby(SimStatusRec* s)

//-----------------------------------------------------------------------------
void CSimulation::fillStatusRecRetrogradeImpact(SimStatusRec* s)
{
  fillStatusRecRendezvous(s);
} // void fillStatusRecRetrogradeImpact(SimStatusRec* s)

//-----------------------------------------------------------------------------
void CSimulation::fillStatusRecActProblem(SimStatusRec* s)
{
  const int min_rows = 30, fw[] = {20,30,20,20};
  const int dim(sizeof(fw)/sizeof(fw[0]));
  const std::string col_names[dim] = {"variable","launch","s/c","best"};
  while(s->data.size()<min_rows) s->data.push_back("");
  std::ostringstream os,os2;os << left;
  int row = -1;
  // column header
  for (int i=0;i<dim;++i) os << " " << setw(fw[i]) << col_names[i];
  s->data[++row]= os.str();

  // body name
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "name";
  os << " " << setw(fw[1]) << m_solsys->Name(getInitialBody());
  os << " " << setw(fw[2]) << "spacecraft";
  os << " " << setw(fw[3]) << m_solsys->Name(getTargetBodyId());
  s->data[++row]=os.str();

  double t_best = getStartTime()+m_time_best;

  // central body name
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "central body";
  os << " " << setw(fw[1]) << m_solsys->Name(m_solsys->State(
    getInitialBody(),getStartTime()).getCenterBodyId());
  os << " " << setw(fw[2]) << m_solsys->Name(
    m_spacecraft.getCenterBodyId());
  os << " " << setw(fw[3]) << m_solsys->Name(m_solsys->State(
    getTargetBodyId(),t_best).getCenterBodyId());
  s->data[++row]=os.str();

  // time [MJD]
  os.clear();os.str("");
  os2.clear();os2.str("");os2 << fixed << setprecision(1);
  os2 << " (" << getMinStartTime() << ") " << getStartTime()
      << " (" << getMaxStartTime() << ")";
  os << " " << setw(fw[0]) << "date [MJD]";
  os << ":" << setw(fw[1]) << os2.str();
  os << " " << setw(fw[2]) << fixed << setprecision(2) << t_best;
  os << " " << setw(fw[3]) << fixed << setprecision(2) << t_best;
  s->data[++row]=os.str();

  // flight time [days]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "flight time [days]";
  os << ":" << setw(fw[1]) << " ";
  os2.clear();os2.str("");
  os2 << fixed<<setprecision(2)<<getBestFlightTime()<<"d";
  os2 << ", " << getBestFlightTime()*DAY2JYR<<"y";
  os2 << " (" << fixed<<setprecision(2)<<getIntegrationInterval()<< "d";
  os2 << ", " << fixed<<setprecision(2)<<getIntegrationInterval()*DAY2JYR;
  os2 <<")";
  os << " " << setw(fw[3]) << os2.str();
  s->data[++row]=os.str();

  // semimajor axis [AU]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "a [AU]";
  os << ":" << setw(fw[1]) << showpos << scientific
    << setprecision(4)
    << m_solsys->a(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::a();
  os << " " << setw(fw[3]) << showpos << scientific
    << setprecision(4)
    << m_solsys->a(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // eccentricity
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "e";
  os << ":" << setw(fw[1]) << noshowpos << fixed
    << setprecision(4)
    << m_solsys->e(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::e();;
  os << " " << setw(fw[3]) << m_solsys->e(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // inclination [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "i [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific
    << setprecision(3) << m_solsys->i(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::i();
  os << " " << setw(fw[3]) << m_solsys->i(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // Omega [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "Omega [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->Omega(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::Omega();
  os << " " << setw(fw[3]) << m_solsys->Omega(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // omega [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "omega [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->omega(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::omega();
  os << " " << setw(fw[3]) << m_solsys->omega(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // anomaly M [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "anomaly M [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->M(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::M();
  os << " " << setw(fw[3]) << m_solsys->M(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // distance from cental body r [AU]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "distance r [AU]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->r(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::r();
  os << " " << setw(fw[3]) << m_solsys->r(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // angle phi [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "phi [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->phi(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::phi();
  os << " " << setw(fw[3]) << m_solsys->phi(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // angle theta [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "theta [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->theta(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::theta();
  os << " " << setw(fw[3]) << m_solsys->theta(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // dr/dt [AU/d]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "dr/dt [AU/d]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->dr(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::dr()*m_max_dr_AU;
  os << " " << setw(fw[3]) << m_solsys->dr(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // dphi/dt [rad/d]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "dphi/dt [rad/d]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->dphi(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::dphi();
  os << " " << setw(fw[3]) << m_solsys->dphi(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // dtheta/dt [rad/d]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "dtheta/dt [rad/d]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->dtheta(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::dtheta();
  os << " " << setw(fw[3]) << m_solsys->dtheta(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // X [AU]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "X [AU]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->x(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::x();
  os << " " << setw(fw[3]) << m_solsys->x(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // Y [AU]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "Y [AU]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->y(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::y();
  os << " " << setw(fw[3]) << m_solsys->y(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // Z [AU]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "Z [AU]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->z(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::z();
  os << " " << setw(fw[3]) << m_solsys->z(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // dX/dt [AU/d]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "dX/dt [AU/d]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->dx(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::dx();
  os << " " << setw(fw[3]) << m_solsys->dx(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // dY /dt [AU/dt]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "dY/dt [AU/d]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->dy(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::dy();
  os << " " << setw(fw[3]) << m_solsys->dy(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // dZ/dt [AU/d]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "dZ/dt [AU/d]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->dz(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::dz();
  os << " " << setw(fw[3]) << m_solsys->dz(getTargetBodyId(), t_best);
  s->data[++row]=os.str();

  // blank line
  os.clear();os.str("");s->data[++row]=os.str();

  // max. # of steps (used and absolute)
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "# of steps (max)";
  os << ":" << setw(fw[1]) << noshowpos << dec << m_steps_max_used;
  os << " " << setw(fw[2]) << "# of steps (max,abs)";
  os << ":" << setw(fw[3]) << noshowpos << dec << getMaxStepsCount();
  s->data[++row]=os.str();

  // current accuracy limits (relative distance and velocity)
  os.clear();os.str("");os2.clear();os2.str("");
  const double max_dr = (m_max_dr_AU*m_accuracy_factor);
  os << " " << setw(fw[0]) << "delta_r [AU(km)]";
  os2<<noshowpos<<scientific<<setprecision(3)
   << max_dr <<"("<<max_dr*AU2KM<<")";
  os << ":" << setw(fw[1]) << os2.str();
  os << " " << setw(fw[2]) << "delta_v [AU/d(km/s)]";
  os2.clear();os2.str("");
  os2<<noshowpos<<scientific<<setprecision(1)<<(m_max_dv_AUd*m_accuracy_factor)
    <<"("<<(m_max_dv_kms*m_accuracy_factor)<<")";
  os << ":" << setw(fw[3]) << os2.str();;
  s->data[++row]=os.str();

  // initial excess velocity
  os.clear();os.str("");
  os << "" << setw(31) << "C3          [AU/d](min,cur,max)" << ":";
  os << showpos << scientific << setprecision(2) <<"("<< getMinLaunchVelocity()
    <<","<< m_initial_velocity  << "," << getMaxLaunchVelocity() << ")";
  s->data[++row]=os.str();

  // initial excess velocity vector azimuth angle
  os.clear();os.str("");
  os << "" << setw(31) << "C3 azimuth   [rad](min,cur,max)" << ":";
  os << showpos << scientific << setprecision(2) <<"("<< getMinLaunchAzimuthAngle()
    <<","<< m_initial_azimuth  << "," << getMaxLaunchAzimuthAngle() << ")";
  s->data[++row]=os.str();

  // initial excess velocity vector elevation angle
  os.clear();os.str("");
  os << "" << setw(31) << "C3 elevation [rad](min,cur,max)" << ":";
  os << showpos << scientific << setprecision(2) <<"("<< getMinLaunchElevationAngle()
    <<","<< m_initial_elevation  << "," << getMaxLaunchElevationAngle() << ")";
  s->data[++row]=os.str();

} // void CSimulation::fillStatusRecActProblem(SimStatusRec* s)

//-----------------------------------------------------------------------------
void CSimulation::fillStatusRecGTOC2Problem(SimStatusRec* s)
{
  BodyID id = m_gtoc2_best_bodies[0];

  const int min_rows = 30, fw[] = {20,30,20,20};
  const int dim(sizeof(fw)/sizeof(fw[0]));
  const std::string col_names[dim] = {"variable","launch","s/c","best"};
  while(s->data.size()<min_rows) s->data.push_back("");
  std::ostringstream os,os2;os << left;
  int row = -1;
  // column header
  for (int i=0;i<dim;++i) os << " " << setw(fw[i]) << col_names[i];
  s->data[++row]= os.str();

  //// body name
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "name";
  os << " " << setw(fw[1]) << m_solsys->Name(getInitialBody());
  os << " " << setw(fw[2]) << "spacecraft";
  os << " " << setw(fw[3]) << m_solsys->Name(id);
  s->data[++row]=os.str();

  double t_best = getStartTime()+m_time_best;

  // central body name
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "central body";
  os << " " << setw(fw[1]) << m_solsys->Name(m_solsys->State(
    getInitialBody(),getStartTime()).getCenterBodyId());
  os << " " << setw(fw[2]) << m_solsys->Name(
    m_spacecraft.getCenterBodyId());
  os << " " << setw(fw[3]) << m_solsys->Name(m_solsys->State(
    id,t_best).getCenterBodyId());
  s->data[++row]=os.str();

  // time [MJD]
  os.clear();os.str("");
  os2.clear();os2.str("");os2 << fixed << setprecision(1);
  os2 << " (" << getMinStartTime() << ") " << getStartTime()
      << " (" << getMaxStartTime() << ")";
  os << " " << setw(fw[0]) << "date [MJD]";
  os << ":" << setw(fw[1]) << os2.str();
  os << " " << setw(fw[2]) << fixed << setprecision(2) << t_best;
  os << " " << setw(fw[3]) << fixed << setprecision(2) << t_best;
  s->data[++row]=os.str();

  // spacecraft mass [kg]
  os.clear();os.str("");
  os2.clear();os2.str("");os2 << fixed << setprecision(1);
  os2 << " (" << m_spacecraft.getPropSys()->getMinInitialPropellantMass()<< ")";
  os2 << " " << m_spacecraft.getInitialPropellantMass();
  os2 << " (" << m_spacecraft.getPropSys()->getMaxInitialPropellantMass()<< ")";
  os << " " << setw(fw[0]) << "prop. mass [kg]";
  os << ":" << setw(fw[1]) << os2.str();
  os << " " << setw(fw[2]);
  os << fixed << setprecision(2) << m_best_sc.getPropellantMass();
  s->data[++row]=os.str();

  // flight time [days]
  os.clear();os.str("");os2.clear();os2.str("");
  os << " " << setw(fw[0]) << "flight time [days]";
  os << ":" << setw(fw[1]) << " ";
  os2 << fixed<<setprecision(2)<<getBestFlightTime()<<"d";
  os2 << ", " << getBestFlightTime()*DAY2JYR<<"y";
  os2 << " (" << fixed<<setprecision(2)<<getIntegrationInterval()<< "d";
  os2 << ", " << fixed<<setprecision(2)<<getIntegrationInterval()*DAY2JYR;
  os2 <<")";
  os  << " " << setw(fw[3]) << os2.str();
  s->data[++row]=os.str();

  // semimajor axis [AU]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "a [AU]";
  os << ":" << setw(fw[1]) << showpos << scientific
    << setprecision(4)
    << m_solsys->a(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::a();
  os << " " << setw(fw[3]) << showpos << scientific
    << setprecision(4)
    << m_solsys->a(id, t_best);
  s->data[++row]=os.str();

  // eccentricity
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "e";
  os << ":" << setw(fw[1]) << noshowpos << fixed
    << setprecision(4)
    << m_solsys->e(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::e();;
  os << " " << setw(fw[3]) << m_solsys->e(id, t_best);
  s->data[++row]=os.str();

  // inclination [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "i [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific
    << setprecision(3) << m_solsys->i(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::i();
  os << " " << setw(fw[3]) << m_solsys->i(id, t_best);
  s->data[++row]=os.str();

  // Omega [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "Omega [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->Omega(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::Omega();
  os << " " << setw(fw[3]) << m_solsys->Omega(id, t_best);
  s->data[++row]=os.str();

  // omega [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "omega [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->omega(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::omega();
  os << " " << setw(fw[3]) << m_solsys->omega(id, t_best);
  s->data[++row]=os.str();

  // anomaly M [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "anomaly M [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->M(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::M();
  os << " " << setw(fw[3]) << m_solsys->M(id, t_best);
  s->data[++row]=os.str();

  // distance from cental body r [AU]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "distance r [AU]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->r(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::r();
  os << " " << setw(fw[3]) << m_solsys->r(id, t_best);
  s->data[++row]=os.str();

  // angle phi [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "phi [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->phi(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::phi();
  os << " " << setw(fw[3]) << m_solsys->phi(id, t_best);
  s->data[++row]=os.str();

  // angle theta [rad]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "theta [rad]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->theta(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::theta();
  os << " " << setw(fw[3]) << m_solsys->theta(id, t_best);
  s->data[++row]=os.str();

  // dr/dt [AU/d]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "dr/dt [AU/d]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->dr(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::dr()*m_max_dr_AU;
  os << " " << setw(fw[3]) << m_solsys->dr(id, t_best);
  s->data[++row]=os.str();

  // dphi/dt [rad/d]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "dphi/dt [rad/d]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->dphi(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::dphi();
  os << " " << setw(fw[3]) << m_solsys->dphi(id, t_best);
  s->data[++row]=os.str();

  // dtheta/dt [rad/d]
  os.clear();os.str("");
  os << " " << setw(fw[0]) << "dtheta/dt [rad/d]";
  os << ":" << setw(fw[1]) << showpos << scientific << setprecision(3)
    << m_solsys->dtheta(getInitialBody(), getStartTime());
  os << " " << setw(fw[2]) << m_best_sc.CState::dtheta();
  os << " " << setw(fw[3]) << m_solsys->dtheta(id, t_best);
  s->data[++row]=os.str();

  // initial excess velocity
  os.clear();os.str("");
  os << "" << setw(31) << "C3          [AU/d](min,cur,max)" << ":";
  os << showpos << scientific << setprecision(2) <<"("<< getMinLaunchVelocity()
    <<","<< m_initial_velocity  << "," << getMaxLaunchVelocity() << ")";
  s->data[++row]=os.str();

  // initial excess velocity vector azimuth angle
  os.clear();os.str("");
  os << "" << setw(31) << "C3 azimuth   [rad](min,cur,max)" << ":";
  os << showpos << scientific << setprecision(2) <<"("<< getMinLaunchAzimuthAngle()
    <<","<< m_initial_azimuth  << "," << getMaxLaunchAzimuthAngle() << ")";
  s->data[++row]=os.str();

  // initial excess velocity vector elevation angle
  os.clear();os.str("");
  os << "" << setw(31) << "C3 elevation [rad](min,cur,max)" << ":";
  os << showpos << scientific << setprecision(2) <<"("<< getMinLaunchElevationAngle()
    <<","<< m_initial_elevation  << "," << getMaxLaunchElevationAngle() << ")";
  s->data[++row]=os.str();

} // void fillStatusRecGTOC2Problem(SimStatusRec* s)

//-----------------------------------------------------------------------------
void CSimulation::fillStatusRecTest1(SimStatusRec* s) 
{
  fillStatusRecRendezvous(s);
} // void fillStatusRecTest1(SimStatusRec* s)

//-----------------------------------------------------------------------------
void CSimulation::fillStatusRecTest2(SimStatusRec* s)
{
  fillStatusRecRendezvous(s);
} // void fillStatusRecTest1(SimStatusRec* s)

//-----------------------------------------------------------------------------
void CSimulation::fillStatusRecTest3(SimStatusRec* s)
{
  fillStatusRecRendezvous(s);
} // void fillStatusRecTest1(SimStatusRec* s)

//-----------------------------------------------------------------------------
void CSimulation::fillStatusRecTest4(SimStatusRec* s)
{
  fillStatusRecRendezvous(s);
} // void fillStatusRecTest1(SimStatusRec* s)

//-----------------------------------------------------------------------------
void CSimulation::fillStatusRecTest5(SimStatusRec* s)
{
  fillStatusRecRendezvous(s);
} // void fillStatusRecTest1(SimStatusRec* s)

//-----------------------------------------------------------------------------
void CSimulation::fillStatusRecTest6(SimStatusRec* s)
{
  fillStatusRecRendezvous(s);
} // void fillStatusRecTest1(SimStatusRec* s)

//-----------------------------------------------------------------------------
void CSimulation::fillStatusRecTest7(SimStatusRec* s)
{
  fillStatusRecRendezvous(s);
} // void fillStatusRecTest1(SimStatusRec* s)

//-----------------------------------------------------------------------------
void CSimulation::fillStatusRecTest8(SimStatusRec* s)
{
  fillStatusRecRendezvous(s);
} // void fillStatusRecTest1(SimStatusRec* s)

//-----------------------------------------------------------------------------
void CSimulation::fillStatusRecTest9(SimStatusRec* s)
{
  fillStatusRecRendezvous(s);
} // void fillStatusRecTest1(SimStatusRec* s)

//-----------------------------------------------------------------------------
void CSimulation::fillStatusRec(SimStatusRec* s)
{
  if (s==0) return;
  switch (tt)
  {
  case RENDEZVOUS:fillStatusRecRendezvous(s);break;
  case FLYBY:fillStatusRecFlyby(s);break;
  case RIMPACT:fillStatusRecRetrogradeImpact(s);break;
  case ACTPROBLEM:fillStatusRecActProblem(s);break;
  case SPIPROBLEM:fillStatusRecActProblem(s);break;
  case GTOC2TRAJ:fillStatusRecGTOC2Problem(s);break;
  case TRAJ_TEST1:fillStatusRecTest1(s);break;
  case TRAJ_TEST2:fillStatusRecTest2(s);break;
  case TRAJ_TEST3:fillStatusRecTest3(s);break;
  case TRAJ_TEST4:fillStatusRecTest4(s);break;
  case TRAJ_TEST5:fillStatusRecTest5(s);break;
  case TRAJ_TEST6:fillStatusRecTest6(s);break;
  case TRAJ_TEST7:fillStatusRecTest7(s);break;
  case TRAJ_TEST8:fillStatusRecTest8(s);break;
  case TRAJ_TEST9:fillStatusRecTest9(s);break;
  } // switch (tt)
} // void fillStatusRec(SimStatusRec* s)


//-----------------------------------------------------------------------------
void CSimulation::initEOM(
  CSC* sc
 ,RK4*&                  integrator_RK4
 ,RKF54*&                integrator_RKF54
 ,DOPRI8*&               integrator_DOPRI8
 ,DE*&                   integrator_DESG)
{
  // get translational equations of motion
  EOM_RK_func pEOM_RK = sc->getPropSys()->getTEOM();

  // initialize Runge-Kutta 4 integrator
  integrator_RK4=new RK4(pEOM_RK, 8, sc);
  Assure(integrator_RK4, "CSimulation::InitEOM(.....)=>(RK4):");

  // initialize Runge-Kutta-Fehlberg 4(5) integrator
  integrator_RKF54=new RKF54(pEOM_RK, 8, sc);
  Assure(integrator_RKF54, "CSimulation::InitEOM(.....)=>(RKF54):");

  // initialize Dormand-Prince 7(8) integrator
  integrator_DOPRI8=new DOPRI8(pEOM_RK, 8, sc);
  Assure(integrator_DOPRI8, "CSimulation::InitEOM(.....)=>(DOPRI8):");

  // initialize 4th order Runge-Kutta multistep integrator
  integrator_DESG=new DE(pEOM_RK, 8, sc);
  Assure(integrator_DESG, "CSimulation::InitEOM(.....)=>(DESG):");

} // end CSimulation::InitEOM


//-----------------------------------------------------------------------------
void CSimulation::initEOM(CSC* sc)
{
  // set type of translational equations of motion
  if (isHeliocentric()) sc->getPropSys()->setEomCoordinateType(POLAR);
  if (isBodyCentric()) sc->getPropSys()->setEomCoordinateType(CARTESIAN);

  const EOM_RK_func pEOM_RK = sc->getPropSys()->getTEOM();

  setIntegrator(new RK4(pEOM_RK, getIntegrationVectorDimension(), sc));
  setIntegrator(new RKF54(pEOM_RK, getIntegrationVectorDimension(), sc));
  setIntegrator(new DOPRI8(pEOM_RK, getIntegrationVectorDimension(), sc));
  setIntegrator(new DE(pEOM_RK, getIntegrationVectorDimension(), sc));

} // end CSimulation::InitEOM


//-----------------------------------------------------------------------------
void CSimulation::freeEOM(void)
{
  SafeDelete(m_SpacecraftOrbit_RK4);
  SafeDelete(m_SpacecraftOrbit_RKF54);
  SafeDelete(m_SpacecraftOrbit_DOPRI8);
  SafeDelete(m_SpacecraftOrbit_DESG);
} // end CSimulation::FreeEOM


//-----------------------------------------------------------------------------
// !Attention!: Don't try to replace calls to this function by something
// like "SCTraj traj = m_spacecraft_trajectory", because this will copy
// all elements of m_spacecraft_trajectory, but only the m_steps_best make
// sense.
void CSimulation::copyCurrentTrajectory(SCTraj& traj)
{
#if defined COMMONCPP
  SimpleGuard guard(m_sc_traj_lock); // lock trajectory data
#endif //COMMONCPP
  traj.resize(m_steps_best+1);        // resize argument to new trajectory size
  for (SCTraj::iterator i     = m_spacecraft_trajectory.begin(),
                        i_end = i+traj.size(),
                       i2     = traj.begin();
       i < i_end;++i,++i2)
      *i2 = *i;
 } // void CSimulation::copyCurrentTrajectory(SCTraj& traj)


//-----------------------------------------------------------------------------
void CSimulation::saveLastTrajectory(void)
{
#if defined COMMONCPP
  SimpleGuard guard1(m_last_saved_trajectory_lock);     // lock data structures
#endif //COMMONCPP
  // copy current s/c trajectory
  copyCurrentTrajectory(m_last_saved_trajectory);

  // determine launch date [MJD]
  if (m_last_saved_trajectory.size()>0)
    m_last_saved_start_time = m_last_saved_trajectory[0].MJD();
  else m_last_saved_start_time = 0;

  // save time steps in separate array
  // first lock because of multithreading
#if defined COMMONCPP
  SimpleGuard guard2(m_last_saved_time_steps_lock);
#endif //COMMONCPP

  // resize array
  m_last_saved_time_steps.resize(m_last_saved_trajectory.size());
  if (m_last_saved_time_steps.size()>0)
  {
    // and copy data
    for (uint i=0;i<m_last_saved_time_steps.size();++i)
    {
      m_last_saved_time_steps[i] =
        (m_last_saved_trajectory[i].MJD()-m_last_saved_start_time);
    }
  }
  // empty dependend data structures
  //m_last_saved_clock_angle_approx_values.resize(0);
  //m_last_saved_cone_angle_approx_values.resize(0);
  //m_last_saved_trajectory_approx.resize(0);
  m_last_saved_clock_angle_history.resize(0);
  m_last_saved_cone_angle_history.resize(0);
  m_last_saved_thrust_history.resize(0);
} // CSimulation::saveLastTrajectory(void)


//-----------------------------------------------------------------------------
void CSimulation::saveLastTrajectory( CIndividual& ind )
{
  m_last_saved_individual = ind;
  switch (m_sim_type)
  {
    case STANDARD:
      m_last_saved_start_velocity =
       ind.getMappedAllele(ind.getLength()-4
       , getMinLaunchVelocity(), getMaxLaunchVelocity());

      m_last_saved_start_azimuth =
       ind.getCircMappedAllele(ind.getLength()-3
       , getMinLaunchAzimuthAngle(), getMaxLaunchAzimuthAngle());

      m_last_saved_start_elevation = ind.getMappedAllele(ind.getLength()-2
       , getMinLaunchElevationAngle(), getMaxLaunchElevationAngle());
      break;
    case APPROXIMATION:
    default:
      m_last_saved_start_velocity  = 0.0;
      m_last_saved_start_azimuth   = 0.0;
      m_last_saved_start_elevation = 0.0;
  } // switch (m_sim_type)
  saveLastTrajectory();
} // CSimulation::saveLastTrajectory(const CIndividual&)


//-----------------------------------------------------------------------------
void CSimulation::copyLastSavedTrajectory(SCTraj& traj)
{
#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_trajectory_lock);  // aquire lock
#endif //COMMONCPP
  traj = m_last_saved_trajectory;
} // CSimulation::copyLastSavedTrajectory(SCTraj&)


//-----------------------------------------------------------------------------
void CSimulation:: calcLastSavedTrajectoryApprox(void)
{
#if defined COMMONCPP
  SimpleGuard guard1(m_last_saved_clock_angle_approx_values_lock);
  SimpleGuard guard2(m_last_saved_cone_angle_approx_values_lock);
  SimpleGuard guard3(m_last_saved_clock_angle_approx_values_dot_lock);
  SimpleGuard guard4(m_last_saved_cone_angle_approx_values_dot_lock);
#endif //COMMONCPP

  // do calculation if there's no s/c trajectory approximation
  //if (m_last_saved_trajectory_approx.size() == 0)
  //{
  //  calcLastTrajectoryClockAngleApproximation();
  //  calcLastTrajectoryConeAngleApproximation();
    //const unsigned int i(tellp());
    //SafeDelete(m_last_saved_trajectory_approx_fitness[i]);
  //  m_last_saved_trajectory_approx_fitness[i] = evaluateSteeringAngleApprox(
  //    m_last_saved_individual
  //  , m_last_saved_clock_angle_approx_values[i]
  //  , m_last_saved_cone_angle_approx_values[i]
  //  , m_last_clock_angle_approx_time_steps
  //  , m_last_saved_trajectory_approx[i]
  //  , m_last_saved_trajectory_approx_steps[i]
  //  , m_last_saved_start_time);
  //}
}


//-----------------------------------------------------------------------------
void CSimulation::copyLastSavedTrajectoryApprox(SCTraj& traj)
{
#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_trajectory_approx_lock);     // aquire lock
#endif //COMMONCPP
  traj = m_last_saved_trajectory_approx;
}


//-----------------------------------------------------------------------------
void CSimulation::getTrajectoryData(
  const SCTraj& sc_traj
 , CTrajectory& traj
 , const unsigned int sz
 ) const
{
  traj.resize(sz, MathLib::Vector(6));
  if (traj.size() == sz)
  {
    for (uint i=0;i<traj.size();++i) traj[i].Set(sc_traj[i].CState::Cart());
  }
  else Warning("CSimulation::getTrajectoryData: error during resize of target data structure!\n");
} // CSimulation::getTrajectoryData(const SCTraj& , CTrajectory& )


//-----------------------------------------------------------------------------
void CSimulation::getTrajectoryTargetData(const SCTraj& sc_traj
                                         ,const BodyID  target
                                         ,CTrajectory&  traj)
{
  traj.resize(sc_traj.size(), MathLib::Vector(6));
  if (traj.size() == sc_traj.size() && m_solsys != 0)
  {
    for (uint i=0;i<traj.size();++i)
      traj[i].Set(m_solsys->Cart(target, sc_traj[i].MJD()));
  }
} // void CSimulation::getTrajectoryTargetData(...)


//-----------------------------------------------------------------------------
Trajectory CSimulation::getLastTrajectoryData(void)
{
#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_trajectory_lock);      // lock trajectory data
#endif // COMMONCPP
  Trajectory traj = new CTrajectory();
  if (traj != 0)
  {
    getTrajectoryData(m_last_saved_trajectory, *traj
      , m_last_saved_trajectory.size());
    if (traj->empty() == true) SafeDelete(traj);
  }
  else
    Warning("CSimulation::getLastTrajectoryData(void):\
 could not allocate memory for CTrajectory object!\n");
  return traj;
} // CSimulation::getLastTrajectoryData(void)


//-----------------------------------------------------------------------------
Trajectory CSimulation::getLastTrajectoryApproxData(void)
{
#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_trajectory_approx_lock); // aquire lock
#endif //COMMONCPP

  //calcLastSavedTrajectoryApprox(); // calculate new approximation if needed
  Trajectory traj = new CTrajectory();
  if (traj != 0)
  {
    getTrajectoryData(m_last_saved_trajectory_approx/*[tellg()]*/, *traj
      ,m_last_saved_trajectory_approx_steps/*[tellg()]*/);
    if (traj->empty() == true) SafeDelete(traj);
  }
  else
    Warning("CSimulation::getLastTrajectoryApproxData(void):\
 could not allocate memory for CTrajectory object!\n");
  return traj;
} // Trajectory CSimulation::getLastTrajectoryApproxData(void)


//-----------------------------------------------------------------------------
Trajectory CSimulation::getTrajectoryTargetData(void)
{
#if defined COMMONCPP
  SimpleGuard guard(m_sc_traj_lock);                  // lock trajectory data
#endif // COMMONCPP

  BodyID id;
  if (getTrajectoryType() != GTOC2TRAJ)
  {
    id = getTargetBodyId();
  }
  else
  {
    id = m_gtoc2_best_bodies[0];
  }

  Trajectory traj = new CTrajectory();
  if (traj != 0)
  {
    getTrajectoryTargetData(m_spacecraft_trajectory, id, *traj);
    if (traj->empty() == true) SafeDelete(traj);
  }
  else
    Warning("CSimulation::getTrajectoryTargetData(void):\
 could not allocate memory for CTrajectory object!\n");
  return traj;
} // Trajectory CSimulation::getTrajectoryTargetData(void)


//-----------------------------------------------------------------------------
Trajectory CSimulation::getLastTrajectoryTargetData(void)
{
#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_trajectory_lock); // lock last trajectory data
#endif // COMMONCPP

  BodyID id;
  if (getTrajectoryType() != GTOC2TRAJ)
  {
    id = getTargetBodyId();
  }
  else
  {
    id = m_gtoc2_best_bodies[0];
  }

  Trajectory traj = new CTrajectory();
  if (traj != 0)
  {
    getTrajectoryTargetData(m_last_saved_trajectory, id, *traj);
    if (traj->empty() == true) SafeDelete(traj);
  }
  else
  {
    // Warning Message!
  }
  return traj;
} // Trajectory CSimulation::getLastTrajectoryTargetData(void)


//-----------------------------------------------------------------------------
void CSimulation::extractDataFromReferenceTrajectory(void)
{
#if defined COMMONCPP
  SimpleGuard g(m_ref_sc_traj_lock);
#endif //COMMONCPP
  if (m_ref_sc_traj.size() > 0)
  {
    m_ref_time.resize(m_ref_sc_traj.size());
    m_ref_clock_angle.resize(m_ref_sc_traj.size());
    m_ref_cone_angle.resize(m_ref_sc_traj.size());
    for (unsigned int i = 0;i<m_ref_sc_traj.size();++i)
    {
      m_ref_time[i]=m_ref_sc_traj[i].MJD()-m_ref_sc_traj[0].MJD();
      m_ref_clock_angle[i]=m_ref_sc_traj[i].getClockAngle();
      m_ref_cone_angle[i]=m_ref_sc_traj[i].getConeAngle();
    }
  }
} // void CSimulation::extractDataFromReferenceTrajectory(void)

//-----------------------------------------------------------------------------
Trajectory CSimulation::getReferenceTrajectoryTargetData(void)
{
#if defined COMMONCPP
  SimpleGuard guard(m_ref_sc_traj_lock); // lock ref trajectory data
#endif // COMMONCPP
  Trajectory traj = new CTrajectory();
  if (traj != 0)
  {
    getTrajectoryTargetData(m_ref_sc_traj, getTargetBodyId(), *traj);
    if (traj->empty() == true) SafeDelete(traj);
  }
  else
  {
    // Warning Message!
  }
  return traj;
} // Trajectory CSimulation::getReferenceTrajectoryTargetData(void)


//-----------------------------------------------------------------------------
int CSimulation::getBasicClockAngleData(
      const    SCTraj& traj
     ,double*& data_array
     ,const    int n_max_points
    )
{
  double* angle_history = 0;
  int array_size = (n_max_points > 0 ? n_max_points:traj.size());

  data_array = 0;
  if (array_size  > 0 ) {

    angle_history = new double[array_size];

    if (angle_history != 0)  // if allocation did not fail
    {
      double max_val = 0.0;  // maximum and minimum values are not used right
      double min_val = 0.0;  // now, but may be in the future to provide
      // users of the returned array with the information of maximal and
      // minimal values in array to that they don't need to find out for
      // themselves

      double sum0 = 0.0;    // the the summation variables are needed to
      double sumPI = 0.0;   // determine, whether more points are concentrated
      // around -PI or zero to possibly transform data to [0;2PI) scope
      register double temp;
      register double temp2;

      for (int i=0;i < array_size; ++i)
      {
        temp = (traj[i]).getClockAngle();

        angle_history[i] = temp;           // set array value

        temp2 = temp > 0 ?temp:temp+TWOPI; // for maximum/minimum determination

        sum0  +=(temp=fabs(temp));
        sumPI +=PI-temp;

        max_val = std::max(max_val,temp2);
        min_val = std::min(min_val,temp2);
      }

      if (sumPI < sum0)   // is most of the data values in the vincinity
      {                   // of +/- PI ?
        //If true, negativ data must be transformed to prevent jumping of
        // values from -(PI+eps) to +(PI+eps) in graphical output

        for (int i=0;i<array_size;++i)
          if (angle_history[i] <= 0.0)
            angle_history[i]+=TWOPI;

      } // if (sumPI < sum0)

    } // if (angle_history != 0)
    else array_size = 0;

  } // if ((getCurrentStepCount() > 0) && (m_spacecraft_trajectory != 0))
  data_array = angle_history;
  return array_size;
} // int CSimulation::getBasicClockAngleData(...)

//-----------------------------------------------------------------------------
void CSimulation::getClockAngleData(
       const SCTraj&     traj
      ,      CDataArray& data
      ,const int max_n_items)
{
  double* angle_history = 0;  // double precision array of angle values

  if (m_clock_angle_history == true ||
      m_clock_angle_history_approximation == true) {

    // get angle history from simulation object
    int array_size = getBasicClockAngleData(
                       traj
                      ,angle_history, max_n_items
                     );

    if (angle_history != 0 && array_size > 0)
    {
      data.resize(array_size);

      if (data.size() == array_size) {

        // copy data values
        for (int i = 0; i < array_size; ++i)
          data[i] = angle_history[i];

      }
      delete[] angle_history;  // free temporary array of double values
    }
  } //if (m_clock_angleHistory == true)
} // void CSimulation::getClockAngleData(...)


//-----------------------------------------------------------------------------
DataArray CSimulation::getClockAngleData(void)
{
#if defined COMMONCPP
  SimpleGuard guard(m_sc_traj_lock);  // lock trajectory data
#endif // COMMONCPP

  DataArray data = new CDataArray();
  if (data != 0) {
    getClockAngleData(m_spacecraft_trajectory, *data, getCurrentStepCount());
    if (data->size() == 0) SafeDelete(data);
  }
  return data;
} // DataArray CSimulation::getClockAngleData(void)


//-----------------------------------------------------------------------------
DataArray CSimulation::getLastTrajectoryClockAngleData(void)
{
#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_trajectory_lock);
#endif // COMMONCPP
  DataArray data_array = 0;

  if (m_clock_angle_history == true)
  {
    if (m_last_saved_trajectory.empty() != true)
    {
      if (m_last_saved_clock_angle_history.empty() == true)
        getClockAngleData(m_last_saved_trajectory
          ,m_last_saved_clock_angle_history
        );

      if (m_last_saved_clock_angle_history.empty() != true)
      {
        data_array = new CDataArray(m_last_saved_clock_angle_history);

        if (data_array == 0) {
          ostringstream os;
          os << "CSimulation::getClockAngleHistory(void):(" << __FILE__;
          os << ":" << __LINE__ << ") could not allcate memory for clock \
            angle data!\n";
          Warning(os.str());
        }
      } else Warning("CSimulation::getLastTrajectoryClockAngleData(void): \
 no clock angle data could be obtained from last saved trajectory!\n");
    }
    else
    {
      ostringstream os;
      os << "CSimulation::getClockAngleHistory(void):(" << __FILE__;
      os << ":" << __LINE__ << ") Last saved trajectory is empty!\n";
      Warning(os.str());
    }
  } // if (m_clock_angle_history == true)
  return data_array;
} // DataArray CSimulation::getLastTrajectoryClockAngleData(void)


//-----------------------------------------------------------------------------
DataArray CSimulation::getLastTrajectoryThrustData(void)
{
  const std::string method_name("DataArray CSimulation::getLastTrajectory\
ThrustData(void)");

#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_trajectory_lock);
#endif // COMMONCPP
  DataArray data_array = 0;
  if (m_thrust_history == true)
  {
    if (m_last_saved_trajectory.empty() != true)
    {
      if (m_last_saved_thrust_history.empty() == true)
        getBasicThrustData(m_last_saved_trajectory,
          m_last_saved_thrust_history);

      if (m_last_saved_thrust_history.empty() != true)
      {
        data_array = new CDataArray(m_last_saved_thrust_history);

        if (data_array == 0) {
          ostringstream os;
          os << method_name << ":(" << __FILE__ <<":" << __LINE__
             << ") could not allcate memory for clock angle data!\n";
          Warning(os.str());
        }
      } else Warning(method_name + std::string(": no thrust data could be obtained from \
last saved trajectory!\n"));
    }
    else
    {
      ostringstream os;
      os << method_name + ":(" << __FILE__ << ":" << __LINE__
         << ") Last saved trajectory is empty!\n";
      Warning(os.str());
    }
  } // if (m_thrust_history == true)
  return data_array;
} //DataArray CSimulation::getLastTrajectoryThrustData(void)


//-----------------------------------------------------------------------------
DataArray CSimulation::getReferenceTrajectoryClockAngleData(void)
{
#if defined COMMONCPP
  SimpleGuard guard(m_ref_sc_traj_lock);
#endif // defined COMMONCPP

  DataArray data_array = 0;

  if (m_ref_sc_traj.empty() != true)
  {

    data_array = new CDataArray();
    if (data_array != 0)
    {
      getBasicClockAngleData(m_ref_sc_traj,*data_array);
      postprocessAngleData(*data_array);
    } else
    {
      ostringstream os;
      os << "CSimulation::getReferenceTrajectoryClockAngleData(void):(" << __FILE__;
      os << ":" << __LINE__ << ") could not allcate memory for clock \
        angle data!\n";
      Warning(os.str());
    }
  }
  else
  {
    ostringstream os;
    os << "CSimulation::getReferenceTrajectoryClockAngleData(void):(" << __FILE__;
    os << ":" << __LINE__ << ") Reference trajectory is empty!\n";
    Warning(os.str());
  }
  return data_array;
} // DataArray CSimulation::getReferenceTrajectoryClockAngleData(void)


//-----------------------------------------------------------------------------
DataArray CSimulation::getReferenceTrajectoryConeAngleData(void)
{
#if defined COMMONCPP
  SimpleGuard guard(m_ref_sc_traj_lock);
#endif // defined COMMONCPP

  DataArray data_array = 0;

  if (m_ref_sc_traj.empty() != true)
  {

    data_array = new CDataArray();
    if (data_array != 0)
    {
      getBasicConeAngleData(m_ref_sc_traj,*data_array);
    } else
    {
      ostringstream os;
      os << "CSimulation::getReferenceTrajectoryConeAngleData(void):(" << __FILE__;
      os << ":" << __LINE__ << ") could not allcate memory for clock \
        angle data!\n";
      Warning(os.str());
    }
  }
  else
  {
    ostringstream os;
    os << "CSimulation::getReferenceTrajectoryConeAngleData(void):(" << __FILE__;
    os << ":" << __LINE__ << ") Reference trajectory is empty!\n";
    Warning(os.str());
  }
  return data_array;
} // DataArray CSimulation::getReferenceTrajectoryConeAngleData(void)


// returns a copy of approximated clock angle values
//-----------------------------------------------------------------------------
DataArray CSimulation::getClockAngleApproximation(void)
{

  CDataArray original_data;
  DataArray  approx_data = 0;
  CDataArray approx_data_d;

#if defined COMMONCPP
  SimpleGuard guard(m_sc_traj_lock);  // lock trajectory data
#endif // COMMONCPP

  if (m_spacecraft_trajectory.empty() != true)
  {

    getClockAngleData(m_spacecraft_trajectory, original_data, getCurrentStepCount());

    if (original_data.empty() != true)
    {
      approx_data   = new CDataArray();

      if (approx_data != 0 )
      {
        calcClockAngleApproximation(
          original_data
         ,m_last_saved_time_steps
         ,m_last_clock_angle_approx_time_steps
         ,*approx_data
         ,approx_data_d
         ,m_clock_angle_approx_poly_degree
         ,m_clock_angle_approx_seg_size
         );

      }
      else
      {
        Warning(string("CSimulation::getClockAngleApproximation(void): \
 Allocation of data array for approximated angle values failed!\n"));
      }

    }
  } //  if (m_spacecraft_trajectory.empty() != 0)
  return approx_data;
}; // DataArray CSimulation::getClockAngleApproximation(void)


//-----------------------------------------------------------------------------
DataArray CSimulation::getReferenceTrajectoryThrustData(void)
{
  const std::string method_name("CSimulation::getReferenceTrajectoryConeAngle\
Data(void)");

#if defined COMMONCPP
  SimpleGuard guard(m_ref_sc_traj_lock);
#endif // defined COMMONCPP

  DataArray data_array = 0;

  if (m_ref_sc_traj.empty() != true)
  {

    data_array = new CDataArray();
    if (data_array != 0)
    {
      getBasicThrustData(m_ref_sc_traj,*data_array);
    } else
    {
      ostringstream os;
      os << method_name << ":(" << __FILE__ << ":" << __LINE__
         << ") could not allcate memory for thrust data!\n";
      Warning(os.str());
    }
  }
  else
  {
    ostringstream os;
    os << method_name << ":(" << __FILE__ << ":" << __LINE__
       << ") Reference trajectory is empty!\n";
    Warning(os.str());
  }
  return data_array;
} // DataArray CSimulation::getReferenceTrajectoryThrustData(void)


//------------------------------------------------------------------------------
Trajectory CSimulation::getReferenceTrajectoryData(void)
{
#if defined COMMONCPP
  SimpleGuard g(m_ref_sc_traj_lock);      // lock trajectory data
#endif // COMMONCPP
  return getTrajectoryData(m_ref_sc_traj);
} // Trajectory CSimulation::getReferenceTrajectoryData(void)


//------------------------------------------------------------------------------
void CSimulation::calcLastTrajectoryClockAngleApproximation(void)
{
  if (m_clock_angle_history_approximation == true)
  {
    // check for available clock angle data and get it if possible
    if (m_last_saved_trajectory.empty() != true &&
        m_last_saved_clock_angle_history.empty() == true)
    {
      // get clock angle values from trajectory
      getClockAngleData(
         m_last_saved_trajectory
        ,m_last_saved_clock_angle_history
      );
    }

    // check for available clock angle approximation data and get it if possible
    if (m_last_saved_clock_angle_history.empty() != true &&
        m_last_saved_clock_angle_approx_values.empty() == true)
    {
//      const unsigned int i(tellg());
      // get clock angle approximation
      calcClockAngleApproximation(
        m_last_saved_clock_angle_history
       ,m_last_saved_time_steps
       ,m_last_clock_angle_approx_time_steps
       ,m_last_saved_clock_angle_approx_values/*[i]*/
       ,m_last_saved_clock_angle_approx_values_dot/*[i]*/
       ,m_clock_angle_approx_poly_degree
       ,m_clock_angle_approx_seg_size
      );
    }
  } // if (m_clock_angle_history_approximation == true)
} // void CSimulation::calcLastTrajectoryClockAngleApproximation(void)

//-----------------------------------------------------------------------------
DataArray CSimulation::getLastTrajectoryClockAngleApproximation(void)
{
#if defined COMMONCPP
  SimpleGuard guard1(m_last_saved_trajectory_lock);
  SimpleGuard guard2(m_last_saved_clock_angle_history_lock);
  SimpleGuard guard3(m_last_saved_clock_angle_approx_values_lock);
  SimpleGuard guard4(m_last_saved_clock_angle_approx_values_dot_lock);
#endif //COMMONCPP

  DataArray  approx_data = 0;
  /*calcLastTrajectoryClockAngleApproximation();*/
  if (m_last_saved_clock_angle_approx_values/*[tellg()]*/.empty() != true)
  {
    approx_data = new CDataArray(m_last_saved_clock_angle_approx_values/*[tellg()]*/);
  } else Warning("CSimulation::getLastTrajectoryClockAngleApproximation(void):\
 data structure of last clock angle approximation is empty!\n");

  return approx_data;
} // DataArray CSimulation::getLastTrajectoryClockAngleApproximation(void)


//-----------------------------------------------------------------------------
void CSimulation::calcApproximation(
  const CDataArray& data
 ,const CDataArray& time
 ,const unsigned int degree
 ,const unsigned int seg_size
 ,const double resolution_factor
 ,      CDataArray& approx_time
 ,      CDataArray& approx_data
 ,      CDataArray& approx_data_derivative
)
{
  if (data.size()==0 || data.size()!=time.size()) return;
  valarray<double> input(data.size());      CopyDataStructure(data, input);
  valarray<double> input_time(time.size()); CopyDataStructure(time, input_time);
  valarray<unsigned int> segmentation,degrees;
  SetUpApproximationArrays(input, segmentation, degrees, seg_size, degree);

  vector<Polynomial> polynomials;
  const PolyApproxType poly_type(ZERO_START);
  valarray<double> v_approx_data;

  CalcPolyApprox(input_time, input, segmentation, degrees, polynomials, v_approx_data);
  //CalculateApproximationPolynomials(input, input_time, segmentation, degrees,
  //  polynomials, poly_type);

  valarray<double> v_approx_time;
  valarray<double> v_approx_data_d;

  //CalcPolyApproxValues(input_time,segmentation, polynomials, v_approx_time
  //  , v_approx_data, resolution_factor, poly_type);

  CalcPolyApproxData(polynomials, input_time, segmentation
    , resolution_factor, v_approx_time, v_approx_data);

  //CalcPolyApproxFirstDerivative(input_time, segmentation, polynomials, v_approx_time
  //  , v_approx_data_d, resolution_factor, poly_type);

  CalcPolyApproxFirstDerivative(polynomials, input_time, segmentation
    , resolution_factor, v_approx_time, v_approx_data_d);

  if (v_approx_data.size() != v_approx_data_d.size()) return;
  CopyDataStructure(v_approx_time, approx_time);
  CopyDataStructure(v_approx_data, approx_data);
  CopyDataStructure(v_approx_data_d, approx_data_derivative);

} // CSimulation::calcApproximation


#undef NEW_CODE
#define NEW_CODE
//-----------------------------------------------------------------------------
void CSimulation::calcClockAngleApproximation(
    const CDataArray& data                  // array of original values
   ,const CDataArray& time                  // array of time steps
   ,      CDataArray& new_time              // array of new time steps
   ,      CDataArray& new_data              // array of new values
   ,      CDataArray& new_data_d            // array of new values first der.
   ,const int         degree
   ,const int         n_seg_points)
{

  if (data.empty() != true)
  {

    bool m_clock_angle_history_near_pi = false;

    for (uint i = 0; i < data.size();++i)
    {
      if (data[i] > PI)
      {
        m_clock_angle_history_near_pi = true;
        break;
      }
    }
#if defined NEW_CODE

    try
    {
      calcApproximation(data, time, degree, n_seg_points, 2.0
        , new_time, new_data, new_data_d);
    }
    catch (std::exception& exc)
    {
      Output(exc.what());
    }

#else
      int     n_segments    = 0;
      double* poly_params   = 0;

      try
      {
        poly_params =
          CalculateApproximationPolynomials(degree
                                          ,data
                                          ,data.size()
                                          ,n_segments
                                          ,n_seg_points);
      }
      catch (exception& exc)
      {
        Output(exc.what());
      }


      if (poly_params != 0)
      {
        new_data.resize(data.size());
        new_data_d.resize(data.size());
        if (new_data.size()   == data.size() &&
            new_data_d.size() == data.size())
        {
          int rest_points = data.size();

          for (int segment = 0;segment < n_segments; ++segment) // segment loop
          {
            const int n_poly_params = degree+1;
            Polynomial polynom(degree, poly_params+(segment*n_poly_params));

            const int points = ((rest_points < (2.0*n_seg_points-1)) && rest_points > n_seg_points)?
                                rest_points:
                                n_seg_points;
            double x_inc = 1.0/(points-1);
            double x=0;

            register double angle;
            if (m_clock_angle_history_near_pi == true)
            {
              for (int i=0;i<points;++i,x+=x_inc)
              {
                angle = polynom.getValue(x);
                angle = bounded(angle, 0.0, TWOPI);
                new_data[getCurrentStepCount()-rest_points+i] = angle;
                new_data_d[data.size()-rest_points+i] =
                  polynom.firstDerivative(x);
              }
            } // if (m_clock_angle_history_near_pi == true)
            else
            {
              for (int i=0;i<points;++i,x+=x_inc)
              {
                angle = polynom.getValue(x);
                angle = bounded(angle, -PI,PI);
                new_data[data.size()-rest_points+i] = angle;
                new_data_d[data.size()-rest_points+i] =
                  polynom.firstDerivative(x);
              }
            }

            rest_points-=(points-1);
          } // for(int i = 0;i < segments;++i) // segment loop

        } // if (new_data.size()   == data.size() &&
        delete[] poly_params;  // free temporary array of approx. polynom params
      } // if (poly_params != 0)
#endif
  } //if (data.empty() != true)
  else Warning("CSimulation::calcClockAngleApproximation(...): input data array is empty!\n");
#undef NEW_CODE
} // CSimulation::calcClockAngleApproximation


//-----------------------------------------------------------------------------
int CSimulation::getBasicConeAngleData(
  const SCTraj&  traj
 ,      double*& data_array
 ,      const    int n_max_points
 )
{
  double* angle_history = 0;
  int array_size = (n_max_points > 0 ? n_max_points:traj.size());

  data_array = 0;
  if (array_size > 0 )
  {

    angle_history = new double[array_size];

    if (angle_history != 0)  // if allocation did not fail
    {
      for (int i=0;i < array_size; ++i)
      {
        angle_history[i] = (traj[i]).getConeAngle();
      }
    } // if (angle_history != 0)
    else array_size = 0;

  } // if (array_size > 0))
  data_array = angle_history;
  return array_size;
} // int CSimulation::getBasicConeAngleData(


// returns a copy of approximated cone angle values
//-----------------------------------------------------------------------------
DataArray CSimulation::getConeAngleApproximation(void)
{
  CDataArray original_data;
  DataArray  approx_data = 0;
  CDataArray approx_data_d;

#if defined COMMONCPP
  SimpleGuard guard(m_sc_traj_lock);  // lock trajectory data
#endif // COMMONCPP

  double upper_angle_bound(0.0);

  if (m_spacecraft_trajectory.empty() != true)
  {

    if (m_spacecraft_trajectory[0].isSolarSail() == true)
    {
      upper_angle_bound = HALFPI-HALFPI_MODIFICATION;
    }
    else
    {
      upper_angle_bound = PI;
    }

    getConeAngleData(m_spacecraft_trajectory, original_data, getCurrentStepCount());

    if (original_data.empty() != true)
    {
      approx_data = new CDataArray();
      if (approx_data != 0)
      {
        calcConeAngleApproximation(
          original_data
         ,m_last_saved_time_steps
         ,m_last_cone_angle_approx_time_steps
         ,*approx_data
         ,approx_data_d
         ,m_cone_angle_approx_poly_degree
         ,m_cone_angle_approx_seg_size
         ,upper_angle_bound
         );

      }
      else
      {
        Warning(string("CSimulation::getConeAngleApproximation(void): \
 Allocation of data array for approximated angle values failed!\n"));
      }

    }
  } //  if (m_spacecraft_trajectory.empty() != 0)
  return approx_data;
} // DataArray CSimulation::getConeAngleApproximation(void)


#undef NEW_CODE
#define NEW_CODE
//-----------------------------------------------------------------------------
void CSimulation::calcConeAngleApproximation(
    const CDataArray& data
   ,const CDataArray& time
   ,      CDataArray& new_time
   ,      CDataArray& new_data
   ,      CDataArray& new_data_d            // array of new values first der.
   ,const int         degree
   ,const int         n_seg_points
   ,const double      upper_bound)
{

  if (data.empty() != true)
  {

#if defined NEW_CODE

    try
    {
      calcApproximation(data, time, degree, n_seg_points, 2.0, new_time,
        new_data, new_data_d);
    }
    catch (exception& exc)
    {
      Output(exc.what());
    }

#else
      int     n_segments    = 0;
      double* poly_params   = 0;

      try
      {
        poly_params =
          CalculateApproximationPolynomials(degree
                                          ,data
                                          ,data.size()
                                          ,n_segments
                                          ,n_seg_points);
      }
      catch (exception& exc)
      {
        Output(exc.what());
      }

      if (poly_params != 0)
      {
        new_data.resize(data.size());
        new_data_d.resize(data.size());
        if (new_data.size()   == data.size() &&
            new_data_d.size() == data.size())
        {
          int rest_points = data.size();

          for (int segment = 0;segment < n_segments; ++segment) // segment loop
          {
            const int n_poly_params = degree+1;
            Polynomial polynom(degree, poly_params+(segment*n_poly_params));

            const int points = ((rest_points < (2.0*n_seg_points-1)) && rest_points > n_seg_points)?
                                rest_points:
                                n_seg_points;
            double x_inc = 1.0/(points-1);
            double x=0;

            register double angle;
            for (int i=0;i<points;++i,x+=x_inc)
            {
              angle = bounded(polynom.getValue(x),0.0,upper_bound);
              angle = bounded(angle, 0.0, TWOPI);
              new_data[data.size()-rest_points+i] = angle;
              new_data_d[data.size()-rest_points+i] =
                polynom.firstDerivative(x);
            }

            rest_points-=(points-1);
          } // for(int i = 0;i < segments;++i) // segment loop

        } // if (new_data.size() > 0)
        delete[] poly_params;  // free temporary array of approx. polynom params
      } // if (poly_params != 0)
#endif // NEW_CODE
  } //if (data.empty() != true)

#undef NEW_CODE
}


//-----------------------------------------------------------------------------
void CSimulation::getConeAngleData(
  const SCTraj&     traj
 ,      CDataArray& data
 ,const int         max_n_items)
{
  double* angle_history = 0;  // double precision array of angle values

  if (m_cone_angle_history == true  ||
      m_cone_angle_history_approximation == true) {

    // get angle history from simulation object
    int array_size = getBasicConeAngleData(traj, angle_history, max_n_items);

    if (angle_history != 0 && array_size > 0)
    {
      data.resize(array_size);

      if (data.size() == array_size) {

        // copy data values
        for (int i = 0; i < array_size; ++i)
          data[i] = angle_history[i];

      }
      delete[] angle_history;  // free temporary array of double values
    }
  } //if (m_clock_angleHistory == true)
} // void CSimulation::getConeAngleData(...)


//-----------------------------------------------------------------------------
#if defined COMMONCPP
DataArray CSimulation::getConeAngleData(void)
{

  SimpleGuard guard(m_sc_traj_lock);  // lock trajectory data

  DataArray java_data = new CDataArray();
  if (java_data != 0) 
  {
    getConeAngleData(m_spacecraft_trajectory,*java_data,getCurrentStepCount());
    if (java_data->empty() == true)
    {
      delete java_data;
      java_data = 0;
    }
  }
  return java_data;

} // DataArray CSimulation::getConeAngleData(void)
#endif


//-----------------------------------------------------------------------------
DataArray CSimulation::getLastTrajectoryConeAngleData(void)
{
#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_trajectory_lock);
#endif // COMMONCPP

  DataArray data_array = 0;

  if (m_last_saved_trajectory.empty() != true)
  {

    if (m_last_saved_cone_angle_history.empty() == true)
      getConeAngleData(
         m_last_saved_trajectory
        ,m_last_saved_cone_angle_history
      );

    if (m_last_saved_cone_angle_history.empty() != true)
    {
      data_array = new CDataArray(m_last_saved_cone_angle_history);

      if (data_array == 0) {
        ostringstream os;
        os << "CSimulation::getLastTrajectoryConeAngleData(void):(" << __FILE__;
        os << ":" << __LINE__ << ") could not allcate memory for clock \
          angle data!\n";
        Warning(os.str());
      }
    }
  }
  else
  {
    ostringstream os;
    os << "CSimulation::getLastTrajectoryConeAngleData(void):(" << __FILE__;
    os << ":" << __LINE__ << ") Last saved trajectory is empty!\n";
    Warning(os.str());
  }
  return data_array;
} // DataArray CSimulation::getLastTrajectoryConeAngleData(void)



//-----------------------------------------------------------------------------
void CSimulation::calcLastTrajectoryConeAngleApproximation(void)
{
  if (m_cone_angle_history_approximation == true)
  {
    // check for available cone angle data and get it if possible
    if (m_last_saved_trajectory.empty() != true &&
        m_last_saved_cone_angle_history.empty() == true)
    {
      // get clock angle values from trajectory
      getConeAngleData(
         m_last_saved_trajectory
        ,m_last_saved_cone_angle_history
      );
    }

    // check for available clock angle approximation data and get it if possible
    if (m_last_saved_cone_angle_history.empty() != true &&
        m_last_saved_cone_angle_approx_values.empty() == true)
    {
      double upper_angle_bound(0.0);
      if (m_last_saved_trajectory[0].isSolarSail() == true)
      {
        upper_angle_bound = HALFPI-HALFPI_MODIFICATION;
      }
      else
      {
        upper_angle_bound = PI;
      }

      // get clock angle approximation
//      const unsigned int i(tellp());
      calcConeAngleApproximation(
        m_last_saved_cone_angle_history
       ,m_last_saved_time_steps
       ,m_last_cone_angle_approx_time_steps
       ,m_last_saved_cone_angle_approx_values/*[i]*/
       ,m_last_saved_cone_angle_approx_values_dot/*[i]*/
       ,m_cone_angle_approx_poly_degree
       ,m_cone_angle_approx_seg_size
       ,upper_angle_bound
      );
    }
  }
} // CSimulation::calcLastTrajectoryConeAngleApproximation


//-----------------------------------------------------------------------------
DataArray CSimulation::getLastTrajectoryConeAngleApproximation(void)
{
#if defined COMMONCPP
  SimpleGuard guard1(m_last_saved_trajectory_lock);
  SimpleGuard guard2(m_last_saved_cone_angle_history_lock);
  SimpleGuard guard3(m_last_saved_cone_angle_approx_values_lock);
  SimpleGuard guard4(m_last_saved_cone_angle_approx_values_dot_lock);
#endif // COMMONCPP
  DataArray  approx_data = 0;
  //calcLastTrajectoryConeAngleApproximation();

  //const unsigned int i(tellg());
  if (m_last_saved_cone_angle_approx_values/*[i]*/.empty() != true)
  {
    approx_data = new CDataArray(m_last_saved_cone_angle_approx_values/*[i]*/);
  }

  return approx_data;
} // DataArray CSimulation::getLastTrajectoryConeAngleApproximation(void)


//-----------------------------------------------------------------------------
void CSimulation::initDESSolver( CSC&    sc
                                ,DESSolverType solver_type
                                ,Vector& ys
                                ,const   double t
                                ,RK4*    integrator_RK4
                                ,RKF54*  integrator_RKF54
                                ,DOPRI8* integrator_DOPRI8
                                ,DE*     integrator_DESG
                                )
{
   switch (solver_type) 
   {
     case _RK4:
       ys=sc.PolarStatePropellantMassSRDVector();
       break;

     case _RKF54:
       ys=sc.PolarStatePropellantMassSRDVector();
       if (integrator_RKF54 != 0)
       {
         integrator_RKF54->init(N_MAX_ITGR_CALLS, getMaxRelIntegrationError(), getMaxAbsIntegrationError());
         integrator_RKF54->resetStopFunc();
       } // if (integrator_RKF54 != 0)
       break;

     case _DOPRI8:
       ys=sc.PolarStatePropellantMassSRDVector();
       if (integrator_DOPRI8 != 0)
         integrator_DOPRI8->Init(N_MAX_ITGR_CALLS, getMaxRelIntegrationError(), getMaxAbsIntegrationError());
       break;

     case _DESG:
       ys=sc.PolarStatePropellantMassSRDVector();
       if (integrator_DESG != 0)
         integrator_DESG->Init(t, getMaxRelIntegrationError(), getMaxAbsIntegrationError());
       break;

     default:
       //FatalError("CSimulation::InitDESSolver(): unknown spacecraft type");
       break;
   }

} // end CSimulation::InitDESSolver

//-----------------------------------------------------------------------------
void CSimulation::initDESSolver( Vector&ys
                                ,const double t
                                )
{

  initDESSolver(
    m_spacecraft
   ,getIntegratorType()
   ,ys
   ,t
   ,m_SpacecraftOrbit_RK4
   ,m_SpacecraftOrbit_RKF54
   ,m_SpacecraftOrbit_DOPRI8
   ,m_SpacecraftOrbit_DESG);

} // end CSimulation::InitDESSolver

// ----------------------------------------------------------------------------
RK_STATE CSimulation::integrate(
  CSC&          sc
 ,Vector &      ys
 ,double &      t
 ,const double  h
 ,DESSolverType des_solver
 ,RK4*          integr_RK4
 ,RKF54*        integr_RKF54
 ,DOPRI8*       integr_DOPRI8
 ,DE*           integr_DESG
 )
{
  RK_STATE rks = RK_DONE;
  switch (des_solver)
  {
    case _RK4:
      integr_RK4->Step(t, ys, h);
      break;
    case _RKF54:
      rks = integr_RKF54->Step(t, ys, h);
      break;
    case _DOPRI8:
      rks = integr_DOPRI8->Step(t, ys, h);
      break;
    case _DESG:
      integr_DESG->Integ(t+=h, ys);
      break;
    default:
      ostringstream os;
      os << __FILE__ << ", L"<<__LINE__<<" Unknown DES solver!\n";
      FatalError(os.str());
      return RK_DONE; // dummy
  }
  sc.setPropellantMassSRDVector(ys, t);
  return rks;
} // end CSimulation::Integrate



// ----------------------------------------------------------------------------
RK_STATE CSimulation::integrate(Vector & ys, double & t, const double h)
// ----------------------------------------------------------------------------
{
  return integrate(m_spacecraft,ys,t,h,getIntegratorType(), m_SpacecraftOrbit_RK4,
    m_SpacecraftOrbit_RKF54,m_SpacecraftOrbit_DOPRI8,m_SpacecraftOrbit_DESG);

} // end CSimulation::Integrate


// ----------------------------------------------------------------------------
void CSimulation::integrationPreProcessing(void)
{
  switch (m_sim_type)
  {
  case STANDARD:

    switch (tt)
    {
    case RENDEZVOUS:

      switch (getTargetStateType())
      {
      case ORBITAL_ELEMENTS_STATE:
        integrationPreProcessingTOES();
        break;
      default:break;
      } // switch (getTargetStateType())

      break;
    default:break;
    } // switch (tt)

    break;

  default:
    break;
  } // switch (m_sim_type)
} // void CSimulation::integrationPreProcessing(void)


// ----------------------------------------------------------------------------
void CSimulation::integrationPreProcessingTOES(void)
{
  const double min_dev = 1.0e-10;
  m_initial_proximity = 0.0;
  double current_val;

  for (unsigned i = SEMIMAJORAXIS; i<=PERIAPSIS; ++i)
  {
    const KEPLERIAN_ELEMENT el = KEPLERIAN_ELEMENT(i);

    if (useTargetOrbitalElement(el) == true)
    {

      const double max_delta = m_data.getTargetOrbitalElementMaxDeviation(el);
      const double target_val = m_data.getTargetOrbitalElement(el);

      switch (el)
      {
      case SEMIMAJORAXIS:current_val=getCurrentSC().CState::a();break;
      case ECCENTRICITY:current_val=getCurrentSC().CState::e();break;
      case INCLINATION:current_val=getCurrentSC().CState::i();break;
      case NODE:current_val=getCurrentSC().CState::Omega();break;
      case PERIAPSIS:current_val=getCurrentSC().CState::omega();break;
      case ANOMALY:current_val=getCurrentSC().CState::M();break;
      default:break;
      }

      getTOEInitialElementDeviation()[i] = std::max(min_dev,
        fabs(CalcScalingFactor(el, current_val, target_val)*
        CalcDistanceFunction(el, current_val, target_val, max_delta)));

      //initial_elements_deviation[i] =
      //  Global::max(1.0e-10, fabs(
      //  CalcScalingFactor(el, m_data.getInitialOrbitElement(el), m_data.getTargetOrbitalElement(el))
      //  *CalcDistanceFunction(el
      //  , m_data.getInitialOrbitElement(el)
      //  , m_data.getTargetOrbitalElement(el)
      //  , m_data.getTargetOrbitalElementMaxDeviation(el))));

      m_initial_proximity +=
        (getTOEInitialElementDeviation()[el]*
        ((getTargetOrbitalElementProxWeights())[el]));
    }
  }
  m_initial_proximity/=(getTargetOrbitalElementProxWeights()).sum();

  if (m_best_initial_proximity>m_initial_proximity)
  {
    m_best_initial_proximity = m_initial_proximity;
  }

} // void CSimulation::integrationPreProcessingTOES(void)


// ----------------------------------------------------------------------------
void CSimulation::initSensitivityFactors(void)
{
  m_sfa=m_solsys->SensitivityFactor(getTargetBodyId(), 0, m_max_dr_km,
   m_max_dv_kms);
  m_sfe=m_solsys->SensitivityFactor(getTargetBodyId(), 1, m_max_dr_km,
   m_max_dv_kms);
  m_sfi=m_solsys->SensitivityFactor(getTargetBodyId(), 2, m_max_dr_km,
   m_max_dv_kms);
  m_sfh=10.0*m_sfi;

  if (m_max_dv_kms>=1.0) // arbitrary rule of thumb
  {
    m_sfh/=1000.0; // arbitrary rule of thumb
    m_sfa/=1000.0; // arbitrary rule of thumb
    m_sfe/=1000.0; // arbitrary rule of thumb
  }
}; // end CSimulation::InitSensitivityFactors


// ----------------------------------------------------------------------------
void CSimulation::calcBestValuesGTOC2(void)
{
  m_gtoc2_best_bodies = m_gtoc2_temp_best_bodies;
  m_best_dr=getCurrProximity().GetNormedDistanceProximity(getMaxRelDistance());
  m_best_dv=getCurrProximity().GetNormedVelocityProximity(getMaxRelVelocity());
  m_best_mp = calcPropMassFitness(
    getCurrentSC().getInitialPropellantMass(),
    getCurrentSC().getPropellantMass());
  m_best_prox_sum = m_temp_prox_sum;
} // void CSimulation::calcBestValuesTest0(void)

// ----------------------------------------------------------------------------
void CSimulation::calcBestValuesTest1(void)
{
} // void CSimulation::calcBestValuesTest0(void)

// ----------------------------------------------------------------------------
void CSimulation::calcBestValuesTest2(void)
{
} // void CSimulation::calcBestValuesTest0(void)

// ----------------------------------------------------------------------------
void CSimulation::calcBestValuesTest3(void)
{
} // void CSimulation::calcBestValuesTest0(void)

// ----------------------------------------------------------------------------
void CSimulation::calcBestValuesTest4(void)
{
} // void CSimulation::calcBestValuesTest0(void)

// ----------------------------------------------------------------------------
void CSimulation::calcBestValuesTest5(void)
{
} // void CSimulation::calcBestValuesTest0(void)

// ----------------------------------------------------------------------------
void CSimulation::calcBestValuesTest6(void)
{
} // void CSimulation::calcBestValuesTest0(void)

// ----------------------------------------------------------------------------
void CSimulation::calcBestValuesTest7(void)
{
} // void CSimulation::calcBestValuesTest0(void)

// ----------------------------------------------------------------------------
void CSimulation::calcBestValuesTest8(void)
{
} // void CSimulation::calcBestValuesTest0(void)

// ----------------------------------------------------------------------------
void CSimulation::calcBestValuesTest9(void)
{
} // void CSimulation::calcBestValuesTest0(void)


// ----------------------------------------------------------------------------
void CSimulation::calcBestValuesRendezvous(void)
{
  double de,dh;

  switch (getTargetStateType())
  {
  case ORBITAL_ELEMENTS_STATE:
    m_best_dr = m_best_dv = 1.0;
    setTOEBestElementDeviation(getTOEElementDeviation());
    break;
  default:
    m_best_dr = getCurrProximity().GetNormedDistanceProximity(m_max_dr_AU);
    m_best_dv = getCurrProximity().GetNormedVelocityProximity(m_max_dv_AUd);
    m_best_da = m_sfa *
      fabs(m_spacecraft.CState::a() - m_solsys->a(getTargetBodyId(), m_spacecraft.MJD()));

    m_best_de = m_sfe *
      fabs(m_spacecraft.CState::e()-m_solsys->e(getTargetBodyId(), m_spacecraft.MJD()));

    de = fabs(1.0-Dot(m_spacecraft.CState::PVector()
      , m_solsys->PVector(getTargetBodyId(), m_spacecraft.MJD())));

    m_best_de = m_best_de* std::max(0.1, de);

    dh = fabs(1.0-Dot(m_spacecraft.CState::RVector()
      , m_solsys->RVector(getTargetBodyId(), m_spacecraft.MJD())));
    m_best_dh = m_sfh*dh;
    break;
  }

  // propellant mass fitness calculation
  if (getCurrentSC().isSolarSail() == false)
  {
    m_best_mp = calcPropMassFitness(
      getCurrentSC().getInitialPropellantMass(),
      getCurrentSC().getPropellantMass());
  } else m_best_mp=0.0;
} // void CSimulation::calcBestValuesRendezvous(void)


// ----------------------------------------------------------------------------
void CSimulation::calcBestValuesFlyby(void)
{
  m_best_dr = m_best_proximity.GetNormedDistanceProximity(m_max_dr_AU);
  m_best_dv = m_best_proximity.GetNormedVelocityProximity(m_max_dv_AUd);

  // propellant mass fitness calculation
  if (getCurrentSC().isSolarSail() == false)
  {
    m_best_mp = calcPropMassFitness(
      getCurrentSC().getInitialPropellantMass(),
      getCurrentSC().getPropellantMass());
  } else m_best_mp=0.0;
} // void calcBestValuesFlyby(void)


// ----------------------------------------------------------------------------
void CSimulation::calcBestValuesACT(
  BodyID        target_body
 ,CSC&          sc
 ,CProximity&   proximity
 ,CSolarSystem* sol_sys
 ,double&       best_dr
 ,double        max_dr_AU
 ,double&       best_act
 ,double&       best_mp
)
{
  best_dr=proximity.GetNormedDistanceProximity(max_dr_AU);
  best_act = AUD2KMS*AUD2KMS*sc.getMass()*
    proximity.GetActVelocity(sc, sol_sys->State(target_body, sc.MJD()));

  // propellant mass fitness calculation
  if (sc.isSolarSail() == true) best_mp=0.0;
  else best_mp=calcPropMassFitness(sc.getInitialPropellantMass(),
    sc.getPropellantMass());

} // void CSimulation::calcBestValuesACT(...)

// ----------------------------------------------------------------------------
void CSimulation::calcBestValuesSPI(void)
{
} // void CSimulation::calcBestValuesSPI(...)




// -----------------------------------------------------------------------------
void CSimulation::calcBestValuesRetrogradeImpact(
   BodyID          target_body
  ,CSC&            sc
  ,CProximity&     proximity
  ,CSolarSystem*   sol_sys
  ,double&         best_dr
  ,double          max_dr_AU
  ,double&         best_dv
  ,double          max_dv_AUd
  ,double&         best_impv
  ,double&         best_rimph
  ,double&         best_rimpe
  ,double&         best_mp
  )
{
  m_best_sc = sc;

  best_dr=proximity.GetNormedDistanceProximity(
    max_dr_AU * m_accuracy_factor);
  best_dv=proximity.GetNormedVelocityProximity(
    max_dv_AUd * m_accuracy_factor);

  Vector hvec1 = sc.CState::RVector();
  Vector hvec2 = sol_sys->RVector(target_body, sc.MJD());

  // for retrograde impact
  double rimphdot = Dot(hvec1, -hvec2);
  double rimpdhtmp = 0.5*(1.0-rimphdot);
  double rimpdh   = 0.5*(pow(rimpdhtmp,5)+rimpdhtmp);
  best_rimph = 1000.0*rimpdh;

  if (getOptGoalType() != RETROGRADE_IMPACT_CRANKING)
  {
    best_impv = AUD2KMS*proximity.GetImpactVelocity(sc, sol_sys->State(target_body, sc.MJD()));

    Vector Pvec1 = sc.CState::PVector();
    Vector Pvec2 = sol_sys->PVector(target_body, sc.MJD());

    double rimpedot = Dot(Pvec1, Pvec2);
    double rimpdetmp = 0.5*(1.0-rimpedot);
    double rimpde   = 0.5*(pow(rimpdetmp,5)+rimpdetmp);
    best_rimpe = 1000.0*rimpde;

    // propellant mass fitness calculation
    if (sc.isSolarSail() == true) best_mp=0.0;
    else best_mp=calcPropMassFitness(sc.getInitialPropellantMass(),
      sc.getPropellantMass());
  }
}
//

// -----------------------------------------------------------------------------
void CSimulation::calcBestValues(
   BodyID          target_body
  ,CSC&            sc
  ,CProximity&     proximity
  ,CSolarSystem*   sol_sys
  ,double&         best_dr
  ,double          max_dr_AU
  ,double&         best_dv
  ,double          max_dv_AUd
  ,double&         best_da
  ,double&         best_de
  ,double&         best_di
  ,double&         best_dO
  ,double&         best_do
  ,double&         best_dM
  ,double&         best_dh
  ,double&         best_impv
  ,double&         best_rimph
  ,double&         best_rimpe
  ,double&         best_act
  ,double&         best_spi_a
  ,double&         best_spi_i
  ,double&         best_spi_e
  ,double&         best_mp
  ,double&         ga_credits
  ,bool            ga_perfomed // gravity assist performed ?
  ,double          sf_a        // sensitivity factor of semi major axis
  ,double          sf_e        // sensitivity factor of eccentricity
  ,double          sf_h        // sensitivity factor of angular momentum
)
{
  switch (tt)
  {

  case ACTPROBLEM:
    calcBestValuesACT(target_body, sc, proximity, sol_sys, best_dr, max_dr_AU
      , best_act, best_mp);
    break;

  case RIMPACT:
    calcBestValuesRetrogradeImpact(target_body, sc, proximity, sol_sys, best_dr
      , max_dr_AU, best_dv, max_dv_AUd, best_impv, best_rimph, best_rimpe
      , best_mp);    
    break;

  default:
    switch (getTargetStateType())
    {
    case ORBITAL_ELEMENTS_STATE:
      best_dr = best_dv = 1.0;
      break;
    default:
      best_dr=proximity.GetNormedDistanceProximity(
        max_dr_AU * m_accuracy_factor);
      best_dv=proximity.GetNormedVelocityProximity(
        max_dv_AUd * m_accuracy_factor);
      break;
    }
    break;
  } // switch (tt)

  if ( (tt != ACTPROBLEM) && (tt != RIMPACT) )
  {
    best_impv = AUD2KMS*proximity.GetImpactVelocity(sc, sol_sys->State(target_body, sc.MJD()));

    Vector evec1=sc.CState::PVector();
    Vector evec2=sol_sys->PVector(target_body, sc.MJD());

    double de=fabs(1.0-Dot(evec1, evec2));

    switch (getTargetStateType())
    {
    case BODY_STATE:
    case BODY_ORBIT_STATE:
    case FREE_ORBIT_STATE:
      best_da = sf_a *
        fabs(sc.CState::a()-sol_sys->a(target_body, sc.MJD()));
      best_de = sf_e *
        fabs(sc.CState::e()-sol_sys->e(target_body, sc.MJD()));
      best_de=best_de* std::max(0.1, de);
      break;
    case ORBITAL_ELEMENTS_STATE:
      break;
    default:
      
      break;
    }

    Vector hvec1 = sc.CState::RVector();
    Vector hvec2 = sol_sys->RVector(target_body, sc.MJD());

    double dh = fabs(1.0-Dot(hvec1, hvec2));
    best_dh = sf_h*dh;

    // for retrograde impact
    double rimphdot = Dot(hvec1, -hvec2);
    double rimpdhtmp = 0.5*(1.0-rimphdot);
    double rimpdh   = 0.5*(pow(rimpdhtmp,5)+rimpdhtmp);
    best_rimph = 10000.0*rimpdh;

    Vector Pvec1 = sc.CState::PVector();
    Vector Pvec2 = sol_sys->PVector(target_body, sc.MJD());

    double rimpedot = Dot(Pvec1, Pvec2);
    double rimpdetmp = 0.5*(1.0-rimpedot);
    double rimpde   = 0.5*(pow(rimpdetmp,5)+rimpdetmp);
    best_rimpe = 10000.0*rimpde;

    // for solar polar imager
    best_spi_a = sc.CState::a() - 0.48;
    best_spi_i = 3.0 * (abs(abs(sc.CState::i()) - 75.0*DEG2RAD));
    best_spi_e = sc.CState::e();

    // propellant mass fitness calculation
    if (sc.isSolarSail() == true) best_mp=0.0;
    else best_mp=calcPropMassFitness(sc.getInitialPropellantMass(),
      sc.getPropellantMass());

    // gravity assist check (currently switched off)
    if (ga_perfomed == true) ga_credits=0.0;
    else ga_credits=0.0;
  }

}; // end CSimulation::CalcBestValues


// ----------------------------------------------------------------------------
// private
void CSimulation::calcBestValues()
{
  m_steps_best = getCurrentStepCount();
  m_time_best = getFlightTime();
  setBestProximity(getCurrProximity());
  m_best_sc = m_spacecraft;

  switch (tt)
  {
  case RENDEZVOUS:calcBestValuesRendezvous();break;
  case FLYBY:calcBestValuesFlyby();break;

  case RIMPACT:
    calcBestValuesRetrogradeImpact(getTargetBodyId(), m_best_sc, m_best_proximity,
      m_solsys, m_best_dr, m_max_dr_AU, m_best_dv, m_max_dv_AUd, m_best_impv,
      m_best_rimph, m_best_rimpe, m_best_mp);
    break;

  case ACTPROBLEM:
    calcBestValuesACT(getTargetBodyId(), getCurrentSC(), m_best_proximity, m_solsys
      , m_best_dr, m_max_dr_AU, m_best_act, m_best_mp);
    break;

  case SPIPROBLEM:calcBestValuesSPI();break;
  case GTOC2TRAJ:calcBestValuesGTOC2();break;
  case TRAJ_TEST1:calcBestValuesTest1();break;
  case TRAJ_TEST2:calcBestValuesTest2();break;
  case TRAJ_TEST3:calcBestValuesTest3();break;
  case TRAJ_TEST4:calcBestValuesTest4();break;
  case TRAJ_TEST5:calcBestValuesTest5();break;
  case TRAJ_TEST6:calcBestValuesTest6();break;
  case TRAJ_TEST7:calcBestValuesTest7();break;
  case TRAJ_TEST8:calcBestValuesTest8();break;
  case TRAJ_TEST9:calcBestValuesTest9();break;
  }

  //calcBestValues(
  //  getTargetBodyId()
  // ,m_spacecraft
  // ,getCurrProximity()
  // ,m_solsys
  // ,m_best_dr
  // ,m_max_dr_AU
  // ,m_best_dv
  // ,m_max_dv_AUd
  // ,m_best_da
  // ,m_best_de
  // ,m_best_di
  // ,m_best_dO
  // ,m_best_do
  // ,m_best_dM
  // ,m_best_dh
  // ,m_best_impv
  // ,m_best_rimph
  // ,m_best_rimpe
  // ,m_best_act
  // ,m_best_spi_a
  // ,m_best_spi_i
  // ,m_best_spi_e
  // ,m_best_mp
  // ,m_ga_credits
  // ,m_ga_performed
  // ,m_sfa
  // ,m_sfe
  // ,m_sfh
  //);
}; // end CSimulation::CalcBestValues

// ----------------------------------------------------------------------------
inline CEAFitness* setFitnessValuesTest(CEAFitness* fit, double t)
{
  CInTranceFitness* f = 0;
  f = dynamic_cast<CInTranceFitness*>(fit);

  if (f == 0) FatalError("CSimulation::setFitnessValuesTest(CEAFitness*, \
 double) => dynamic cast failed!");
  return f;
} // inline setFitnessValuesTest(CEAFitness* fit, double t)

// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValuesTest0(CEAFitness* fit, double t)
{
  return setFitnessValuesTest(fit, t);
}

// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValuesTest1(CEAFitness* fit, double t)
{
  return setFitnessValuesTest(fit, t);
}

// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValuesTest2(CEAFitness* fit, double t)
{
  return setFitnessValuesTest(fit, t);
}

// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValuesTest3(CEAFitness* fit, double t)
{
  return setFitnessValuesTest(fit, t);
}

// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValuesTest4(CEAFitness* fit, double t)
{
  return setFitnessValuesTest(fit, t);
}

// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValuesTest5(CEAFitness* fit, double t)
{
  return setFitnessValuesTest(fit, t);
}

// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValuesTest6(CEAFitness* fit, double t)
{
  return setFitnessValuesTest(fit, t);
}

// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValuesTest7(CEAFitness* fit, double t)
{
  return setFitnessValuesTest(fit, t);
}

// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValuesTest8(CEAFitness* fit, double t)
{
  return setFitnessValuesTest(fit, t);
}

// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValuesTest9(CEAFitness* fit, double t)
{
  return setFitnessValuesTest(fit, t);
}

// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValuesACT(CInTranceFitness* fit)
{
  fit->setFitness(TIME_FITNESS,100.0*(1.0-m_time_best/getIntegrationInterval()));
  fit->setFitness(POS_FITNESS,log10(1.0/(m_best_dr<EPS?EPS:m_best_dr)));
  fit->setFitness(MPROP_FITNESS,m_best_mp);
  fit->setFitness(ACT_FITNESS, m_best_act);
  fit->setFitness(MIXED_FITNESS, -m_best_prox_sum);
  return fit;
} // CEAFitness* CSimulation::setFitnessValuesACT(CEAFitness* fit, double t)


// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValuesGTOC2(CInTranceFitness* fit)
{
  fit->setFitness(TIME_FITNESS,100.0*(1.0-m_time_best/getIntegrationInterval()));
  fit->setFitness(POS_FITNESS,log10(1.0/(m_best_dr<EPS?EPS:m_best_dr)));
  fit->setFitness(VEL_FITNESS,log10(1.0/(m_best_dv<EPS?EPS:m_best_dv)));
  fit->setFitness(MPROP_FITNESS,m_best_mp);

  if ( (m_best_dr<=1.0) && (m_best_dv<=1.0) )
  {
    fit->setFitness(MIXED_FITNESS, m_best_mp/ std::max(1.0, m_time_best));
  } // if ( (m_best_dr<=1.0) && (m_best_dv<=1.0) )
  else
  {
    fit->setFitness(MIXED_FITNESS, m_best_prox_sum);
  }
  return fit;
} // CEAFitness* setFitnessValuesGTOC2(CInTranceFitness* fit)


// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValuesOrbElements(CInTranceFitness* fit)
{
  const std::string method_name = 
    "CEAFitness* CSimulation::setFitnessValuesOrbElements(CIntranceFitness*\
    fit)";
  const std::string err_msg1 = method_name + 
    std::string(": invalid optimization goal!\n");

  double tmp = 0.0;
  bool match = true;

  // time fitness
  fit->setFitness(TIME_FITNESS, (getIntegrationInterval()-m_time_best));

  //semimajor axis
  if (m_data.useTargetOrbitalElement(SEMIMAJORAXIS)==true)
  {
    m_best_da = getTOEBestElementDeviation()[SEMIMAJORAXIS];
    fit->setFitness(A_FITNESS,log10(1.0/(m_best_da<EPS?EPS:m_best_da)));
    tmp += std::min(1.0,fit->getFitness(A_FITNESS));
    if (m_best_da>1.0) match = false;
  }

  // eccentricity
  if (m_data.useTargetOrbitalElement(ECCENTRICITY)==true)
  {
    m_best_de = getTOEBestElementDeviation()[ECCENTRICITY];
    fit->setFitness(E_FITNESS,log10(1.0/(m_best_de<EPS?EPS:m_best_de)));
    tmp += std::min(1.0, fit->getFitness(E_FITNESS));
    if (m_best_de>1.0) match = false;
  }

  // inclination
  if (m_data.useTargetOrbitalElement(INCLINATION)==true)
  {
    m_best_di = getTOEBestElementDeviation()[INCLINATION];
    fit->setFitness(I_FITNESS,log10(1.0/(m_best_di<EPS?EPS:m_best_di)));
    tmp += std::min(1.0, fit->getFitness(I_FITNESS));
    if (m_best_di>1.0) match = false;
  }

  // node
  if (m_data.useTargetOrbitalElement(NODE)==true)
  {
    m_best_dO = getTOEBestElementDeviation()[NODE];
    fit->setFitness(NODE_FITNESS,log10(1.0/(m_best_dO<EPS?EPS:m_best_dO)));
    tmp += std::min(1.0, fit->getFitness(NODE_FITNESS));
    if (m_best_dO>1.0) match = false;
  }

  // periapsis
  if (m_data.useTargetOrbitalElement(PERIAPSIS)==true)
  {
    m_best_do = getTOEBestElementDeviation()[PERIAPSIS];
    fit->setFitness(PERI_FITNESS,log10(1.0/(m_best_do<EPS?EPS:m_best_do)));
    tmp += std::min(1.0, fit->getFitness(PERIAPSIS));
    if (m_best_do>1.0) match = false;
  }

  fit->setFitness(MPROP_FITNESS, m_best_mp);

  if (/*m_best_element_change_rest_times.sum()*/match == true)
  {
    switch (getOptGoalType())
    {
    case TRANSFER_TIME:
      fit->setFitness(MIXED_FITNESS, 
          log10(getBestRemainingFlightTime()) 
        - log10(getMaxFlightTime()/ std::max(1,getMaxStepsCount()))
        + tmp/100.0);
      break;
    case PROPELLANT_MASS:
      fit->setFitness(MIXED_FITNESS, m_best_mp);
      break;
    default:
      FatalError(err_msg1);
      break;
    }
  } 
  else fit->setFitness(MIXED_FITNESS, -m_best_prox_sum);

  return fit;
} // CEAFitness* CSimulation::setFitnessValuesOrbElements(CIntranceFitness* fit)


// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValuesFlyBy(CEAFitness* f)
{
  f->setFitness(TIME_FITNESS, 100.0*(1.0-m_time_best/getIntegrationInterval()));
  f->setFitness(POS_FITNESS, log10(1.0/(m_best_dr<EPS?EPS:m_best_dr)));
  f->setFitness(VEL_FITNESS, log10(1.0/(m_best_dv<EPS?EPS:m_best_dv)));
  f->setFitness(MPROP_FITNESS, m_best_mp);
  return f;
} // void setFitnessValuesFlyBy(void)


// ----------------------------------------------------------------------------
CEAFitness* CSimulation::setFitnessValues(CEAFitness* fit, double time_slice)
{
  using namespace std;
  using namespace constraints;

  CInTranceFitness* f = 0;
  f = dynamic_cast<CInTranceFitness*>(fit);

  if (f == 0) FatalError("CSimulation::setFitnessValues(CEAFitness*, double)\
) => zero pointer fitness object!");

  if (time_slice != HUGE_VAL)
  {
    vector<CConstraint*>::iterator constraint = m_constraints.begin();
    while (constraint != m_constraints.end())
    {
      if ( !(*constraint)->isMatched() && (*constraint)->isHardConstraint() )
      {
        f->reset(MIN_FITNESS);
        return f;
      }
      ++constraint;
    }
    if (tt==ACTPROBLEM) return setFitnessValuesACT(f);

    else if ( (tt==RENDEZVOUS) && (getTargetStateType()==ORBITAL_ELEMENTS_STATE))
      return setFitnessValuesOrbElements(f);

    else if ( tt == GTOC2TRAJ ) setFitnessValuesGTOC2(f);

    else if ( tt==FLYBY ) return setFitnessValuesFlyBy(f);
    else 
    {
      f->setFitness(TIME_FITNESS, 100.0*(1.0-m_time_best/getIntegrationInterval()));
      f->setFitness(POS_FITNESS, log10(1.0/(m_best_dr<EPS?EPS:m_best_dr)));
      f->setFitness(VEL_FITNESS, log10(1.0/(m_best_dv<EPS?EPS:m_best_dv)));
      f->setFitness(A_FITNESS, log10(1.0/(m_best_da<EPS?EPS:m_best_da)));
      f->setFitness(E_FITNESS, log10(1.0/(m_best_de<EPS?EPS:m_best_de)));
      f->setFitness(I_FITNESS, log10(1.0/(m_best_di<EPS?EPS:m_best_di)));
      f->setFitness(NODE_FITNESS, log10(1.0/(m_best_dO<EPS?EPS:m_best_dO)));
      f->setFitness(PERI_FITNESS, log10(1.0/(m_best_do<EPS?EPS:m_best_do)));
      f->setFitness(M_FITNESS, log10(1.0/(m_best_dM<EPS?EPS:m_best_dM)));
      f->setFitness(H_FITNESS, log10(1.0/(m_best_dh<EPS?EPS:m_best_dh)));
      f->setFitness(IMPV_FITNESS, m_best_impv);
      f->setFitness(RIMPH_FITNESS, log10(1.0/(m_best_rimph<EPS?EPS:m_best_rimph)));
      f->setFitness(RIMPE_FITNESS, log10(1.0/(m_best_rimpe<EPS?EPS:m_best_rimpe)));
      f->setFitness(MPROP_FITNESS, m_best_mp);
      f->setFitness(GA_FITNESS, m_ga_credits);
      f->setFitness(DEFL_FITNESS, m_best_deflection);
      f->setFitness(SF_ACCURACY, m_accuracy_factor);
      f->setFitness(SPI_A_FITNESS, log10(1.0/(m_best_spi_a<EPS?EPS:m_best_spi_a)));
      f->setFitness(SPI_I_FITNESS, log10(1.0/(m_best_spi_i<EPS?EPS:m_best_spi_i)));
      f->setFitness(SPI_E_FITNESS, log10(1.0/(m_best_spi_e<EPS?EPS:m_best_spi_e)));
    }
  } // if (time_slice != HUGE_VAL)
  else f->reset(MIN_FITNESS);
  return f;
}; // end CSimulation::setFitnessValues


// ----------------------------------------------------------------------------
//private
void CSimulation::setLaunchState(
 const double time,
 const double velocity,
 const double azimuth,
 const double elevation)
{
  switch (getInitialStateType())
  {
  case BODY_STATE:
  case VECTOR_STATE:
    m_spacecraft.setLaunchState(
      time,
      getInitialBody(),
      velocity,
      azimuth,
      elevation);
    m_data.setInitialOrbitElement(SEMIMAJORAXIS, m_spacecraft.a(time));
    m_data.setInitialOrbitElement(ECCENTRICITY, m_spacecraft.e(time));
    m_data.setInitialOrbitElement(INCLINATION, m_spacecraft.i(time));
    m_data.setInitialOrbitElement(NODE, m_spacecraft.Omega(time));
    m_data.setInitialOrbitElement(PERIAPSIS, m_spacecraft.omega(time));
    m_data.setInitialOrbitElement(ANOMALY, m_spacecraft.M(time));
    break;
  case FREE_ORBIT_STATE:
    m_spacecraft.setCenterBodyId(getInitCentralBody());
    m_spacecraft.setCenterState(m_solsys->getBody(getInitCentralBody()));
    m_spacecraft.setLaunchState(
      time,
      m_data.getInitialOrbitElement(SEMIMAJORAXIS),
      m_data.getInitialOrbitElement(ECCENTRICITY),
      m_data.getInitialOrbitElement(INCLINATION),
      m_data.getInitialOrbitElement(NODE),
      m_data.getInitialOrbitElement(PERIAPSIS),
      m_data.getInitialOrbitElement(ANOMALY),
      velocity,
      azimuth, elevation);
    break;
  default:
    FatalError("void CSimulation::setLaunchState( const double time,\
const double velocity, const double azimuth, const double elevation):\
unknown initial state type!\n");
    break;
  }
  m_launch_sc = m_spacecraft;
}; // end CSimulation::setLaunchState


// ----------------------------------------------------------------------------
void CSimulation::updateProximityTargetBody(void)
{
  getCurrProximity().SetProximity(m_solsys->State(getTargetBodyId()
    , m_spacecraft.MJD()), m_spacecraft);
}

// ----------------------------------------------------------------------------
void CSimulation::updateProximityGtoc2(void)
{
  const int group = m_gtoc2_group[m_gtoc2_fp];
  CBody* asteroid = 0;
  int id = 0;

  const double t = m_t;
  const double delta_angle = 10.0*DEG2RAD;

  // heliocentric polar angles
  const double sc_phi = getCurrentSC().CState::phi();
  const double sc_theta = getCurrentSC().CState::theta();
  double as_phi=0, as_theta=0;
  double temp;
  m_temp_prox_sum = -HUGE_VAL;
  m_gtoc2_temp_asteroid = 0;
  while ( (asteroid = m_data.getGTOC2Body(group, id)) != 0)//scan all asteroids
  {
    as_phi = asteroid->phi(t);
    as_theta = asteroid->theta(t);

    if (fabs(as_phi-sc_phi) < delta_angle)
    { // phi criterion fulfilled

      if (fabs(as_theta-sc_theta) < delta_angle)
      { // theta criterion fulfilled
        getCurrProximity().SetProximity(*asteroid, m_spacecraft);
        const double dist_prox = getCurrProximity().GetNormedDistanceProximity(getMaxRelDistance());
        const double vel_prox = getCurrProximity().GetNormedVelocityProximity(getMaxRelVelocity());
        //temp = - getCurrProximity().GetNormedProximity(getMaxRelDistance()
        //  , getMaxRelVelocity());
        temp = -std::max(dist_prox, vel_prox);
        if (temp > m_temp_prox_sum)
        { // new temporary best proximity found
          m_temp_prox_sum = temp;
          m_gtoc2_temp_asteroid = asteroid;
          m_gtoc2_temp_best_bodies[m_gtoc2_fp] = asteroid->getId();
        } // if (temp > m_temp_prox)
      } // if (fabs(as_theta-sc_theta) < delta_angle)

    } // if (fabs(as_phi-sc_phi) < delta_angle)

    ++id;
  } // while ( (asteroid = m_data.getGTOC2Body(group, id)) != 0)
} // void CSimulation::updateProximityGtoc2(void)


// ----------------------------------------------------------------------------
void CSimulation::updateProximityTargetBodyOrbit(void)
{
#ifdef development_version
  assert(m_solsys!=0);
#endif // development_version

  const double t = m_spacecraft.MJD();
  const BodyID id = getTargetBodyId();
  const double a = m_solsys->a(id,t);
  const double e = m_solsys->e(id,t);
  const double i = m_solsys->i(id,t);
  const double Om= m_solsys->Omega(id,t);
  const double om= m_solsys->omega(id,t);
  const double nu= (getCurrentSC().CState::phi() - m_solsys->Omega(id, t)- m_solsys->omega(id,t));

  temp_state.SetState(a,e,i,Om,om,nu,KEPLERF, t, getInitCentralBody());
  getCurrProximity().SetProximity(temp_state, getCurrentSC());
} // void CSimulation::updateProximityTargetBodyOrbit(void)


// ----------------------------------------------------------------------------
void CSimulation::updateProximityTargetFreeOrbit(void)
{

  //getCurrProximity().SetProximity(CState(
  //   m_data.getTargetOrbitalElement(SEMIMAJORAXIS)
  //  ,m_data.getTargetOrbitalElement(ECCENTRICITY)
  //  ,m_data.getTargetOrbitalElement(INCLINATION)
  //  ,m_data.getTargetOrbitalElement(NODE)
  //  ,m_data.getTargetOrbitalElement(PERIAPSIS)
  //  ,m_spacecraft.CState::phi() - m_data.getTargetOrbitalElement(NODE)
  //  - m_data.getTargetOrbitalElement(PERIAPSIS)
  //  , KEPLERF)
  //  , m_spacecraft);
  temp_state.SetState(m_data.getTargetOrbitalElement(SEMIMAJORAXIS)
    ,m_data.getTargetOrbitalElement(ECCENTRICITY)
    ,m_data.getTargetOrbitalElement(INCLINATION)
    ,m_data.getTargetOrbitalElement(NODE)
    ,m_data.getTargetOrbitalElement(PERIAPSIS)
    ,m_spacecraft.CState::phi() - m_data.getTargetOrbitalElement(NODE)
    - m_data.getTargetOrbitalElement(PERIAPSIS)
  ,KEPLERF, getCurrentSC().MJD(), getInitCentralBody());
  getCurrProximity().SetProximity(temp_state, getCurrentSC());
} // void CSimulation::updateProximityTargetFreeOrbit(void)


// ----------------------------------------------------------------------------
void CSimulation::updateProximityTargetOrbitalElements(void)
{
  getCurrProximity().SetProximity(1.0, 1.0);
  getCurrProximity().setSemimajorAxisProximity(1.0e20);
  getCurrProximity().setEccentricityProximity(1.0e20);
  getCurrProximity().setInclinationProximity(1.0e20);
  if (m_data.useTargetOrbitalElement(SEMIMAJORAXIS)==true)
    getCurrProximity().setSemimajorAxisProximity(m_spacecraft.CState::a()
    , m_data.getTargetOrbitalElement(SEMIMAJORAXIS));
  if (m_data.useTargetOrbitalElement(ECCENTRICITY)==true)
    getCurrProximity().setEccentricityProximity(m_spacecraft.CState::e()
    , m_data.getTargetOrbitalElement(ECCENTRICITY));
  if (m_data.useTargetOrbitalElement(INCLINATION)==true)
    getCurrProximity().setInclinationProximity(m_spacecraft.CState::i()
    , m_data.getTargetOrbitalElement(INCLINATION));
  if (m_data.useTargetOrbitalElement(NODE)==true)
    getCurrProximity().setNodeProximity(m_spacecraft.CState::Omega()
    , m_data.getTargetOrbitalElement(NODE));
  if (m_data.useTargetOrbitalElement(PERIAPSIS)==true)
    getCurrProximity().setPericenterProximity(m_spacecraft.CState::omega()
    , m_data.getTargetOrbitalElement(PERIAPSIS));
  if (m_data.useTargetOrbitalElement(ANOMALY)==true)
    getCurrProximity().setMeanAnomalyProximity(m_spacecraft.CState::M()
    , m_data.getTargetOrbitalElement(ANOMALY));
} // void CSimulation::updateProximityTargetFreeOrbit(void)

// ----------------------------------------------------------------------------
// private
void CSimulation::updateProximity(void)
{
  switch (getTargetStateType())
  {
  case BODY_STATE:
    if (getTrajectoryType()==GTOC2TRAJ)
    {
      updateProximityGtoc2();
    }
    else updateProximityTargetBody();
    break;
  case BODY_ORBIT_STATE:updateProximityTargetBodyOrbit();break;
  case FREE_ORBIT_STATE:updateProximityTargetFreeOrbit();break;
  case ORBITAL_ELEMENTS_STATE:updateProximityTargetOrbitalElements();break;
  default:break;
  }
}; // end CSimulation::UpdateProximity


// ----------------------------------------------------------------------------
void CSimulation::statusMessage(const std::ostringstream& mess)
{
  Output(mess.str());      // output through global output object
  writeMessageToReportFile(mess); // and report file
} // void CSimulation::statusMessage(const std::ostringstream& mess)


// ----------------------------------------------------------------------------
// private
void CSimulation::addToSpacecraftTrajectory(const CSC&sc)
{
  m_spacecraft_trajectory[getCurrentStepCount()] = sc;
}; // end CSimulation::AddToSpacecraftTrajectory


// ----------------------------------------------------------------------------
//private
void CSimulation::PrintACTCompetitionLaunchWindows(void)
{
  double time(getMaxStartTime()+getIntegrationInterval()),t_p(time);
  std::vector<double> /*pericenter_times,*/pericenter_angles;
  std::vector< std::vector<double> > start_conjunction,start_opposition;

  double eps(1.0e-5),t_tmp(0.0);                                              // time accuracy, temporary variable
  std::ostringstream ss("");
  while ((t_p-m_solsys->OrbitPeriod(getTargetBodyId(),time) > getMinStartTime()))
  {
    m_ACT_pericenter_times.push_back(t_p = calcPerihelionTime(getTargetBodyId(), t_p, m_solsys, eps));
    pericenter_angles.push_back(m_solsys->phi(getTargetBodyId(), t_p));
    t_p-= m_solsys->OrbitPeriod(getTargetBodyId(),time)*0.98;
  }
  double t_start(getMinStartTime());
  for (unsigned int i=0;i<m_ACT_pericenter_times.size();++i)
  {
    t_start=getMinStartTime();
    start_conjunction.push_back(std::vector<double>());
    while (t_start<m_ACT_pericenter_times[i] && t_start<=getMaxStartTime())
    {
      start_conjunction[i].push_back(t_tmp=calcAngleTime(getInitialBody(), pericenter_angles[i],t_start,m_solsys,eps));
      t_start+=m_solsys->OrbitPeriod(getInitialBody(), t_start)*0.9999;
    }
  }
  for (unsigned int i=0;i<m_ACT_pericenter_times.size();++i)
  {
    t_start = getMinStartTime();
    start_opposition.push_back(std::vector<double>());
    while (t_start<m_ACT_pericenter_times[i] && t_start<=getMaxStartTime())
    {
      start_opposition[i].push_back(t_tmp=calcAngleTime(getInitialBody(), pericenter_angles[i]+PI,t_start,m_solsys,eps));
      t_start+=m_solsys->OrbitPeriod(getInitialBody(), t_start)*0.9999;
    }
  }

  // output
  double t_transfer(0.0);
  for (unsigned int i=0;i<m_ACT_pericenter_times.size();++i)
  {
    ss<<"\narrival at "<<m_ACT_pericenter_times[i]<<"MJD ("<<fixed << setprecision(3)
      <<pericenter_angles[i]*RAD2DEG<<"deg):\n";
    Output(ss.str());ss.seekp(0,ios_base::beg);ss.str("");
    ss<<setw(28)<<right<<"conjunction launch"<<"|"<<setw(28) << "opposition launch\n";
    ss<<"---------------------------------------------------------\n";
    for (unsigned int j=0;j<start_conjunction[i].size(),j<start_opposition[i].size();++j)
    {
      t_transfer=(m_ACT_pericenter_times[i]-(start_conjunction[i])[j]);
      if ( (t_transfer>0.0) && (t_transfer<getIntegrationInterval()))
      {
        ss << setw(10) << fixed << setprecision(3) << (start_conjunction[i])[j] << "MJD ("
          << setw(10) << fixed << setprecision(3) <<(m_ACT_pericenter_times[i]- (start_conjunction[i])[j])
          <<"d) | ";
      }
      t_transfer=(m_ACT_pericenter_times[i]-(start_opposition[i])[j]);
      if ( (t_transfer>0.0) && (t_transfer<getIntegrationInterval()))
      {
        if (ss.str().length()==0) ss << setw(28) << " " << "| ";
        ss << setw(10) << fixed << setprecision(3) << (start_opposition[i])[j] << "MJD ("
          << setw(10) << fixed << setprecision(3) <<(m_ACT_pericenter_times[i]- (start_opposition[i])[j])
          <<"d)";
      }
      if (ss.str().length()>0)
      {
        ss << "\n";Output(ss.str());ss.seekp(0,ios_base::beg);ss.str("");
      }
    }
  }
  Output("\n");
}

// -----------------------------------------------------------------------------
void CSimulation::setBodyInitialState(void)
{
} // void setBodyInitialState(void)

// -----------------------------------------------------------------------------
void CSimulation::setOrbitInitialState(void) throw (std::logic_error)
{
  const std::string method_name = 
    "void CSimulation::setOrbitInitialState(void)";
  const std::string err_msg_01 = 
    " couldn't insert initial state into solar system object!\n";

  const std::string body_name = "initial_orbit";

  const double body_mass = 1.0;

  const BodyID center_body = getInitCentralBody();

  const double mean_motion = CalcMeanMotion(
    m_data.getInitialOrbitElement(SEMIMAJORAXIS)
    , m_solsys->GM_AUd(center_body));

  CBody* initial_body = 0;    // an artificial start body being the launch state

  initial_body = new CBody(body_name, body_mass, getMinStartTime(),
    m_data.getInitialOrbitElement(SEMIMAJORAXIS),
    m_data.getInitialOrbitElement(ECCENTRICITY),
    m_data.getInitialOrbitElement(INCLINATION),
    m_data.getInitialOrbitElement(NODE),
    m_data.getInitialOrbitElement(PERIAPSIS),
    m_data.getInitialOrbitElement(ANOMALY),
    0.0,0.0,0.0,0.0,0.0,                                  // the rates of change
    mean_motion
    , KEPLERM
    , m_solsys->getBody(center_body));

  BodyID temp_id = 0;
  temp_id = m_solsys->insertBody(initial_body, true, 100);

  if (temp_id != 0)
  {
    setInitialBody(temp_id);
  }
  else throw std::logic_error(method_name+err_msg_01);

} // void setOrbitInitialState(void)

// -----------------------------------------------------------------------------
void CSimulation::setVectorInitialState(void)
{
  const std::string method_name = 
    "void CSimulation::setOrbitInitialState(void)";
  const std::string err_msg_01 = 
    " couldn't insert initial state into solar system object!\n";

  const std::string body_name = "initial_orbit";

  const double body_mass = 1.0;

  const BodyID center_body_id = getInitCentralBody();
  CBody* center_body = m_solsys->getBody(center_body_id);

  CBody* initial_body = 0;   // an artificial body object being the launch state
    
  initial_body = new CBody(body_name, body_mass
    , getMinStartTime()
    , m_data.getInitPosVec()(0)
    , m_data.getInitPosVec()(1)
    , m_data.getInitPosVec()(2)
    , m_data.getInitVelVec()(0)
    , m_data.getInitVelVec()(1)
    , m_data.getInitVelVec()(2)
    ,0,0,0,0,0,0
    ,CARTESIAN
    ,center_body);
  initial_body->setCenterBodyId(center_body_id);

  BodyID temp_id = 0;
  temp_id = m_solsys->insertBody(initial_body, true, 100);

  if (temp_id != 0)
  {
    setInitialBody(temp_id);
  }
  else throw std::logic_error(method_name+err_msg_01);
} // void setVectorInitialState(void)

// ----------------------------------------------------------------------------
void CSimulation::setInitialState(void) throw (std::invalid_argument)
{
  const std::string method_name = 
    "void CSimulation::setInitialState(void) throw (std::invalid_argument)";

  const std::string err_msg_01 = " unknown initial state type!\n";

  switch (getInitialStateType())
  {
  case BODY_STATE:
    setBodyInitialState();
    break;
  case FREE_ORBIT_STATE:
    setOrbitInitialState();
    break;
  case VECTOR_STATE:
    setVectorInitialState();
    break;
  default:
    throw std::invalid_argument(method_name + err_msg_01);
    break;
  }
} // void setInitialState(void)

// ----------------------------------------------------------------------------
void CSimulation::setTargetStateBody(void)
{
  if (getTrajectoryType() ==  GTOC2TRAJ)
  {
    setTargetBodyId(m_data.getGTOC2Body(3, 45)->getId());
  }
} // void setTargetStateElements(void)

// ----------------------------------------------------------------------------
void CSimulation::setTargetStateBodyOrbit(void)
{
} // void setTargetStateElements(void)

// ----------------------------------------------------------------------------
void CSimulation::setTargetStateOrbit(void)
{
} // void setTargetStateElements(void)

// ----------------------------------------------------------------------------
void CSimulation::setTargetStateElements(void)
{
  const std::string method_name = 
    "void CSimulation::setTargetStateElements(void)";
  const std::string err_msg_01 = 
    " couldn't insert target state into solar system object!\n";

  CBody* target_body = 0;    // used for orbital parameters target state
  m_best_initial_proximity = HUGE_VAL;
  const std::string body_name = "target_orbit";
  const double body_mass = 1.0;
  const BodyID center_body_id = getInitCentralBody();
  CBody* center_body = m_solsys->getBody(center_body_id);
  const double mean_motion = 
    CalcMeanMotion(m_data.getTargetOrbitalElement(SEMIMAJORAXIS)
    , m_solsys->GM_AUd(center_body_id));

  target_body = new CBody(
    body_name
    ,body_mass
    , MJD_J2000,
    m_data.getTargetOrbitalElement(SEMIMAJORAXIS),
    m_data.getTargetOrbitalElement(ECCENTRICITY),
    m_data.getTargetOrbitalElement(INCLINATION),
    m_data.getTargetOrbitalElement(NODE),
    m_data.getTargetOrbitalElement(PERIAPSIS),
    m_data.getTargetOrbitalElement(ANOMALY),
    0.0,0.0,0.0,0.0,0.0
    ,mean_motion
    ,KEPLERM
    ,center_body);

  target_body->setCenterBodyId(center_body_id);

  setTargetBodyId(m_solsys->insertBody(target_body, true, 100));

  if (getTargetBodyId()==0) std::logic_error(method_name+err_msg_01);

} // void setTargetStateElements(void)

// ----------------------------------------------------------------------------
void CSimulation::setTargetStateReserved01(void)
{
} // void setTargetStateReserved01(void)

// ----------------------------------------------------------------------------
void CSimulation::setTargetState(void) throw (std::invalid_argument)
{
  const std::string method_name = "void CSimulation::setTargetState(void)";
  const std::string err_msg_01 = ": unknown target state type!\n";
   // read target body or orbit
   
   switch (getTargetStateType())
   {
   case BODY_STATE:
     setTargetStateBody();
     break;
   case BODY_ORBIT_STATE:
     setTargetStateBodyOrbit();
     break;
   case FREE_ORBIT_STATE:
     setTargetStateOrbit();
     break;
   case ORBITAL_ELEMENTS_STATE:
     setTargetStateElements();
     break;
   case RESERVED_01:
     setTargetStateReserved01();
     break;
   default:
     throw std::invalid_argument(method_name + err_msg_01);
     break;
   }

} // void setTargetState(void) throw (std::invalid_argument)

// ----------------------------------------------------------------------------
void CSimulation::initSpacecraft(const CInputSC& scp_input)
throw (std::logic_error)
{
  const std::string method_name =
    "void CSimulation::initSpacecraft(const CInputSC&)";
  const std::string err_msg01 = 
    ": couldn't insert spacecraft into solar system!\n";

  m_spacecraft.setPropSys(scp_input.psystype);

  if (m_spacecraft.isSolarSail()) 
  {
    m_spacecraft.setLightness(scp_input.sailcraft_performance);
    m_spacecraft.setOpticalModel(
      scp_input.ideal_or_real_sail,
      0.0, // srd
      scp_input.max_temperature,
      scp_input.reflection_coefficient_new,
      scp_input.specular_reflection_coefficient_new,
      scp_input.emission_coefficient_front_new,
      scp_input.emission_coefficient_back_new,
      scp_input.reflection_coefficient_aged,
      scp_input.specular_reflection_coefficient_aged,
      scp_input.emission_coefficient_front_aged,
      scp_input.emission_coefficient_back_aged,
      scp_input.reflection_coefficient_halflife,
      scp_input.specular_reflection_coefficient_halflife,
      scp_input.emission_coefficient_front_halflife,
      scp_input.emission_coefficient_back_halflife,
      scp_input.non_lambertian_coefficient_front,
      scp_input.non_lambertian_coefficient_back);
  }
  if (m_spacecraft.isNEP())
  {
    m_spacecraft.setSCMass(scp_input.sc_mass);
    m_spacecraft.setInitialPropellantMassInterval(
      scp_input.min_mass_initial_propellant,
      scp_input.max_mass_initial_propellant);
    m_spacecraft.setMaxThrust(scp_input.thrust_max);
    m_spacecraft.setIsp(scp_input.isp);
    m_spacecraft.setThrottleType(scp_input.throttletype);
  } else if (m_spacecraft.isSEP())
  {
    m_spacecraft.setSCMass(scp_input.sc_mass);
    m_spacecraft.setInitialPropellantMassInterval(
      scp_input.min_mass_initial_propellant,
      scp_input.max_mass_initial_propellant);
    m_spacecraft.setCharacteristicPower(scp_input.char_power);
    m_spacecraft.setNumberOfThrusters(scp_input.thrusters);
    m_spacecraft.setThrottleType(scp_input.throttletype);
  }

  m_spacecraft.setCenterBodyId(getInitCentralBody());
  m_spacecraft.setCenterState(m_solsys->getBody(getInitCentralBody()));

  if (m_solsys->insertBody(&m_spacecraft, true, m_solsys->size()+10) == 0)
    throw std::logic_error(method_name + err_msg01);

} // void initSpacecraft(void)

// ----------------------------------------------------------------------------
// public
CSimulation::CSimulation(
const CInputSimulation& sim_input
,const CInputSC&         scp_input
,      CNC*              nc
,const string& rep_file_name) 
throw (std::invalid_argument, std::logic_error)
:m_data(sim_input)
,m_nc(nc)
,m_should_abort(0)
,m_rep(rep_file_name)
,m_sim_type(STANDARD)
,m_solsys(sim_input.m_solsys)
,m_spacecraft(scp_input, m_solsys)
,m_spacecraft_trajectory(0)
,m_curr_DSSC_max_angle(sim_input.getDSSCMaxStepAngle())
,m_steps(0)
,m_steps_best(0)
,m_time_best(0.0)
,m_time_max_dynamic(getIntegrationInterval())
,m_last_saved_individual()
,m_ga_body(NULL)
,m_initial_delements(6)
,m_initial_mp(0.0)
,m_initial_velocity(0.0), m_initial_azimuth(0.0), m_initial_elevation(0.0)
,m_max_dr_AU(sim_input.max_dr), m_max_dv_AUd(sim_input.max_dv)
,m_sfa(0.0)
,m_sfe(0.0)
,m_sfi(0.0)
,m_sfh(0.0)
,m_ga_performed(false)
,m_ga_credits(0.0)
,m_best_dr(0.0)
,m_best_dv(0.0)
,m_best_da(0.0)
,m_best_de(0.0)
,m_best_dh(0.0)
,m_best_impv(0.0)
,m_best_rimph(0.0)
,m_best_rimpe(0.0)
,m_best_act(0.0)
,m_best_spi_a(0.0)
,m_best_spi_i(0.0)
,m_best_spi_e(0.0)
,m_best_mp(0.0)
,m_best_deflection(0.0)
,m_best_s(6)
,m_temp_prox_sum(-HUGE)
,m_impact_time(sim_input.impact_time)
,m_impact_efficiency(sim_input.impact_efficiency)
,m_deflection(0.0)
,m_SpacecraftOrbit_RK4(NULL)
,m_SpacecraftOrbit_RKF54(NULL)
,m_SpacecraftOrbit_DOPRI8(NULL)
,m_SpacecraftOrbit_DESG(NULL)
,m_clock_angle_history(false)
,m_clock_angle_history_approximation(false)
,m_clock_angle_history_size(0)
,m_cone_angle_history(false)
,m_cone_angle_history_approximation(false)
,m_thrust_history(false)
,m_status(0)
,m_last_saved_start_time(0.0)
,m_last_saved_start_velocity(0.0)
,m_last_saved_start_azimuth(0.0)
,m_last_saved_start_elevation(0.0)
,m_last_saved_trajectory_approx()
,m_last_saved_trajectory_approx_steps()
,m_last_saved_clock_angle_approx_values()
,m_last_saved_clock_angle_approx_values_dot()
,m_last_saved_cone_angle_approx_values()
,m_last_saved_cone_angle_approx_values_dot()
,m_approx_put_index(0)
,m_approx_get_index(0)
,m_max_clock_jump(sim_input.max_clock_jump)
,m_approx_res_factor(sim_input.resolution_factor)
,m_last_saved_thrust_history(0)
,m_last_saved_clock_angle_history(0)
,m_last_saved_cone_angle_history(0)
{
  const std::string method_name = "CSimulation::CSimulation(...)";
  const std::string err_msg01 = " invalid ref. solar system object!\n";

  const unsigned int sz = getMaxStepsCount()+1;
  m_last_saved_trajectory .reserve(sz);
  m_last_saved_trajectory.resize(sz);
  m_last_saved_thrust_history.reserve(sz);
  m_last_saved_thrust_history.resize(sz);
  m_last_saved_clock_angle_history.reserve(sz);
  m_last_saved_clock_angle_history.resize(sz);
  m_last_saved_cone_angle_history.reserve(sz);
  m_last_saved_cone_angle_history.resize(sz);

  setAccuracyFactor(DFLT_ACCURACY_FACTOR);

  // set global parameters
  tt = getTrajectoryType();
  goal = getOptGoalType();

  if (m_spacecraft.CState::getRefSolSys() == 0)
    throw std::logic_error(method_name + err_msg01);

  // create spacecraft trajectory
  m_spacecraft_trajectory.resize(getMaxStepsCount()+1);

  // in case a new instance of CSimulation has been created, we need
  // new sensitivity factors (because of a static variable we do it here)
  init_sensitivity_factors = true;

  setInitialState();
  setTargetState();
  temp_state.setCenterBodyId(getInitCentralBody());
  temp_state.setCenterState(m_solsys->getBody(getInitCentralBody()));

  initSpacecraft(scp_input);
  
  initGravityAssistBodies();

  initConstraints();

  // set maximum final conditions at target body or orbit in AUd
  m_max_dr_km=m_max_dr_AU*AU2KM;
  m_max_dv_kms=m_max_dv_AUd*AUD2KMS;

  checkOptimizationGoal();

  // 1st ESA ACT Global Trajectory Optimization Competition
  if ((getOptGoalType()==ACT_DEFLECTION) && (getTrajectoryType()==ACTPROBLEM))
    PrintACTCompetitionLaunchWindows();

  // 2nd Global Trajectory Optimization Problem
  if ( getTrajectoryType() == GTOC2TRAJ )
  {
    m_gtoc2_fp = 0;
    m_gtoc2_group.resize(4,0);
    for (unsigned int i=0;i<4;++i) m_gtoc2_group[i] = (3-i);
    m_gtoc2_best_bodies.resize(4,0);
    m_gtoc2_temp_best_bodies.resize(4,0);
  } // if ( getTrajectoryType() == GTOC2TRAJ )


  initIntegrationVector();
  initEOM(&getCurrentSC());   // initialize problem-specific equations of motion
} // end CSimulation::CSimulation


// ----------------------------------------------------------------------------
void CSimulation::initConstraints(void)
{
  using namespace constraints;

  // set trajectory integration constraint of:

  // minimum distance from central body
  MIN_CENTRAL_BODY_CONSTRAINT.reset(new CMinCentralDistConstraint(this));
  m_constraints.push_back(MIN_CENTRAL_BODY_CONSTRAINT.get());

  // maximum solar sail temperature (in case of solar sail type s/c)
  if (getCurrentSC().isSolarSail() == true)
  {
    MAX_SAIL_TEMP_CONSTRAINT.reset(new CMaxSailTempConstraint(this));
    m_constraints.push_back(MAX_SAIL_TEMP_CONSTRAINT.get());
  }

  // proximity to target state
  if (tt != ACTPROBLEM && tt != GTOC2TRAJ)
  {
    TRGT_PROX_CONSTRAINT.reset(new CProximityConstraint(this));
    m_constraints.push_back(TRGT_PROX_CONSTRAINT.get());
  }

#ifdef CONSEP
  CONSEP_CONSTRAINT.reset(new CCONSEPConstraint(this));
  m_constraints.push_back(CONSEP_CONSTRAINT.get());
#endif // CONSEP

  // GTOC2 - specific problem constraint
  if (tt == GTOC2TRAJ)
  {
    GTOC2_CONSTRAINT.reset(new CGTOC2Constraint(this));
    m_constraints.push_back(GTOC2_CONSTRAINT.get());
  } // if (tt == GTOC2TRAJ)

} // void initConstraints(void)


// ----------------------------------------------------------------------------
void CSimulation::initGravityAssistBodies(void)
{
  if (m_data.number_of_ga_bodies != 0) 
  {
    m_ga_body=new bool[NEPTUNE+1];
    for (int i=MERCURY; i<=NEPTUNE; ++i)
    {
      m_ga_body[i] = m_data.ga_body[i];
    }
    getCurrentSC().GravityAssistSteeringOn();
  } else 
  { 
    getCurrentSC().GravityAssistSteeringOff(); 
  }
} // void initGravitiyAssistBodies(void)


// ----------------------------------------------------------------------------
void CSimulation::checkOptimizationGoal(void) throw (std::invalid_argument)
{
  if ( getCurrentSC().isSolarSail() && (getOptGoalType() == PROPELLANT_MASS) )
    throw (std::invalid_argument("CSimulation::checkOptimizationGoal(): \
propellant mass minimization not possible for solar sail"));
} // void checkOptimzationGoal(void)


// ----------------------------------------------------------------------------
void CSimulation::calcProposedAccuracy(void) 
{

  double a=m_solsys->a(getTargetBodyId(), m_spacecraft.MJD());

  double e=m_solsys->e(getTargetBodyId(), m_spacecraft.MJD());

  const double dt = (getIntegrationInterval()/getMaxStepsCount());

  cout<<endl<<"length of one integration step: " <<dt<<" days" <<endl;

  cout<<"at target body's perihel:" <<endl;

  double v=sqrt(GM_SUN_AUD/a*(1+e)/(1-e));
  cout<<"  velocity of target body: " <<v*AUD2KMS<<" km/s" <<endl;

  double dx=v*dt*AU2KM;
  cout<<"  distance of one integration step: " <<dx<<" km" <<endl;

  if (m_spacecraft.isSolarSail()) {
    const double lambda=m_spacecraft.getLightness()/pow(a*(1-e), 2.0);

    const double ac=lambda*GM_SUN_AUD;
    cout<<"  maximum acceleration: " <<ac*AU/86400.0/86400.0*1000.0<<
     " mm/s^2" <<endl;

    const double dv=ac*sin(0.61548)*cos(0.61548)*cos(0.61548)*dt;
    cout<<"  maximum \"breaking\" DeltaV during one integration step: " <<dv*
     AUD2KMS*1000.0<<" m/s" <<endl;

    const double dE=v*dv+0.5*dv*dv;
    cout<<"  maximum \"breaking\" DeltaE during one integration step: " <<dE*
     AUD2KMS*AUD2KMS<<" km^2/s^2" <<endl;
  }

  cout<<"at target body's aphel:" <<endl;

  v=sqrt(GM_SUN_AUD/a*(1-e)/(1+e));
  cout<<"  velocity of target body: " <<v*AUD2KMS<<" km/s" <<endl;

  dx=v*dt*AU2KM;
  cout<<"  distance of one integration step: " <<dx<<" km" <<endl;

  if (m_spacecraft.isSolarSail()) {
    const double lambda=m_spacecraft.getLightness()/pow(a*(1+e), 2.0);

    const double ac=lambda*GM_SUN_AUD;
    cout<<"  maximum acceleration: " <<ac*AU/86400.0/86400.0*1000.0<<
     " mm/s^2" <<endl;

    const double dv=ac*sin(0.61548)*cos(0.61548)*cos(0.61548)*dt;
    cout<<"  maximum \"breaking\" DeltaV during one integration step: " <<dv*
     AUD2KMS*1000.0<<" m/s" <<endl;

    const double dE=v*dv+0.5*dv*dv;
    cout<<"  maximum \"breaking\" DeltaE during one integration step: " <<dE*
     AUD2KMS*AUD2KMS<<" km^2/s^2" <<endl;
  }

  double E_target=0.5*v*v-GM_SUN_AUD/(a*(1+e));

  a=m_solsys->a(getInitialBody(), m_spacecraft.MJD());

  e=m_solsys->e(getInitialBody(), m_spacecraft.MJD());

  v=sqrt(GM_SUN_AUD/a*(1-e)/(1+e));
  double E_launch=0.5*v*v-GM_SUN_AUD/(a*(1+e));

  cout<<"specific energy of launch body: " <<E_launch*AUD2KMS*AUD2KMS<<
   " km^2/s^2" <<endl;
  cout<<"specific energy of target body: " <<E_target*AUD2KMS*AUD2KMS<<
   " km^2/s^2" <<endl;
  cout<<"specific energy difference: " <<fabs(E_target-E_launch)*AUD2KMS*
   AUD2KMS<<" km^2/s^2" <<endl;

}; // end CSimulation::CalcProposedAccuracy


// ----------------------------------------------------------------------------
void CSimulation::performGravityAssist(Vector&ys, Vector&x, Vector&v, double&t)
{
  bool print_ga_performed=false;

  for (int gabody=MERCURY; gabody<=NEPTUNE; gabody++) 
  {

    if (m_ga_body[gabody]) {
      if (Norm(m_solsys->Cart(gabody, t).Slice(0, 2)
              -m_spacecraft.CState::Cart().Slice(0, 2))<m_solsys->RSOI(gabody, t)) {
         // set neurocontroller input

         m_spacecraft.getSteering().setNCInput(
          m_nc, m_spacecraft,
          m_solsys->State(getTargetBodyId(), t),
          &m_spacecraft, &m_solsys->State(gabody, t)
         );

         Vector input_ncp(m_spacecraft.getCtrlInputDim());
         int k;
         for (k=0; k<m_spacecraft.getCtrlInputDim(); ++k)
           input_ncp(k)=m_spacecraft.getCtrlInput(k);

         // calculate neurocontroller output

         Vector nc_output(m_spacecraft.getCtrlOutputDim());
         m_nc->propagate(input_ncp, nc_output);

         for (k=0; k<m_spacecraft.getCtrlOutputDim(); ++k)
           m_spacecraft.setCtrlOutput(k, nc_output(k));

         double ii, jj;
         if (m_spacecraft.isSolarSail()) {
           ii= std::max(std::min(nc_output(3), 1.0), -1.0);
           jj= std::max(std::min(nc_output(4)/10.0, 1.0), -1.0);
         } else {
           ii= std::max(std::min(nc_output(4), 1.0), -1.0);
           jj= std::max(std::min(nc_output(5)/10.0, 1.0), -1.0);
         }
         double delta=sqrt(ii*ii+jj*jj);
         if (delta>1.0) {
           ii/=delta;
           jj/=delta;
         }
         if (delta>0.99) {
           ii*=0.99; //!!!qd!!!
           jj*=0.99; //!!!qd!!!
         }

         // calculate energy w.r.t. the gravity assist body

         CState delta_state(m_spacecraft.CState::Cart()-m_solsys->Cart(gabody, t));
         Vector delta_r =  m_spacecraft.CState::Cart().Slice(0,2)
                         - m_solsys->Cart(gabody, t).Slice(0,2);
         Vector delta_v =  m_spacecraft.CState::Cart().Slice(3,5)
                         - m_solsys->Cart(gabody, t).Slice(3,5);

         double norm_delta_r=Norm(delta_r);
         double norm_delta_v=Norm(delta_v);
         double energy=
          0.5*norm_delta_v*norm_delta_v-
          m_solsys->State(gabody, t).getGM(AUD)/norm_delta_r;

         if (energy>0.0) {
           cout<<" " <<ii<<" " <<jj<<" ";
           m_solsys->GravityAssistAnalytical(gabody, m_spacecraft, ii, jj);
           t=m_spacecraft.MJD();
           x=Vector(m_spacecraft.CState::r(), m_spacecraft.CState::phi(),
            m_spacecraft.CState::theta());
           v=Vector(m_spacecraft.CState::dr(), m_spacecraft.CState::dphi(),
            m_spacecraft.CState::dtheta());
           ys=m_spacecraft.PolarStatePropellantMassSRDVector();
           m_ga_performed=true;
           print_ga_performed=true;

           cout<<" -";
           switch (gabody) {
           case MERCURY:
           cout<<"M";
           break;
           case VENUS:
           cout<<"V";
           break;
           case EARTH:
           cout<<"E";
           break;
           case MARS:
           cout<<"M";
           break;
           case JUPITER:
           cout<<"J";
           break;
           case SATURN:
           cout<<"S";
           break;
           case URANUS:
           cout<<"U";
           break;
           case NEPTUNE:
           cout<<"N";
           break;
           default:
           break;
           } // end switch

         } // end if

      } // end if

      if (print_ga_performed) {
        cout<<"GA performed" <<endl;
        print_ga_performed=false;
      }

    } // end if
  } // end for

}; // end CSimulation::PerformGravityAssist


const std::string mn_001 = "void CSimulation::performSteering(void)";
const std::string err_msg_001_01 = mn_001 + ": unknown simulation type!\n";
const std::string err_msg_001_02 = mn_001 + ": unknown target state type!\n";
const std::string err_msg_001_03 = mn_001 + ": approximation not implemented!\n";
//-----------------------------------------------------------------------------
void CSimulation::performSteering(void)
{
  if (NativeBridgeProbe::enabled())
  {
    if (NativeBridgeProbe::mercury_profile())
    {
      NativeBridgeProbe::require(m_sim_type==STANDARD && tt==RENDEZVOUS
        && getTargetStateType()==BODY_STATE && getTargetBodyId()==MERCURY
        && isHeliocentric() && m_data.number_of_ga_bodies==0,
        "Mercury bridge requires standard heliocentric BODY_STATE rendezvous without gravity assists");
      NativeBridgeProbe::require(getIntegrationInterval()==600.0 && getMaxStepsCount()==600
        && getMinStartTime()==52654.5 && getMaxStartTime()==52654.5
        && getInitialStateType()==BODY_STATE && getInitialBody()==EARTH,
        "Mercury bridge requires original 600-day/600-step Earth departure at MJD 52654.5");
      NativeBridgeProbe::require(getMinLaunchVelocity()==0.0 && getMaxLaunchVelocity()==0.0
        && getMinLaunchAzimuthAngle()==0.0 && getMaxLaunchAzimuthAngle()==0.0
        && getMinLaunchElevationAngle()==0.0 && getMaxLaunchElevationAngle()==0.0,
        "Mercury bridge requires original fixed zero launch excess velocity and angles");
      NativeBridgeProbe::require(getIntegratorType()==_RKF54 && !useDAC()
        && NativeBridgeProbe::accepts_mercury_precision(NativeBridgeProbe::probe().profile,
             getMaxRelIntegrationError(),getMaxAbsIntegrationError())
        && getMinSolarDistance()==0.2 && m_spacecraft.isSolarSail() && m_spacecraft.isIdealSail(),
        "Mercury bridge requires profile-approved RKF54 tolerances, original solar-distance bound and ideal-sail model");
    }
    else
    {
      NativeBridgeProbe::require(m_sim_type == STANDARD && tt == ACTPROBLEM
        && isHeliocentric() && m_data.number_of_ga_bodies == 0,
        "ideal-sail diagnostic requires standard heliocentric ACT fixture without gravity assists");
      NativeBridgeProbe::require(getIntegrationInterval()>0.0 && getIntegrationInterval()<=10.0,
        "ideal-sail diagnostic requires a native horizon of at most ten days");
    }
    NativeBridgeProbe::require(!getStatusDSSC() && !getDynamicIntegrationIntervalReduction(),
      "ideal-sail diagnostic requires fixed native macro-step policy");
    NativeBridgeProbe::require(getMinStartTime()==getMaxStartTime()
      && getMinLaunchVelocity()==getMaxLaunchVelocity()
      && getMinLaunchAzimuthAngle()==getMaxLaunchAzimuthAngle()
      && getMinLaunchElevationAngle()==getMaxLaunchElevationAngle(),
      "ideal-sail diagnostic requires fixed departure epoch and launch conditions");
    if (NativeBridgeProbe::mercury_profile()) NativeBridgeProbe::note_mercury_profile_checked();
    if (NativeBridgeProbe::replaying())
    {
      // Preserve CSV topology metadata only. These are explicit placeholders,
      // not recomputed NN outputs; bridge physical exports contain no NN data.
      m_spacecraft.setCtrlInputDim(m_nc->getNumberOfInputNeurons());
      m_spacecraft.setCtrlOutputDim(m_nc->getNumberOfOutputNeurons());
      m_spacecraft.getCtrlInputVec().Set(0.0);
      m_spacecraft.getCtrlOutputVec().Set(0.0);
      NativeBridgeProbe::import_command(m_spacecraft);
      m_spacecraft.UpdateAccDeltaVec(0.0);
      NativeBridgeProbe::record_command(m_spacecraft);
      return;
    }
  }
  switch (m_sim_type)
  {

  case STANDARD:
    switch (getTargetStateType())
    {

    case ORBITAL_ELEMENTS_STATE:
      // update target state object
      m_solsys->State(getTargetBodyId(), m_t).SetState(
        (m_data.useTargetOrbitalElement(SEMIMAJORAXIS)?
        m_data.getTargetOrbitalElement(SEMIMAJORAXIS):m_spacecraft.CState::a()),
        (m_data.useTargetOrbitalElement(ECCENTRICITY)?
        m_data.getTargetOrbitalElement(ECCENTRICITY):m_spacecraft.CState::e()),
        (m_data.useTargetOrbitalElement(INCLINATION)?
        m_data.getTargetOrbitalElement(INCLINATION):m_spacecraft.CState::i()),
        (m_data.useTargetOrbitalElement(NODE)?
        m_data.getTargetOrbitalElement(NODE): m_spacecraft.CState::Omega()),
        (m_data.useTargetOrbitalElement(PERIAPSIS)?
        m_data.getTargetOrbitalElement(PERIAPSIS):m_spacecraft.CState::omega()),
        (m_data.useTargetOrbitalElement(ANOMALY)?
        m_data.getTargetOrbitalElement(ANOMALY): m_spacecraft.CState::M()),
        KEPLERM, m_t, m_spacecraft.getCenterBodyId());

      m_spacecraft.performSteering(m_nc
        , m_solsys->State(getTargetBodyId(), m_t));

      m_spacecraft.UpdateAccDeltaVec(0.0);
      break;
    case BODY_STATE:

    case BODY_ORBIT_STATE:

    case FREE_ORBIT_STATE:
      m_spacecraft.performSteering(m_nc
        , m_solsys->State(getTargetBodyId(), m_t));
      m_spacecraft.UpdateAccDeltaVec(0.0);
      break;
    default:
      FatalError(err_msg_001_02);
      break;
    } // switch (getTargetStateType())
    break;

  case APPROXIMATION:
    FatalError(err_msg_001_03);
    break;

  default:
    FatalError(err_msg_001_01);
    break;
  } // switch (m_sim_type)
  NativeBridgeProbe::record_command(m_spacecraft);
} // void CSimulation::performSteering(void)


//-----------------------------------------------------------------------------
void CSimulation::writeOutputFiles(const string& csv_file,
                                   const string& dat_file,
                                   const string& gesop_file,
                                   const string& wrl_file,
                                   const string& ang_file)
{
  Output(std::string("... writing simulation output files\n\n"));
  const unsigned int w(5);
  ostringstream oss("");

  bool csv_output(true);
  bool dat_output(true);
  bool gesop_output(true);
  bool vrml_output(true);
  bool angle_output(true);

  switch (m_sim_type)
  {
  case STANDARD: break;
  case APPROXIMATION:break;
  default:break;
  }

  // comma seperated value file
  if (csv_output == true)
  {
    if (!csv_file.empty())
    {
      oss.str(" ");oss<<setw(w)<<left<<"CSV"<<":";Output(oss.str());
      switch (m_sim_type)
      {
      case STANDARD:
      case APPROXIMATION:
      default:
        PrintCSV(csv_file.c_str());
        switch (tt)
        {
        case ACTPROBLEM:writeACTCsv(csv_file);break;
        case GTOC2TRAJ:writeGTOC2Csv(csv_file);break;
        default:break;
        }
        Output("ok\n");
        break;
      }
    } else Output(" empty file name string (csv - file)\n");
  }


  // DAT (TECPLOT) file
  if (dat_output == true)
  {
    if (!dat_file.empty())
    {
      oss.str(" ");oss.seekp(0);oss<<setw(w)<<left<<"DAT"<<":";
      Output(oss.str());
      PrintDAT(dat_file.c_str());Output("ok\n");
    }
    else Output(" empty file name string (dat - file)\n");
  }


  // GESOP file
  if (gesop_output == true)
  {
    if (!gesop_file.empty())
    {
      oss.str(" ");oss.seekp(0);oss<<setw(w)<<left<<"GESOP"<<":";
      Output(oss.str());
      PrintGESOP(gesop_file.c_str());Output("ok\n");
    }
    else Output(" empty file name string (GESOP - file)\n");
  }


  // VRML file
  if (vrml_output == true)
  {
    if (!wrl_file.empty())
    {
      oss.str(" ");oss.seekp(0);oss<<setw(w)<<left<<"VRML"<<":";
      Output(oss.str());
      switch (m_sim_type)
      {
      case STANDARD:
      case APPROXIMATION:
      default:
        PrintVRML(wrl_file.c_str());
        Output("ok\n");
        break;
      }
    }
    else Output(" empty file name string (wrl - file)\n");
  }

  // steering angles output
  if (angle_output == true)
  {
    if (!ang_file.empty())
    {
      oss.str(" ");oss.seekp(0);oss<<setw(w)<<left<<"ANGLS"<<":";
      Output(oss.str());
      PrintANG(ang_file.c_str());Output("ok\n");
    }
    else Output(" empty file name string (ang - file)\n");
  }

}; // end CSimulation::WriteOutputFiles



//-----------------------------------------------------------------------------
void CSimulation::evaluateInterpolated(void)
{
  int substep_factor = 1;

  int output_points = getMinOutputPoints();

  while (output_points > getMaxStepsCount()) {
    substep_factor*=2;
    output_points/=2;
  }

  // copy old trajectory to local variable
  SCTraj old_trajectory;
  old_trajectory.swap(m_spacecraft_trajectory);

  m_spacecraft_trajectory.resize((substep_factor * getMaxStepsCount()) + 1);

  initEOM(&old_trajectory[0]);

  m_steps = -1;

  Vector ys(8), x(3), v(3);

  int i,j;
  for (i=0; i <= m_steps_best; i++) {

    m_spacecraft=old_trajectory[i];

    double t=m_spacecraft.MJD();

    double h = (getIntegrationInterval()/getMaxStepsCount()/substep_factor);

    initDESSolver(ys, t);

    incCurrentStepCounter();

    addToSpacecraftTrajectory(m_spacecraft);

    if (i < m_steps_best) { // not last step

      for (j = 1; j < substep_factor; j++) {

        integrate(ys, t, h);
        m_spacecraft.UpdateAccDeltaVec(h);

        incCurrentStepCounter();
        addToSpacecraftTrajectory(m_spacecraft);

      } // end for j
    } // end if

  } // end for i

  m_steps_best *= substep_factor;

  freeEOM();

} // end CSimulation::EvaluateInterpolated




//-----------------------------------------------------------------------------
CEAFitness* CSimulation::evaluate(
                                   CIndividual& solution
                                  ,const bool final
                                  ,const bool monitor
                                  ,const double time)
{
  switch (getSimType())
  {
  case STANDARD:
    evaluateStandardSolution(solution, final, monitor, time);
    return 0;
    break;
  case APPROXIMATION:
    evaluateApproximationSolution(solution,final,monitor,time);
    return 0;
    break;
  default:return 0;//new CInTranceFitness;
  }
}

// some helper template functions

//-----------------------------------------------------------------------------
template<class T>
inline double abs_dev(const T val, const T ref_val)
{
  return double(val)-double(ref_val);
}

//-----------------------------------------------------------------------------
template<class T>
inline double rel_dev(const T val, const T ref_val)
{
  return abs_dev(val, ref_val)/double(ref_val);
}

//-----------------------------------------------------------------------------
template<class T>
inline double calcRatioSubFitness(const double offset, const T val
                                , const T ref_val)
{
  return offset-(double(val)/double(ref_val));
}

//-----------------------------------------------------------------------------
template<class T>
inline double calcRatioSubFitness2(const T val, const T min, const T max)
{
  if (val<=max && val>=min && max!=min)
  return (double(val)-double(min))/(double(max)-double(min));
  else return(0.0);
}


//-----------------------------------------------------------------------------
template<class T>
inline void bound(T& what, const double lo, const double up)
{
  for (unsigned long int i = 0; i<what.size(); ++i) 
    what[i]=bounded(what[i],lo,up);
}


//-----------------------------------------------------------------------------
/*!
 * \brief Counts sign changes in given data structure.
 */
template<class T>
unsigned int CountSignChanges(const T& d)
{
  if (d.size() > 2)
  {
    unsigned int scc(0);           // counter of sign changes in first derivative
    bool upwards(d[0]<d[1]?true:false);// true if value history is ascending
    T::value_type last(d[1]);
    for (unsigned int i = 2; i<(d.size()-1);++i)
    {
      if ((upwards==true && (d[i]-d[i-1])<0) || (upwards==false && (d[i]-d[i-1])>0))
      {
        ++scc;
        upwards = (upwards==true?false:true);
        last = d[i];
      }
    }
    return scc;
  }
  else return 1;
}


//-----------------------------------------------------------------------------
/*!
 *
 */
template<class T>
  unsigned int CountSignChangesDerivative(const T& d)
{
  if (d.size()>1)
  {
    unsigned int sc(0);
    bool pos(d[0]>0?true:false);
    for (unsigned int i=1;i<d.size()-1;++i)
    {
      if ((pos==true && d[i]<0) || (pos==false && d[i]>0))
      {
        ++sc;
        pos = (pos==true?false:true);
      }
    }
    return sc;
  }
  else return 0;
}


//-----------------------------------------------------------------------------
template<class R, class S, class T, class U, class V>
void PrepareArray(const R& x  , const S& y  , const T& s
                ,       U& n_x,       V& n_y)
{
  if (x.size()==y.size() && x.size()>0 && s.size()<x.size())
  {
    n_x.resize(s.size()+1, 0);n_y.resize(n_x.size(),0);
    n_x[0] = x[0]; n_y[0]=y[0];           // set first values
    for (uint i=0;i<s.size();++i)
    {
      if (s[i]>x.size()) break;   // check for correct indices
      n_x[i+1] = x[s[i]];
      n_y[i+1] = y[s[i]];
    }
  }
}



//-----------------------------------------------------------------------------
inline double MapNcOutput(
 const CNC* nc,
 const double min_val,
 const double max_val,
       double nc_o)
{
  double ret_val(min_val);
  switch (nc->getActivationFunction())
  {
    case SIGMOID:
      break;
    case TANH:
    case SIMPLE:
      nc_o = (nc_o+0.5)*2.0;
      break;
    case SIGMOIDIDENTITY:
    case SIMPLEIDENTITY:
    case TANHIDENTITY:
      return bounded(nc_o, min_val, max_val);
      break;
    default:
      Output("setSteeringMoments(....): unknown activation function!\n");
      break;
  }
  return nc_o*(max_val-min_val)+min_val;
}

//------------------------------------------------------------------------------
void applyApproxNC1(
       CNC&                                  nc,
       const std::vector<double>&            time,
       const std::vector<double>::size_type  i_i1,
 const std::vector<double>::size_type        i_i2,
 const std::vector<double>&                  data,
       std::vector<uint>&                    new_seg,
       std::vector<double>&                  new_data,
 const int                                   i_o1,
 const double                                min_val,
 const double                                max_val)
{
  const int s(std::min( time.size(), (i_i2-i_i1+1) ) );
  if ( (s>2) && (s<=int(data.size())) && (nc.getNumberOfOutputNeurons()==2) )
  {
    // init control variables
    const int n(nc.getNumberOfInputNeurons()), m(n/2);

    int i_l_in(0), i_r_in(m),i_l_nc(m),i_r_nc(std::min(m+s, n)-1);
    const double init_val(-1.0);
    int i(0), di(0), i_store(i_o1),h(0);
    Vector input(n), output(2);
    input.Set(init_val);
    new_seg.clear();
    new_data.clear();

    for (int p=0;p < (i_r_in-i_l_in); ++p)   // set input vector
      input[i_l_nc+p]=data[i_l_in+p];
    for (;;)
    {
      nc.propagate(input, output);               // apply neuronal net
      di = int((i_r_nc-m)*MapNcOutput(&nc, 0.0, 1.0, output[0]));
      di = std::max(10, di);
      i += di;
      i= std::min(i,s-1);

      if ( (i_store+10) > int(new_seg/*time*/.capacity()) )
      {
        int new_size( int(std::max(50, i_store)*1.10) );
        new_seg/*time*/.reserve(new_size);
        new_data.reserve(new_size);
      }
      new_seg.push_back(i);//time[i];
      new_data.push_back(MapNcOutput(&nc, min_val, max_val, output[1]));
      ++i_store;

      if ( i == (s-1) ) break;
      i_l_in = std::max(     0, (i-m) );
      i_r_in = std::min( (s-1), (i+m) );
      h = (m - (i_r_in-i));

      while (h>0)
      {
        input(n-h) = init_val;--h;
      }

      i_l_nc = std::max(i_l_nc-di, 0);
      i_r_nc = m + (i_r_in - i)-1;

      for (int p=0;p < (i_r_in-i_l_in); ++p)   // set input vector
      input[i_l_nc+p]=data[i_l_in+p];
    } //for (;;)

  } //   if (i_i1 < (i_i2-1))

} // void applyApproxNC1(....)


//------------------------------------------------------------------------------
void applyApproxNC2(
       CNC&                                  nc,
 const std::vector<double>&                  time,
 const std::vector<double>::size_type        i_i1,
 const std::vector<double>::size_type        i_i2,
 const std::vector<double>&                  data,
       std::vector<uint>&                    new_seg,
       std::vector<double>&                  new_data,
 const int                                   i_o1,
 const double                                min_val,
 const double                                max_val)
{
  const int min_seg_size(5);
        unsigned int input_idx(i_i1);
  const int last_input_idx(i_i2);
  const unsigned int n_input_values(std::min( time.size(), (i_i2-i_i1+1) ) );
  unsigned int input_values(0);
  const unsigned int max_input_values(nc.getNumberOfInputNeurons());
        unsigned int segment(0);
  const double init_val(0.0);            // nc input vector initialization value
  unsigned int last_clock_seg_end(0);
  int i(i_i1),di(0);
  int i_store(0);
  Vector input(max_input_values);                     // neural net input vector
  Vector output(nc.getNumberOfOutputNeurons());   // neural net output vector

  // determination of clock/cone angle knot points (for spline interpolation)
  // and modification values with neural net
  while (i < last_input_idx)
  {
      input_values = std::min(n_input_values-i,max_input_values);
      input.Set(init_val);                               // reset input vector

      for (uint k=0;k<input_values;++k) input[k] = data[i+k];

      nc.propagate(input, output);                  // applicate neural net

      // calculate next knot point
      di = 1+uint(MapNcOutput(&nc, 1.0,double(input_values), output[0]));
      i += di;

      if ( (i_store+10) > int(new_seg.capacity()) )
      {
        int new_size( int(std::max(50, i_store)*1.10) );
        new_seg.reserve(new_size);
        new_data.reserve(new_size);
      }
      new_seg.push_back(i);
      new_data.push_back(MapNcOutput(&nc, min_val, max_val, output[1]));
      ++i_store;

  } // while (max(clock_input_idx, cone_input_idx) < last_clock_input_idx)

  if (i>last_input_idx)
  {
    if ((new_seg[i_store-1]-new_seg[i_store-2])<min_seg_size)
    {
      new_seg.pop_back();
      new_data.pop_back();
      new_seg.back() = last_input_idx;
    }
  }

}

// Check for fitness object and set new one if check fails. Usually this is
// the case for freshly created individuals which haven't been evaluated yet.
//-----------------------------------------------------------------------------
inline bool CheckFitnessObjectOfApproxIndividual(CIndividual& individual)
{
  using namespace std;
  using namespace approx;

  if (individual.hasFitnessObject() != true)
  {
    try
    {
      individual.setFitnessObject(new CApproxFitness);
    }
    catch (std::bad_alloc)
    {
      const string method_name("bool CheckFitnessObjectOfApproxIndividual(CIndividual& individual):");
      const string warn_msg("solution has no fitness object and allocation of a new one failed!\n");
      Warning(method_name + warn_msg);
      return false;
    }
  }
  return true;
}



//-----------------------------------------------------------------------------
CEAFitness* CSimulation::evaluateApproximationSolution(
         CIndividual& solution
  ,const bool         final
  ,const bool         monitor
  ,const double       time)
{
  using namespace approx;

  if (CheckFitnessObjectOfApproxIndividual(solution) == false) return 0;

  // load approximator (neural net)
  m_nc->load(solution.getSolution(), 0, m_nc->getNumberOfParameters()-1);
  // load approximator 2 (neural net)
  m_nc2->load(solution.getSolution(), m_nc->getNumberOfParameters(),
    m_nc->getNumberOfParameters()+m_nc2->getNumberOfParameters()-1);

  // maximum # of clock angle value segments
  const unsigned int max_clock_segs(m_ref_clock_angle.size()-1);
  // maximum # of cone angle value segments
  const unsigned int max_cone_segs(m_ref_cone_angle.size()-1);

  // check for equal input array dimensions
  if (max_clock_segs != max_cone_segs)
  {
    Warning("CEAFitness* CSimulation::evaluateApproximationSolution(...):\
            unequal numbers of data values");
    return 0; // return empty fitness object
  }

  // output data structures of clock/cone angle segment end indices
  //-------------------------------------------------------------------
  // Indices in original angle data structures are stored within the
  // following two data structures, eg. an index value of 23 at
  // position 5 means that the end of segment 6 is marked by the 23th
  // value of original data set. Thus every index entry at position i
  // must be greater than the index entries at positions i-n (n=1...i).
  vector<unsigned int> clock_segmentation(max_clock_segs);
  vector<unsigned int> cone_segmentation(max_cone_segs);

  // data structures of modification values
  //-------------------------------------------------------------------
  // Besides suggesting segment end indices the neural net approximator
  // also gives a modification value. This value has to be added to
  // the original angle value to obtain the new clock/cone angle value
  // at segment index position. This new angle value will be taken for
  // spline - interpolation.
  vector<double> mod_clock(max_clock_segs, 0.0);     // clock angle modification
  vector<double> mod_cone(max_cone_segs, 0.0);        // coneangle modification

  // maximal number of input data values per scan operation
  const unsigned int max_input_values(m_nc->getNumberOfInputNeurons());

  // total # of input values (clock/cone angle)
  const unsigned int n_input_values = m_ref_clock_angle.size();
  // last input indices to process
  const unsigned int last_clock_input_idx(n_input_values-1);      // clock angle
  const unsigned int last_cone_input_idx(last_clock_input_idx);    // cone angle

  const unsigned int min_seg_points(3);          // min. # of points per segment
  const unsigned int min_seg_size(min_seg_points);       // minimum segment size

  // actual indices in clock/cone angle input data structures
  unsigned int clock_input_idx(0), cone_input_idx(0);
  unsigned int input_values(0);   // current # of input values (single scan op.)

  Vector input(max_input_values);                     // neural net input vector
  Vector output(m_nc->getNumberOfOutputNeurons());   // neural net output vector

  unsigned int clock_segment(0), cone_segment(0);      // actual segment indices
  unsigned int segment_size(0);                            // temporary variable
  unsigned int last_clock_seg_end(0), last_cone_seg_end(0); // last index
  unsigned int max_seg(0);

  const double init_val(0.0);            // nc input vector initialization value

  std::vector<double> data,new_time_steps;
  if (m_nc->getType()>0)
  {
    applyApproxNC1(*m_nc, m_ref_time, 0, m_ref_time.size()-1, m_ref_clock_angle,
      clock_segmentation, m_best_clock_knot_values,0, -PI, THREEPI);

    clock_segment = clock_segmentation.size();

    applyApproxNC1(*m_nc2, m_ref_time, 0, m_ref_time.size()-1, m_ref_cone_angle,
      cone_segmentation, m_best_cone_knot_values,0, 0, HALFPI);

    cone_segment = cone_segmentation.size();
  }
  else
  {

    //applyApproxNC2(*m_nc, m_ref_time, 0, m_ref_time.size()-1, m_ref_clock_angle,
    //  clock_segmentation, m_best_clock_knot_values,0, -PI, THREEPI);

    //clock_segment = clock_segmentation.size();

    //applyApproxNC2(*m_nc2, m_ref_time, 0, m_ref_time.size()-1, m_ref_cone_angle,
    //  cone_segmentation, m_best_cone_knot_values,0, 0, HALFPI);

    //cone_segment = cone_segmentation.size();

    // determination of clock/cone angle knot points (for spline interpolation)
    // and modification values with neural net
    while (clock_input_idx < last_clock_input_idx ||
      cone_input_idx  < last_cone_input_idx)
    {
      // first, clock angle segmentation estimation
      if (clock_input_idx < last_clock_input_idx)           // already at end ?
      {
        // no, so determine number of next input values
        input_values = std::min(n_input_values-clock_input_idx,max_input_values);
        input.Set(init_val);                               // reset input vector

        for (uint i=0;i<input_values;++i)                   // load input vector
          input[i] = m_ref_clock_angle[clock_input_idx+i];

        m_nc->propagate(input, output);                  // applicate neural net

        // calculate next knot point
        segment_size = 1+uint(MapNcOutput(m_nc, 1.0, double(input_values), output[0]));
        clock_segmentation[clock_segment]=(--segment_size)+last_clock_seg_end;

        // calculate knot modification value
        mod_clock[clock_segment] =
          MapNcOutput(m_nc, -TWOPI, TWOPI, output[1]);

        // save last used index in original data structure and
        // new start index for input values calculation
        last_clock_seg_end = clock_input_idx = clock_segmentation[clock_segment++];
      } // if (clock_input_idx < last_clock_input_idx)

      // second, do the same with cone angle values
      if (cone_input_idx < last_cone_input_idx)
      {
        input_values = std::min(n_input_values-cone_input_idx,max_input_values);
        input.Set(init_val); // reset input vector
        for (uint i=0;i<input_values;++i)
          input[i]=m_ref_cone_angle[cone_input_idx+i];
        m_nc2->propagate(input, output);
        segment_size = 1+uint(MapNcOutput(m_nc2, 1.0, double(input_values), output[0]));
        cone_segmentation[cone_segment]=(--segment_size)+last_cone_seg_end;
        mod_cone[cone_segment] =
          MapNcOutput(m_nc2, -HALFPI, HALFPI, output[1]);
        last_cone_seg_end = cone_input_idx = cone_segmentation[cone_segment++];
      } // if (cone_input_idx < last_cone_input_idx)

    } // while (max(clock_input_idx, cone_input_idx) < last_clock_input_idx)

    //in case of uncorrect segmentation force guilty end index
    if (clock_input_idx>last_clock_input_idx)
    {
      if ((clock_segmentation[clock_segment-1]-clock_segmentation[clock_segment-2])<min_seg_size) --clock_segment;
      clock_segmentation[clock_segment-1]=last_clock_input_idx;
    }
    if (cone_input_idx>last_cone_input_idx)
    {
      if ((cone_segmentation[cone_segment-1]-cone_segmentation[cone_segment-2])<min_seg_size) --cone_segment;
      cone_segmentation[cone_segment-1]=last_cone_input_idx;
    }
  }

  if (m_nc->getType() == 0)
  {
    // resize modification value vectors
    mod_clock.resize(clock_segment);
    mod_cone.resize(cone_segment);
  }
  // get application specific fitness object from individual and reset it
  CApproxFitness* f = dynamic_cast<CApproxFitness*>(solution.getFitnessObject());
  f->reset(-HUGE_VAL);

  // set clock/cone angle segmentation subfitness
  const double clock_seg_fit(1.0-calcRatioSubFitness2(clock_segment,uint(2),
    n_input_values));
  f->setFitness(CLOCK_SEG,clock_seg_fit);
  const double cone_seg_fit(1.0-calcRatioSubFitness2(cone_segment,uint(2),
    n_input_values));
  f->setFitness(CONE_SEG,cone_seg_fit);

  // resize knots data structures to correct size
  clock_segmentation.resize(clock_segment);
  cone_segmentation.resize(cone_segment);

  // calculate approximation of steering angles
  // first, clock angle
  std::valarray<unsigned int> segs;
  CopyDataStructure(clock_segmentation,segs);

  ResampleArray(m_ref_time, new_time_steps, m_approx_res_factor); // get new time steps

  // get knot time and modified angle value from reference data structures and
  // prepare spline-interpolation input data structures
  PrepareArray(m_ref_time,m_ref_clock_angle
    ,segs,m_best_clock_knot_time,m_best_clock_knot_values);

  uint i(0);

  if (m_nc->getType() == 0)
  {
    // carry out modification
    for (i=0;i<mod_clock.size();++i) m_best_clock_knot_values[i+1]+=mod_clock[i];
  }

  // clock angle history approximation through spline interpolation of
  // certain, modified knot values taken from original angle history
  SplineApproximation(m_best_clock_knot_time, m_best_clock_knot_values
    , new_time_steps, m_last_saved_clock_angle_approx_values
    , m_last_saved_clock_angle_approx_values_dot);

  // second, cone angle
  CopyDataStructure(cone_segmentation,segs);

  PrepareArray(m_ref_time,m_ref_cone_angle,segs,m_best_cone_knot_time
    ,m_best_cone_knot_values);

  if (m_nc->getType()==0)
    for (i=0;i<mod_cone.size();++i) m_best_cone_knot_values[i+1]+=mod_cone[i];

  ResampleArray(m_ref_time, new_time_steps, m_approx_res_factor);
  SplineApproximation(m_best_cone_knot_time, m_best_cone_knot_values, new_time_steps
    , m_last_saved_cone_angle_approx_values
    , m_last_saved_cone_angle_approx_values_dot);

  // set clock/cone angle sign change sub fitness of approximation as a
  // measurement for smoothness of new steering angle history
  const double clock_scc(calcRatioSubFitness(1.0
    ,CountSignChangesDerivative(m_last_saved_clock_angle_approx_values_dot), n_input_values));
  f->setFitness(approx::CLOCK_SC, clock_scc );

  const double cone_scc(calcRatioSubFitness(1.0
    ,CountSignChangesDerivative(m_last_saved_cone_angle_approx_values_dot), n_input_values));
  f->setFitness(approx::CONE_SC, cone_scc);

  // assess steering angles max/min solution qualities
  m_last_saved_trajectory_approx_steps=0;
  i=0;
  double cla(0.0), coa(0.0);
  double max_clock(0.0), max_cone(0.0), min_clock(0.0), min_cone(0.0);
  const unsigned int sz(m_last_saved_clock_angle_approx_values.size());
  while (i<sz)
  {
    cla = m_last_saved_clock_angle_approx_values[i];
    coa = m_last_saved_cone_angle_approx_values[i];
    max_clock = std::max(max_clock, cla);
    min_clock = std::min(min_clock, cla);
    max_cone = std::max(max_cone, coa);
    min_cone = std::min(min_cone, coa);
    ++i;
  }
  double max_clock_fit( max_clock>THREEPI ? (THREEPI/max_clock) : 1.0);
  double min_clock_fit( min_clock<(-HALFPI) ? -HALFPI/min_clock : 1.0);
  double max_cone_fit( max_cone>HALFPI ? (HALFPI/max_cone) : 1.0);
  double min_cone_fit( min_cone<0.0 ? pow(10.0,min_cone) : 1.0);
  // set max/min steering angles fitness values
  f->setFitness(approx::MAX_CLOCK,max_clock_fit);
  f->setFitness(approx::MIN_CLOCK,min_clock_fit);
  f->setFitness(approx::MAX_CONE,max_cone_fit);
  f->setFitness(approx::MIN_CONE,min_cone_fit);

  if ((max_clock_fit+min_clock_fit+max_cone_fit+min_cone_fit)<4.0) return 0;


  // evaluate approximation of steering angles by EOM integration
  evaluateSteeringAngleApprox(
      solution
    , m_last_saved_clock_angle_approx_values
    , m_last_saved_cone_angle_approx_values
    , new_time_steps
    , m_last_saved_trajectory_approx
    , m_last_saved_trajectory_approx_steps
    , m_ref_sc_traj[0].MJD()
    , final
    , monitor
  );

  // set original optimization fitness
  f->setFitness(approx::OPT_FITNESS, m_approx_fitness.getMixedFitness());

  // if it is the final evaluation of best individual then save
  // the resulting s/c trajectory in certain member field to use
  // existing code for data file output
  if (final==true)
  {
    m_spacecraft_trajectory = m_last_saved_trajectory_approx;
    m_steps_best = m_last_saved_trajectory_approx_steps-1;
//    auto_ptr<CTrajectory> traj(getLastInertialRotStateTrajectory());
//    m_ref_sc_rot_traj.resize(traj->size(), Vector(8));
//    for (uint i=0;i<traj->size();++i) m_ref_sc_rot_traj[i] = (*traj)[i];
    m_ref_time.resize(m_steps_best);
    for (int i=0;i<m_steps_best;++i) m_ref_time[i]=new_time_steps[i];
  }

  return 0;  // return final fitness object
}


//------------------------------------------------------------------------------
void PrintCsvVectorHeader(
       std::ofstream& ofs,
 const std::string& prefix,
 const Vector& v)
{
  const int n(v.size());  // # of input neurons
  int i(0);
  std::ostringstream oss("");
  oss << dec << left << setw(0);
  while (i<n)
  {
    if (oss.bad() == true) oss.clear();
    oss.str("");
    oss << prefix << "(" << i << ")";
    if ( i<(n-1) ) oss << ", ";
    ofs << oss.str();
    ++i;
  }
}


//-----------------------------------------------------------------------------
void EOM_NEA_RK(double t, const Vector & y, Vector & yp, void* pAux)
{
  CSolarSystem* s;
  s = static_cast<CSolarSystem*>(pAux);

  Vector rsc(y(0), y(1), y(2));

  Vector rS = s->Cart(SUN     ,t).Slice(0,2) - rsc;
  Vector rR = s->Cart(MERCURY ,t).Slice(0,2) - rsc;
  Vector rV = s->Cart(VENUS   ,t).Slice(0,2) - rsc;
  Vector rE = s->Cart(EARTH   ,t).Slice(0,2) - rsc;
  Vector rM = s->Cart(MARS    ,t).Slice(0,2) - rsc;
  Vector rJ = s->Cart(JUPITER ,t).Slice(0,2) - rsc;
  Vector rA = s->Cart(SATURN  ,t).Slice(0,2) - rsc;
  Vector rU = s->Cart(URANUS  ,t).Slice(0,2) - rsc;
  Vector rN = s->Cart(NEPTUNE ,t).Slice(0,2) - rsc;

  const double muS = 1.32713440018e20*GM_SI2GM_AUD;

  Vector gS = (muS          /pow(Norm(rS),3)) * rS;
  Vector gR = (muS/6023600.0/pow(Norm(rR),3)) * rR;
  Vector gV = (muS/408523.71/pow(Norm(rV),3)) * rV;
  Vector gE = (muS/328900.56/pow(Norm(rE),3)) * rE;
  Vector gM = (muS/3098708.0/pow(Norm(rM),3)) * rM;
  Vector gJ = (muS/1047.3486/pow(Norm(rJ),3)) * rJ;
  Vector gA = (muS/3497.898 /pow(Norm(rA),3)) * rA;
  Vector gU = (muS/22902.98 /pow(Norm(rU),3)) * rU;
  Vector gN = (muS/19412.24 /pow(Norm(rN),3)) * rN;

  Vector g = gS + gR + gV + gE + gM + gJ + gA + gU + gN;

  yp = Vector(y(3), y(4), y(5), g(0), g(1), g(2));
};

//-----------------------------------------------------------------------------
double CSimulation::calcDeflection(
    const double t0
   ,const double t2
   ,const Vector& nea_r0
   ,const Vector& nea_v0
   ,const Vector& sc_v0
   ,const double impact_eff
)
{
  // bigger step size for integration far away from celestial
  const double step_size_1(1.0);
  // smaller step size for integration in vincinity of Earth or Mars
  const double step_size_2(0.01);

  // asteroid velocity vector after impact
  Vector nea_v1(nea_v0+(sc_v0-nea_v0)*impact_eff);

  // initialize integrator objects
  RKF54 orbit_a(&EOM_NEA_RK, 6, m_solsys);
  orbit_a.init(10000, 1e-12, 1e-12);
  RKF54 orbit_b(&EOM_NEA_RK, 6, m_solsys);
  orbit_b.init(10000, 1e-12, 1e-12);

  // begin of integration for both trajectories and initial step size
  double t_a(t0), t_b(t0), h(step_size_1);

  // state vectors of both trajectories
  Vector ys_a(Stack(nea_r0, nea_v0));   // undisturbed one
  Vector ys_b(Stack(nea_r0, nea_v1));   // after hit from s/c

  // position vector of s/c on trajectory 1
  Vector r_a(3);

  while (t_a < t2 && t_b == t_a)
  {
    orbit_a.Step(t_a, ys_a, h);    // integration step on trajectory 1
    orbit_b.Step(t_b, ys_b, h);    // the same for trajectory 2

    // determine step size depending on current distance to Earth and Mars
    h = (Norm(m_solsys->Cart(MARS ,t_a).Slice(0,2)-(r_a=ys_a.Slice(0,2))) > 0.2
      && Norm(m_solsys->Cart(EARTH,t_a).Slice(0,2)-r_a) > 0.1) ? step_size_1 : step_size_2;
  }

  // return final relative distance between s/c on trajectory 1 and 2
  return Norm(r_a - ys_b.Slice(0,2))*AU2KM;//result;
} // double CSimulation::calcDeviation(....)


//-----------------------------------------------------------------------------
void CSimulation::checkSensitivityFactorsInitialization(void)
{
  switch (jt)
  {
  case J_AND:
  case J_OR:
  case J_RIMPACT:
  case J_RIMPACT_AND:
  case J_RIMPACT_OR:
  case J_ACT:
  case J_SPI:
  case J_ORB_ELEMENTS:
  case J_GTOCOMP2:
  case J_TEST_1:
  case J_TEST_2:
  case J_TEST_3:
  case J_TEST_4:
  case J_TEST_5:
  case J_TEST_6:
  case J_TEST_7:
  case J_TEST_8:
  case J_TEST_9: break;
  default:
    if (init_sensitivity_factors == true) initSensitivityFactors();
    init_sensitivity_factors=false;
    break;
  }
} // void CSimulation::checkSensitivityFactorsInitialization(void)

//-----------------------------------------------------------------------------
bool CSimulation::continueStandardEvaluationIntegrationLoop(bool final,
                                                            bool monitor)
{
  using namespace std;
  using namespace constraints;

  if ( (shouldAbort()==true) && (final == false) && (monitor==false) ) 
    return false;
  if (getCurrentStepCount() >= m_steps_max_used)
    return false;
  if ( getFlightTime() >= getIntegrationInterval())
    return false;

  vector<CConstraint*>::iterator constraint = m_constraints.begin();
  while (constraint!= m_constraints.end())
  {
    if ((*constraint)->isMatched() == false) return false;
    ++constraint;
  }
  return true;
} // bool CSimulation::continueStandardEvaluationIntegrationLoop(...)


//-----------------------------------------------------------------------------
void CSimulation::checkDynamicMaxIntegrationTimeReduction(void)
{
  switch (tt)
  {
  case RENDEZVOUS:
    switch (getTargetStateType())
    {
    case BODY_STATE:
    case BODY_ORBIT_STATE:
    case FREE_ORBIT_STATE:
      if (getCurrProximity().CheckProximity(m_max_dr_AU, m_max_dv_AUd))
      {
        m_time_max_dynamic = std::min(m_time_max_dynamic, m_time_best);
      }
      break;
    case ORBITAL_ELEMENTS_STATE:
      if (m_best_prox_sum==0.0)
        m_time_max_dynamic = std::min(m_time_max_dynamic, m_time_best*1.01);
      break;

    default:
      break;
    } // switch (getTargetStateType())
    break;

  case FLYBY:
    switch (getTargetStateType())
    {
    case BODY_STATE:
    case BODY_ORBIT_STATE:
    case FREE_ORBIT_STATE:
      if (getCurrProximity().CheckProximity(m_max_dr_AU, m_max_dv_AUd))
      {
        m_time_max_dynamic = std::min(m_time_max_dynamic, m_time_best*1.01);
      }
      break;

    case ORBITAL_ELEMENTS_STATE:
      break;

    default:
      break;
    } // switch (getTargetStateType())

  case RIMPACT:
    break;

  case ACTPROBLEM:
    break;

  case SPIPROBLEM:
    break;

  default:
    break;
  } // switch (tt)
} // void checkDynamicMaxIntegrationTimeReduction(void)


//-----------------------------------------------------------------------------
double CalcElementChangeTime(const KEPLERIAN_ELEMENT el,
                             const double target_value,
                             CSC& sc,
                             const double t,
                             const double max_diff,
                             const double max_thrust
                             )
{
  double S = 0, d = 0, r = 0, current_value = 0;

  switch (el)
  {
  case SEMIMAJORAXIS: current_value = sc.a(t); break;
  case ECCENTRICITY: current_value = sc.e(t); break;
  case INCLINATION: current_value = sc.i(t); break;
  case NODE: current_value = sc.Omega(t); break;
  case PERIAPSIS: current_value = sc.omega(t); break;
  default:
    assert(false); // raise assertation failed exception for debugging
    break;
  }

  S = CalcScalingFactor(el, current_value, target_value);
  d = CalcDistanceFunction(el, current_value, target_value, max_diff);
  r = CalcMaxElementRateOfChange(el, sc.State(t));

  double factor = (S * fabs(d)/(r*max_thrust));
  return factor;
} // double CalcElementChangeTime(...)


//-----------------------------------------------------------------------------
bool CSimulation::checkForNewBestStepTargetOrbitalElements(void)
{
  double t = 0, ref_t = 0;
  const double t_left = getRemainingFlightTime();

  double max_thrust_acc_factor = m_spacecraft.MaxThrustAcceleration()*MS22AUD2;
  if (m_spacecraft.isSolarSail() == true)
  {
    MathLib::Vector s(m_solsys->Cart(m_spacecraft.getId(), SUN, m_t,AUD));
    const double r_sun_square(pow(s(0),int(2))+pow(s(1),int(2))+pow(s(2),int(2)));
    max_thrust_acc_factor *= pow(1.0/m_spacecraft.r(m_t),int(2));
  }

  const bool time_to_go_tactic = true;

  double current_val = 0.0;
  for (size_t i = SEMIMAJORAXIS; i <= ANOMALY; ++i)
  {
    const KEPLERIAN_ELEMENT e = KEPLERIAN_ELEMENT(i);
    
    if (m_data.useTargetOrbitalElement(e)==true)
    {
      const double max_delta = m_data.getTargetOrbitalElementMaxDeviation(e);
      /*const double target_val = m_data.getTargetOrbitalElement(e);
      switch (e)
      {
      case SEMIMAJORAXIS:current_val = getCurrentSC().CState::a();break;
      case ECCENTRICITY :current_val = getCurrentSC().CState::e();break;
      case INCLINATION  :current_val = getCurrentSC().CState::i();break;
      case NODE         :current_val = getCurrentSC().CState::Omega();break;
      case PERIAPSIS    :current_val = getCurrentSC().CState::omega();break;
      case ANOMALY      :current_val = getCurrentSC().CState::M();break;
      default:break;
      }*/

      //getTOEElementDeviation()[e] =
      //  fabs(CalcScalingFactor(e, current_val, target_val)*
      //  CalcDistanceFunction(e, current_val, target_val, max_delta))
      //  /getTOEInitialElementDeviation()[e];

      //if (time_to_go_tactic == true)
      //{
      //  getTOERemChangeTime()[e] = 
      //    CalcElementChangeTime(e, target_val, m_spacecraft, m_t, max_delta
      //    , max_thrust_acc_factor);
      //}

      switch (e)
      {
        case SEMIMAJORAXIS:
          getTOEElementDeviation()[e]=
            getCurrProximity().getNormedSemimajorAxisProximity(max_delta);
          break;
        case ECCENTRICITY:
          getTOEElementDeviation()[e]=
            getCurrProximity().getNormedEccentricityProximity(max_delta);
          break;
        case INCLINATION:
          getTOEElementDeviation()[e]=
            getCurrProximity().getNormedInclinationProximity(max_delta);
          break;
        case NODE:
          getTOEElementDeviation()[e]=
            getCurrProximity().getNormedNodeProximity(max_delta);
          break;
        case PERIAPSIS:
          getTOEElementDeviation()[e]=
            getCurrProximity().getNormedPericenterProximity(max_delta);
          break;
        case ANOMALY:
          getTOEElementDeviation()[e]=
            getCurrProximity().getNormedMeanAnomalyProximity(max_delta);
          break;
        default:
          break;
      }
    } else getTOEElementDeviation()[e] = 0.0;
  }

  if (time_to_go_tactic == true) 
  {
    int j;
    for (j=SEMIMAJORAXIS;j<=ANOMALY;++j)
    {
      const KEPLERIAN_ELEMENT e = KEPLERIAN_ELEMENT(j);
      if (useTargetOrbitalElement(e)==true)
      {
        double s = 0.0;
        int k;
        for (k=SEMIMAJORAXIS; k<=ANOMALY;++k)
        {
          const KEPLERIAN_ELEMENT f = KEPLERIAN_ELEMENT(k);       
          if (e!=f && (useTargetOrbitalElement(f)==true))
          {
            s += getTOEElementDeviation()[e]/
                std::max(1.0e-3,getTOEElementDeviation()[f]);
          }
        } // for (k=SEMIMAJORAXIS; k<=ANOMALY;++k)
        getTargetOrbitalElementProxWeights()[e] = s;
      } // if (useTargetOrbitalElement(e)==true)
    } // for (j=SEMIMAJORAXIS;j<=ANOMALY;++j)

    for (j=SEMIMAJORAXIS;j<=ANOMALY;++j)
    {
      if (getTOEElementDeviation()[j]>1.0)
      {
       t+=(getTOEElementDeviation()[j]*getTargetOrbitalElementProxWeights()[j]);
      }
    }
    t /= getTargetOrbitalElementProxWeights().sum();

    if (t > 0.0)
    {
      //const double rest_time_elements = 
      //  (getTOERemChangeTime()*getTargetOrbitalElementProxWeights()).sum()
      //  /(getTargetOrbitalElementProxWeights().sum());
 
      t *= (
        + pow(getMaxFlightTime()
         / std::max(getMaxFlightTime()/2,
             std::min(getMaxFlightTime()/2,getRemainingFlightTime())), int(2))
      );
    }

  }
   
  if (t<m_best_prox_sum)
  {
    m_best_prox_sum  = t;
    //m_steps_best     = getCurrentStepCount();
    //m_time_best      = getFlightTime();
    //m_best_proximity = getCurrProximity();
    ////setTOEBestRemChangeTime(getTOERemChangeTime());
    //m_best_sc        = getCurrentSC();
    //setTOEBestElementDeviation(getTOEElementDeviation());

    //if (getCurrentSC().isSolarSail() == false)
    //{
    //  m_best_mp = calcPropMassFitness(
    //    getCurrentSC().getInitialPropellantMass(),
    //    getCurrentSC().getPropellantMass());
    //}
    return true;
  }
  else return false;
} // void CSimulation::checkForNewBestStepTargetOrbitalElements(void)


//-----------------------------------------------------------------------------
bool CSimulation::checkForNewBestStepGTOC2(void)
{
  if (m_temp_prox_sum > m_best_prox_sum)
  {
    getCurrProximity().SetProximity(*m_gtoc2_temp_asteroid, getCurrentSC());
    calcBestValues();
    return true;
  } // if (m_temp_prox>m_best_prox_sum)
  return false;
} // void CSimulation::checkForNewBestStepTest0(void)


//-----------------------------------------------------------------------------
void CSimulation::checkForNewBestStepTest1(void)
{
} // void CSimulation::checkForNewBestStepTest1(void)

//-----------------------------------------------------------------------------
void CSimulation::checkForNewBestStepTest2(void)
{
} // void CSimulation::checkForNewBestStepTest2(void)

//-----------------------------------------------------------------------------
void CSimulation::checkForNewBestStepTest3(void)
{
} // void CSimulation::checkForNewBestStepTest3(void)

//-----------------------------------------------------------------------------
void CSimulation::checkForNewBestStepTest4(void)
{
} // void CSimulation::checkForNewBestStepTest4(void)

//-----------------------------------------------------------------------------
void CSimulation::checkForNewBestStepTest5(void)
{
} // void CSimulation::checkForNewBestStepTest5(void)

//-----------------------------------------------------------------------------
void CSimulation::checkForNewBestStepTest6(void)
{
} // void CSimulation::checkForNewBestStepTest6(void)

//-----------------------------------------------------------------------------
void CSimulation::checkForNewBestStepTest7(void)
{
} // void CSimulation::checkForNewBestStepTest7(void)

//-----------------------------------------------------------------------------
void CSimulation::checkForNewBestStepTest8(void)
{
} // void CSimulation::checkForNewBestStepTest8(void)

//-----------------------------------------------------------------------------
void CSimulation::checkForNewBestStepTest9(void)
{
} // void CSimulation::checkForNewBestStepTest9(void)


//-----------------------------------------------------------------------------
void CSimulation::checkForNewBestStepRetrogradeImpact(void)
{
  m_steps_best = getCurrentStepCount();
  m_time_best = (m_t-getStartTime());
  m_best_proximity = getCurrProximity();
  calcBestValues();
  if (getOptGoalType() == MAXIMUM_DEFLECTION)
  {
    m_best_s = m_spacecraft.CState::Cart();
  }
} // void CSimulation::checkForNewBestStepRetrogradeImpact(void)


//-----------------------------------------------------------------------------
void CSimulation::checkForNewBestStepRendezvous(void)
{
  bool new_best_step = false;
  switch (getTargetStateType())
  {
  case ORBITAL_ELEMENTS_STATE:
    new_best_step = checkForNewBestStepTargetOrbitalElements();
    break;

  default:
    // Body rendezvous stops at a feasible current step; retain that endpoint.
    // Other mission/state routes and same-feasibility ordering stay unchanged.
    if (tt == RENDEZVOUS && getTargetStateType() == BODY_STATE
      && getCurrProximity().CheckProximity(m_max_dr_AU, m_max_dv_AUd) !=
        m_best_proximity.CheckProximity(m_max_dr_AU, m_max_dv_AUd))
    {
      new_best_step = getCurrProximity().CheckProximity(m_max_dr_AU, m_max_dv_AUd);
    }
    else if (getCurrProximity().GetNormedProximity(m_max_dr_AU, m_max_dv_AUd) <
      m_best_proximity.GetNormedProximity(m_max_dr_AU, m_max_dv_AUd))
    {
      new_best_step = true;
    }
    break;
  }
  if (new_best_step == true) calcBestValues();
} // void CSimulation::checkForNewBestStepRendezvous(void)


//-----------------------------------------------------------------------------
void CSimulation::checkForNewBestStepACTProblem(void)
{
  const double delta_r = getCurrProximity().GetDistanceProximity();
  const double m_act = AUD2KMS*AUD2KMS*m_spacecraft.getMass()*
    getCurrProximity().GetActVelocity(m_spacecraft, 
     m_solsys->State(getTargetBodyId(), m_t));
  const double f_r = std::max(0.0, (delta_r-m_max_dr_AU))/m_max_dr_AU;

  //double t = - (m_act - 
  //  pow(Global::max(0.0, (delta_r-m_max_dr_AU))*m_data.act_distance_factor
  //  , m_data.act_exponent));
  double t = -(
      m_act*(1.0/(1.0+f_r)) 
      + (m_act*0.01*(1.0- std::min(1.0,delta_r/m_max_dr_AU)))
    - (f_r/(1.0+f_r))*pow((f_r*m_max_dr_AU*m_data.act_distance_factor),m_data.act_exponent));

  if (t<m_best_prox_sum)
  {
    m_best_prox_sum = t;
    m_best_dr = delta_r;
    calcBestValues();
  }
}


//-----------------------------------------------------------------------------
void CSimulation::checkForNewBestStepSPIProblem(void)
{
  m_steps_best = getCurrentStepCount();
  m_time_best = getFlightTime();
  calcBestValues();
} // void CSimulation::checkForNewBestStepSPIProblem(void)


//-----------------------------------------------------------------------------
bool CSimulation::checkForNewBestStep(void)
{
  switch (tt)  //  trajectory type
  {
  case RIMPACT: checkForNewBestStepRetrogradeImpact();break;
  case ACTPROBLEM: checkForNewBestStepACTProblem();break;
  case SPIPROBLEM: checkForNewBestStepSPIProblem();break;
  case RENDEZVOUS:checkForNewBestStepRendezvous();break;
  case FLYBY:checkForNewBestStepRendezvous();break;
  case GTOC2TRAJ: return checkForNewBestStepGTOC2();break;
  case TRAJ_TEST1:checkForNewBestStepTest1();break;
  case TRAJ_TEST2:checkForNewBestStepTest2();break;
  case TRAJ_TEST3:checkForNewBestStepTest3();break;
  case TRAJ_TEST4:checkForNewBestStepTest4();break;
  case TRAJ_TEST5:checkForNewBestStepTest5();break;
  case TRAJ_TEST6:checkForNewBestStepTest6();break;
  case TRAJ_TEST7:checkForNewBestStepTest7();break;
  case TRAJ_TEST8:checkForNewBestStepTest8();break;
  case TRAJ_TEST9:checkForNewBestStepTest9();break;
  default:
    Warning("void CSimulation::checkForNewBestStep(void): \
            Unknown trajectory type!\n");
    break;
  } // end switch
  return false;
} // void CSimulation::checkForNewBestStep(void)


//-----------------------------------------------------------------------------
double CSimulation::calcDSSCStepSize(void)
{
  double h = 1.0e10;

  // orbit angular velocity control ?
  if (getStatusDSSCStepAngleControl() == true)
  {
    const double oav = sqrt(pow(m_spacecraft.CState::dphi(),int(2))+
      pow(m_spacecraft.CState::dtheta(),int(2)));
    h = getAccuracyFactor()*getCurrMaxDSSCOrbAngle()/ std::max(1.0e-7,oav);
  }

  // total velocity step size control ?
  if (getStatusDSSCStepDistanceControl() == true)
  {
    h = std::min(h, getDSSCMaxStepDistance()/ std::max(1.0e-7
      , m_spacecraft.CState::Velocity()));
  }

  // target approach step size control ?
  if (getStatusDSSCApproachControl()== true)
  {
    double _h=getCurrProximity().GetDistanceProximity()*getDSSCApproachFactor();
    if ( _h <= getDSSCMinStepSize() )
    {
      if (m_SpacecraftOrbit_RKF54 != 0)
        m_SpacecraftOrbit_RKF54->setStopFunc(IntegrationStopFunction, this);
    } // if ( _h <= getDSSCMinStepSize() )
    else
    {
      if (m_SpacecraftOrbit_RKF54!=0) m_SpacecraftOrbit_RKF54->resetStopFunc();
    } // else branch
    h = std::min(h, _h);
  } // if (getStatusDSSCApproachControl()== true)

  return MathLib::bounded(h, 
    getDSSCMinStepSize(), getAccuracyFactor()*getDSSCMaxStepSize());
}


//-----------------------------------------------------------------------------
double CSimulation::calcIntegrationStepWidth(const double h_min, bool final)
{
  double h = 1.0e10; // days

  if (getStatusDSSC() == true)
  {
    h = std::min(calcDSSCStepSize(), getRemainingFlightTime());
  } else h = (getIntegrationInterval()/ std::max(1, m_steps_max_used));

  switch (tt)
  {
  case RENDEZVOUS:
    return h;
    break;
  case FLYBY:
    if (m_data.number_of_ga_bodies == 0
      && m_spacecraft.isSolarSail()
      && m_spacecraft.CState::r() >= SAILJETTISONDISTANCE
      && m_time_max_dynamic-(m_t-getStartTime()) > SAILJETTISONINTEGRATIONSTEPWIDTH
      && (m_solsys->r(getTargetBodyId(), m_t)-m_spacecraft.CState::r() >= 2.0*m_spacecraft.CState::dr()*SAILJETTISONINTEGRATIONSTEPWIDTH
      || 200.0-m_spacecraft.CState::r() >= 2.0*m_spacecraft.CState::dr()*SAILJETTISONINTEGRATIONSTEPWIDTH))
    {
      return SAILJETTISONINTEGRATIONSTEPWIDTH;
    }
    else return h;
    break;
  case RIMPACT:
    return h;
    break;
  case ACTPROBLEM:
    return h;
    break;
  case SPIPROBLEM:
    return h;
    break;
  default:    
    return h;
    break;
  }
} // double CSimulation::calcIntegrationStepWidth(const double h_min,bool final)


//-----------------------------------------------------------------------------
void CSimulation::setFinalTargetStateOrbitalElements(void)
{
  MathLib::Vector elem(6),delem(6);
  const double gm = m_solsys->getBody(getTargetBodyId())->getCenterState()->getGM(AUD);

  if (m_data.useTargetOrbitalElement(SEMIMAJORAXIS))
    elem(SEMIMAJORAXIS) = m_data.getTargetOrbitalElement(SEMIMAJORAXIS);
  else elem(SEMIMAJORAXIS) = m_best_sc.CState::a();

  if (m_data.useTargetOrbitalElement(ECCENTRICITY))
    elem(ECCENTRICITY) = m_data.getTargetOrbitalElement(ECCENTRICITY);
  else elem(ECCENTRICITY) = m_best_sc.CState::e();

  if (m_data.useTargetOrbitalElement(INCLINATION))
    elem(INCLINATION) = m_data.getTargetOrbitalElement(INCLINATION);
  else elem(INCLINATION) = m_best_sc.CState::i();

  if (m_data.useTargetOrbitalElement(NODE))
    elem(NODE) = m_data.getTargetOrbitalElement(NODE);
  else elem(NODE) = m_best_sc.CState::Omega();

  if (elem(INCLINATION)==0.0) elem(NODE)=0.0; // in case of no inclination there's no node

  if (m_data.useTargetOrbitalElement(PERIAPSIS))
    elem(PERIAPSIS) = m_data.getTargetOrbitalElement(PERIAPSIS);
  else elem(PERIAPSIS) = m_best_sc.CState::omega();

  if (elem(ECCENTRICITY)==0.0) elem(PERIAPSIS)=0.0; // no eccentricity -> no periapsis

  if (m_data.useTargetOrbitalElement(ANOMALY))
    elem(ANOMALY) = m_data.getTargetOrbitalElement(ANOMALY);
  else elem(ANOMALY) = m_best_sc.CState::M();

  if (elem(INCLINATION)==0.0) elem(ANOMALY)+=m_best_sc.CState::Omega();
  if (elem(ECCENTRICITY)==0.0) elem(ANOMALY)+=m_best_sc.CState::omega();

  elem(PERIAPSIS) += elem(NODE);
  elem(ANOMALY) += elem(PERIAPSIS);
  delem(ANOMALY) = CalcMeanMotion(elem(SEMIMAJORAXIS), gm);

  m_solsys->getBody(getTargetBodyId())->setKeplerianElements(elem,
    (m_time_best+getStartTime()));
  m_solsys->getBody(getTargetBodyId())->setKeplerianElementsRate(delem,
    (m_time_best+getStartTime()));
} // void CSimulation::setFinalTargetStateOrbitalElements(void)


//-----------------------------------------------------------------------------
void CSimulation::setFinalTargetState(void)
{
  switch (getTargetStateType())
  {
  case ORBITAL_ELEMENTS_STATE:setFinalTargetStateOrbitalElements();break;
  default:break;
  }
} // void CSimulation::setFinalTargetState(void)


//-----------------------------------------------------------------------------
void CSimulation::initIntegrationVector(void)
{
  switch (tt)
  {
  case RENDEZVOUS:
  case FLYBY:
  case RIMPACT:
  case ACTPROBLEM:
  case SPIPROBLEM:
  case GTOC2TRAJ:
    setIntegrationVectorDimension(8);
    if (m_x.size()!=3) m_x.Set(MathLib::Vector(3));
    if (m_v.size()!=3) m_v.Set(MathLib::Vector(3));
    break;
  case TRAJ_TEST1:break;
  case TRAJ_TEST2:break;
  case TRAJ_TEST3:break;
  case TRAJ_TEST4:break;
  case TRAJ_TEST5:break;
  case TRAJ_TEST6:break;
  case TRAJ_TEST7:break;
  case TRAJ_TEST8:break;
  case TRAJ_TEST9:break;
  default:break;
  }
  if (m_ys.size() != getIntegrationVectorDimension())
    m_ys.Set(MathLib::Vector(getIntegrationVectorDimension()));
} // void CSimulation::setIntegrationVectorDimension(void)

//-----------------------------------------------------------------------------
void CSimulation::setRendezvousInitConditions(CIndividual& i, double time)
{
  m_initial_mp = m_initial_velocity = 0.0;
  m_initial_azimuth = m_initial_elevation = 0.0;

  // start time
  if (time == 0.0)
  {
    m_t = i.getMappedAllele(i.getLength()-5
      , getMinStartTime(), getMaxStartTime(), true);
  }
  else m_t = time;

  // initial s/c propellant mass
  if (m_spacecraft.isNEP()||m_spacecraft.isSEP())
  {
    m_initial_mp = 
      i.getMappedAllele(i.getLength()-1, 0.0, 1.0, true);
  }
  m_spacecraft.Init(m_initial_mp);        // initialize s/c with propellant mass

  // initial excess velocity
  m_initial_velocity  = 
    i.getMappedAllele(i.getLength()-4
    , getMinLaunchVelocity(), getMaxLaunchVelocity());

  // initial azimuth angle of excess velocity vector
  m_initial_azimuth   = 
    i.getCircMappedAllele(i.getLength()-3
    , getMinLaunchAzimuthAngle(), getMaxLaunchAzimuthAngle());

  // initial elevation angle of excess velocity vector
  m_initial_elevation = 
    i.getMappedAllele(i.getLength()-2
    , getMinLaunchElevationAngle(), getMaxLaunchElevationAngle());

  setStartTime(m_t);

  m_steps_max_used = getMaxStepsCount();
  m_h = getIntegrationInterval()/(m_steps_max_used);
  m_h_min = getIntegrationInterval()/getMaxStepsCount();

  m_best_prox_sum = INIT_BEST_PROX_SUM;

  if (getTargetStateType()==ORBITAL_ELEMENTS_STATE)
  {
    if (m_delements.size() != N_KEPLERIAN_ELEMENTS)
      m_delements.resize(N_KEPLERIAN_ELEMENTS, 0.0);
    if (m_best_delements.size() != N_KEPLERIAN_ELEMENTS)
      m_best_delements.resize(N_KEPLERIAN_ELEMENTS, 0.0);
    if (m_element_change_rest_times.size() != N_KEPLERIAN_ELEMENTS)
      m_element_change_rest_times.resize(N_KEPLERIAN_ELEMENTS, 0.0);
    if (m_best_element_change_rest_times.size() != N_KEPLERIAN_ELEMENTS)
      m_best_element_change_rest_times.resize(N_KEPLERIAN_ELEMENTS, 0.0);
  }

} // void CSimulation::setRendezvousInitConditions(CIndividual& i, double time)


//-----------------------------------------------------------------------------
void CSimulation::setGTOC2InitConditions(CIndividual& i, double time)
{
  m_initial_mp = m_initial_velocity = 0.0;
  m_initial_azimuth = m_initial_elevation = 0.0;

  // initial s/c propellant mass
  if (m_spacecraft.isNEP()||m_spacecraft.isSEP())
  {
    m_initial_mp = 
      i.getMappedAllele(i.getLength()-1, 0.0, 1.0);
  }
  m_spacecraft.Init(m_initial_mp);        // initialize s/c with propellant mass

  // start time
  if (time == 0.0)
  {
    m_t = i.getMappedAllele(i.getLength()-5
      , getMinStartTime(), getMaxStartTime(), true);
  }
  else m_t = time;
  setStartTime(m_t);

  // initial excess velocity
  m_initial_velocity  = 
    i.getMappedAllele(i.getLength()-4
    , getMinLaunchVelocity(), getMaxLaunchVelocity(), true);

  // initial azimuth angle of excess velocity vector
  m_initial_azimuth   = 
    i.getCircMappedAllele(i.getLength()-3
    , getMinLaunchAzimuthAngle(), getMaxLaunchAzimuthAngle());

  // initial elevation angle of excess velocity vector
  m_initial_elevation = 
    i.getMappedAllele(i.getLength()-2
    , getMinLaunchElevationAngle(), getMaxLaunchElevationAngle(), true);

  m_steps_max_used = getMaxStepsCount();
  m_h = getIntegrationInterval()/(m_steps_max_used);
  m_h_min = getIntegrationInterval()/getMaxStepsCount();

  m_gtoc2_fp = 1; // start at first flight phase

} // void CSimulation::setRendezvousInitConditions(CIndividual& i, double time)


//-----------------------------------------------------------------------------
void CSimulation::setTrajTest1InitConditions(CIndividual& i, double time)
{
} // void CSimulation::setRendezvousInitConditions(CIndividual& i, double time)


//-----------------------------------------------------------------------------
void CSimulation::setTrajTest2InitConditions(CIndividual& i, double time)
{
} // void CSimulation::setRendezvousInitConditions(CIndividual& i, double time)


//-----------------------------------------------------------------------------
void CSimulation::setTrajTest3InitConditions(CIndividual& i, double time)
{
} // void CSimulation::setRendezvousInitConditions(CIndividual& i, double time)


//-----------------------------------------------------------------------------
void CSimulation::setTrajTest4InitConditions(CIndividual& i, double time)
{
} // void CSimulation::setRendezvousInitConditions(CIndividual& i, double time)


//-----------------------------------------------------------------------------
void CSimulation::setTrajTest5InitConditions(CIndividual& i, double time)
{
} // void CSimulation::setRendezvousInitConditions(CIndividual& i, double time)


//-----------------------------------------------------------------------------
void CSimulation::setTrajTest6InitConditions(CIndividual& i, double time)
{
} // void CSimulation::setRendezvousInitConditions(CIndividual& i, double time)


//-----------------------------------------------------------------------------
void CSimulation::setTrajTest7InitConditions(CIndividual& i, double time)
{
} // void CSimulation::setRendezvousInitConditions(CIndividual& i, double time)


//-----------------------------------------------------------------------------
void CSimulation::setTrajTest8InitConditions(CIndividual& i, double time)
{
} // void CSimulation::setRendezvousInitConditions(CIndividual& i, double time)

//----------------------------------------------------------------------------
void CSimulation::setTrajTest9InitConditions(CIndividual& i, double time)
{
} // void CSimulation::setRendezvousInitConditions(CIndividual& i, double time)

//-----------------------------------------------------------------------------
void CSimulation::setStandardEvaluationInitConditions(CIndividual& solution
                                                      , double time)
{
  // determine simulation start time
  switch (tt)
  {
  case RENDEZVOUS:setRendezvousInitConditions(solution, time);break;
  case FLYBY:setFlyByInitConditions(solution, time);break;
  case RIMPACT:setRetrogradeImpactInitConditions(solution, time);break;
  case ACTPROBLEM:setActProblemInitConditions(solution, time);break;
  case SPIPROBLEM:setSpiProblemInitConditions(solution, time);break;
  case GTOC2TRAJ:setGTOC2InitConditions(solution, time);break;
  case TRAJ_TEST1:setTrajTest1InitConditions(solution, time);break;
  case TRAJ_TEST2:setTrajTest2InitConditions(solution, time);break;
  case TRAJ_TEST3:setTrajTest3InitConditions(solution, time);break;
  case TRAJ_TEST4:setTrajTest4InitConditions(solution, time);break;
  case TRAJ_TEST5:setTrajTest5InitConditions(solution, time);break;
  case TRAJ_TEST6:setTrajTest6InitConditions(solution, time);break;
  case TRAJ_TEST7:setTrajTest7InitConditions(solution, time);break;
  case TRAJ_TEST8:setTrajTest8InitConditions(solution, time);break;
  case TRAJ_TEST9:setTrajTest9InitConditions(solution, time);break;
  default:
    FatalError("void CSimulation::setStandardEvaluationInitConditions(...): no handler for this type of trajectory");
    break;
  } //   switch (tt)
} // void CSimulation::setStandardEvaluationInitConditions(CIndividual& solution)



//------------------------------------------------------------------------------
void CSimulation::initStandardEvaluation(CIndividual& solution,
                                         bool final, bool monitor, double time)
{
  checkSensitivityFactorsInitialization();     // calc special factors if needed
  resetCurrentStepCounter();
  m_ga_performed = false;                          // no gravitiy assist is used

  setStandardEvaluationInitConditions(solution, time);

  m_nc->load(solution.getSolution());        // load solution to neurocontroller

  m_spacecraft.setSteeringStrategy(m_nc);               // set steering strategy

  // set initial launch conditions (flight direction, velocity)
  setLaunchState(m_t, m_initial_velocity,m_initial_azimuth,m_initial_elevation);
  m_spacecraft.ResetDeltaV();
  updateProximity();                                 // estimate first proximity
  calcBestValues();             // calculate remaining values of first best step
  setInitProximity(getCurrProximity());
  continueStandardEvaluationIntegrationLoop(false, false);
  performSteering();              // get first steering command from controller

  // store initial launch conditions in trajectory if final evaluation or 
  // monitor evaluation
//  if (monitor || final) addToSpacecraftTrajectory(m_spacecraft);
  
  initDESSolver(m_ys, m_t);        // transformation of CPolarState into Vector 
} // void CSimulation::initStandardEvaluation(...)

//-----------------------------------------------------------------------------
CEAFitness* CSimulation::evaluateStandardSolution(
        CIndividual& solution
 ,const bool final
 ,const bool monitor
 ,const double time
)
{
  const std::string method_name("CSimulation::evaluateStandardSolution(....)");
  // check for fitness object and set a new one if check fails
  // Ususally this is the case for freshly created individuals
  // which haven't been evaluated yet.
  if (solution.hasFitnessObject() != true)
  {
    try
    {
      solution.setFitnessObject(new CInTranceFitness());      
    }
    catch (std::bad_alloc)
    {
      Warning(method_name
        +std::string(": solution has no fitness object and allocation of a new\
                     one failed!\n"));
      return 0;
    }
  } // if (solution.hasFitnessObject() != true)

  // check for abort simulation condition
  if ( (shouldAbort()==true) && (monitor==false) && (final==false) )
  {
    setFitnessValues(solution.getFitnessObject(), HUGE_VAL);
    return 0;
  } // if ( (shouldAbort()==true) && (monitor==false) && (final==false) )

  // problem specific initialization stuff
  initStandardEvaluation(solution, final, monitor, time);
  integrationPreProcessing();          // problem specific pre-integration stuff
  if (NativeBridgeProbe::scheduling())
  {
    // Explicit-endpoint diagnostic only. Init above establishes the same native
    // initial state once. No optimiser-guessed shooting node is ever installed.
    NativeBridgeProbe::require(getIntegratorType()==_RKF54,
      "scheduled replay currently requires checked native RKF54");
    NativeBridgeProbe::require(NativeBridgeProbe::probe().imported.back().time-m_t<=getIntegrationInterval(),
      "scheduled duration exceeds native configured integration horizon");
    if (m_spacecraft_trajectory.size()<NativeBridgeProbe::probe().imported.size())
      m_spacecraft_trajectory.resize(NativeBridgeProbe::probe().imported.size());
    const auto checkPhysicalConstraints = [this]()
    {
      for (std::vector<constraints::CConstraint*>::iterator c=m_constraints.begin(); c!=m_constraints.end(); ++c)
      {
        // This propagation diagnostic has no target-acceptance event. Do not
        // terminate its supplied schedule at an incidental proximity match.
        if (*c==TRGT_PROX_CONSTRAINT.get()) continue;
        NativeBridgeProbe::require((*c)->isMatched(), "native physical path constraint violated during schedule");
      }
    };
    checkPhysicalConstraints();
    while (NativeBridgeProbe::has_next_interval())
    {
      if (monitor || final) addToSpacecraftTrajectory(m_spacecraft);
      const double h=NativeBridgeProbe::begin_interval(m_spacecraft);
      incCurrentStepCounter();
      setIntegrationStepSize(h);
      const RK_STATE result=integrate(m_ys,m_t,h);
      NativeBridgeProbe::finish_interval(m_spacecraft,static_cast<int>(result),result==RK_DONE);
      checkPhysicalConstraints();
      performSteering(); // import the next boundary, including terminal metadata
      updateProximity(); // preserve ordinary diagnostic/reporting bookkeeping
    }
    NativeBridgeProbe::record_endpoint(m_spacecraft);
    if (NativeBridgeProbe::mercury_profile())
    {
      // Evaluate native Mercury ephemeris at the attained endpoint, never at an
      // earlier best-proximity epoch. This diagnostic does not certify rendezvous.
      const Vector target_cart=m_solsys->State(getTargetBodyId(),m_spacecraft.MJD()).Cart();
      NativeBridgeProbe::record_target_endpoint(m_spacecraft,static_cast<int>(getTargetBodyId()),
        m_spacecraft.MJD(),target_cart);
    }
    // Legacy output must represent this explicit final state, not a previously
    // selected point. bridge_endpoint.csv is diagnostic state, not mission acceptance.
    calcBestValues();
    if (monitor || final) addToSpacecraftTrajectory(m_spacecraft);
    evaluationPostprocess();
    setFitnessValues(solution.getFitnessObject(),0.0);
    if (monitor || final) setFinalTargetState();
    return 0;
  }
  setIntegrationStepSize(m_h);                          // set initial step size
  RK_STATE itgr_result;
  // integration loop
  while (continueStandardEvaluationIntegrationLoop(final, monitor) == true)
  {
    if (monitor || final) addToSpacecraftTrajectory(m_spacecraft);   // store sc

    incCurrentStepCounter();                              // update step counter

    // dynamically update integration step size
    setIntegrationStepSize(calcIntegrationStepWidth(m_h, final));

    if (m_data.number_of_ga_bodies !=0) 
      performGravityAssist(m_ys, m_x, m_v, m_t);

    itgr_result = integrate(m_ys, m_t, getIntegrationStepSize());
    if ( (itgr_result != RK_DONE) && (itgr_result != RK_STOP_FUNCTION) )
    {
      setFitnessValues(solution.getFitnessObject(), HUGE_VAL);
      return 0;
    }
    if (getDynamicIntegrationIntervalReduction() == true
      && ( getFlightTime() >= m_time_max_dynamic)) break;

    performSteering();                           // obtain new steering commando
    m_spacecraft.CheckRemainingPropellant(m_h_dyn*DAY2SEC);

    updateProximity();          // update proximity object with current r&v data
    checkForNewBestStep();      // check for a possible new best step & store it

  } // end while

  if (monitor || final) addToSpacecraftTrajectory(getCurrentSC());   // store sc

  checkDynamicMaxIntegrationTimeReduction();

  if ( (final==true) && (getMinOutputPoints()>getMaxStepsCount()) )
  {
    evaluateInterpolated();
  } // if ( (final==true) && (getMinOutputPoints()>getMaxStepsCount()) )

  evaluationPostprocess();
  setFitnessValues(solution.getFitnessObject(), 0.0);     // set fitness object

  if ( (monitor==true) || (final == true) ) setFinalTargetState();
  return 0;
} // evaluateStandardSolution


//-----------------------------------------------------------------------------
void CSimulation::evaluationPostprocess(void)
{
  // the following code will be executed, if a special optimization goal is
  // used and one of the boundary conditions (relative distance) is matched
  if ( (getOptGoalType() == MAXIMUM_DEFLECTION) && (m_best_dr<=1.0) )
  {
    double t0(getStartTime()+m_time_best); // time of s/c hitting the asteroid
    double t2(m_impact_time);
    const Vector nea_s0(m_solsys->Cart(getTargetBodyId(), t0));
    const Vector sc_v0(m_best_s.Slice(3,5));

    m_best_deflection = calcDeflection(
      getStartTime()                         // s/c impact time (on asteroid)
      ,t2                    // original time of asteroid (NEA) impact on target
      ,nea_s0.Slice(0,2)         // asteroid position (t=t0)
      ,nea_s0.Slice(3,5)         // asteroid velocity (t=t0)
      ,sc_v0                     // s/c velocity (t=t0)
      ,m_impact_efficiency     // how much  momentum will be transferred to NEA
      );
  }
  else m_best_deflection = 0.0;

} // CSimulation::evaluationPostprocess(void)

//-----------------------------------------------------------------------------
inline void StoreSC(SCTraj& traj, const CSC& sc, const uint i)
{
  if ((i+10)>traj.size())
    traj.resize(std::max(size_t(traj.size()*1.25),size_t(10)));
  traj[i] = sc;
} //inline void StoreSC(SCTraj& traj, const CSC& sc, const uint i)


//-----------------------------------------------------------------------------
template<class T, class U>
CEAFitness* CSimulation::evaluateSteeringAngleApprox(
        CIndividual&  solution
 ,const U&            clock_angle
 ,const U&            cone_angle
 ,const T&            time
 ,      SCTraj&       traj
 ,      unsigned int& n_steps
 ,      double        start_t          // default value = 0.0
 ,bool  final
 ,bool  monitor
 )
{
#define NEW_CODE
  bool error = false;
  m_n_approx_sc_traj = 0;

  if ((clock_angle.size() == cone_angle.size()) &&
      (cone_angle.size() == time.size()) &&
      (clock_angle.size() > 0))
  {
    const double t_clock(clock_angle[clock_angle.size()-1]);
    const double t_cone(cone_angle[cone_angle.size()-1]);
    const double t_time(time[time.size()-1]);

    int steps_index = 0;
    int max_steps(time.size()-1);
    resetCurrentStepCounter();
    m_steps_best = 0;

    double h = getIntegrationInterval()/max_steps;

    m_ga_performed=false;

    RK4*    integrator_RK4(0);
    RKF54*  integrator_RKF54(0);
    DOPRI8* integrator_DOPRI8(0);
    DE*     integrator_DESG(0);

    if ( (init_sensitivity_factors == true) && (jt != J_AND) && (jt != J_OR) )
    {
      initSensitivityFactors();
      init_sensitivity_factors=false;
    }

    // initialize equations of motion
    initEOM(&m_spacecraft, integrator_RK4, integrator_RKF54, integrator_DOPRI8, integrator_DESG);

    double t(start_t); // set start time

    // Initialize s/c (primary to estimate propellant mass for s/c types
    // that are driven by mass consuming propulsion systems (not solar sail)
    if (m_spacecraft.isSolarSail())
    {
      m_spacecraft.Init(0.0);
    }
    else if (m_spacecraft.isNEP() || m_spacecraft.isSEP() )
    {
      const double mp=solution.getMappedAllele(solution.getLength()-1,0.0,1.0);
      m_spacecraft.Init(mp);
    }
    else
      FatalError("CSimulation::evaluateSteeringAngleApprox(): unknown spacecraft type");

    m_spacecraft.setLaunchState(t, getInitialBody(), m_last_saved_start_velocity
      , m_last_saved_start_azimuth, m_last_saved_start_elevation);
    m_spacecraft.ResetDeltaV();

    // iterator in time array
    unsigned int i_time(0);

    // get new steering angles
    unsigned int i_clock_angle(0), i_cone_angle(0);
    double clock_angle_val(clock_angle[i_clock_angle]);
    AssertAngularRange(clock_angle_val);
    double cone_angle_val(cone_angle[i_cone_angle]);

    // correct them to correct range
    MathLib::AssertAngularRange(clock_angle_val);

    // and use as steering input of s/c
    m_spacecraft.performSteering(clock_angle_val, cone_angle_val);
    m_spacecraft.UpdateAccDeltaVec(0.0);

    unsigned int i_traj(0);

    // set first item in new trajectory
#if defined NEW_CODE
    if (final == true || monitor == true)
#endif // NEW_CODE
     StoreSC(traj, m_spacecraft, 0);// traj[0]=m_spacecraft;

    double start_time = m_spacecraft.MJD();

    // transformation of CPolarState into Vector
    Vector ys(8), x(3), v(3);

    initDESSolver(
       m_spacecraft
      ,getIntegratorType()
      ,ys
      ,t
      ,integrator_RK4
      ,integrator_RKF54
      ,integrator_DOPRI8
      ,integrator_DESG
    );

    updateProximity();
    m_best_proximity = getCurrProximity();
    calcBestValues();

    double h_dyn = h;
    const unsigned int i_cone_angle_end(cone_angle.size()-1);

    // integration
    while ( (i_cone_angle < i_cone_angle_end)
     && !getCurrProximity().CheckProximity(m_max_dr_AU, m_max_dv_AUd)
     && m_spacecraft.CheckSolarDistance(getMinSolarDistance())
     && m_spacecraft.CheckTemperature() )
    {
      h_dyn = time[i_time+1]-time[i_time];

      //if (n_ga_bodies != 0) performGravityAssist(ys, x, v, t);

      // integration of equations of motion
      if (integrate(m_spacecraft, ys, t, h_dyn ,getIntegratorType()
        ,integrator_RK4,    integrator_RKF54
        ,integrator_DOPRI8, integrator_DESG
        ) != RK_DONE)
      {
        SafeDelete(integrator_RK4);SafeDelete(integrator_RKF54);
        SafeDelete(integrator_DOPRI8);SafeDelete(integrator_DESG);
        setFitnessValues(&m_approx_fitness,-HUGE_VAL);
        return 0;
      }

      // go one step further in all arrays
      ++i_traj;
      ++i_clock_angle;
      ++i_cone_angle;
      ++steps_index;
      incCurrentStepCounter();

      // get new steering angles
      clock_angle_val = clock_angle[i_clock_angle];
      AssertAngularRange(clock_angle_val);
      cone_angle_val  = cone_angle[i_cone_angle];

      MathLib::AssertAngularRange(clock_angle_val);              // possibly correct angles

      // perform steering
      m_spacecraft.performSteering(clock_angle_val, cone_angle_val);

      // propellant left check
      m_spacecraft.CheckRemainingPropellant(h_dyn*DAY2SEC);

      // update of acceleration vector
      m_spacecraft.UpdateAccDeltaVec(h_dyn);

      // add s/c to approximation trajectory
#if defined NEW_CODE
      if (final == true || monitor == true)
#endif // NEW_CODE
        StoreSC(traj, m_spacecraft, i_traj);

      updateProximity(); // update proximity object

      // check proximity
      if (getCurrProximity().GetNormedProximity(m_max_dr_AU, m_max_dv_AUd)
        < m_best_proximity.GetNormedProximity(m_max_dr_AU, m_max_dv_AUd)
        || tt == RIMPACT
        || tt == ACTPROBLEM)
      {
        m_steps_best=getCurrentStepCount();
        m_time_best=t-
#if defined NEW_CODE
start_time;
#else //NEW_CODE
traj[0].MJD();
#endif // NEW_CODE
        m_best_proximity=getCurrProximity();
        calcBestValues();
      }

   } // end while

   SafeDelete(integrator_RK4);SafeDelete(integrator_RKF54);
   SafeDelete(integrator_DOPRI8);SafeDelete(integrator_DESG);

   n_steps = m_steps_best+1;

  // the following code will be executed, if a special optimization goal is
  // used and one of the boundary conditions (relative distance) is matched
  if ( (getOptGoalType() == MAXIMUM_DEFLECTION) && (m_best_dr<=1.0) )
  {

    double t0(
#if defined NEW_CODE
      start_time                          // time of s/c hitting target asteroid
#else
      traj[0].MJD()
#endif
      +m_time_best);
    double t2(m_impact_time);
    const Vector nea_s0(m_solsys->Cart(getTargetBodyId(), t0));
    const Vector sc_v0(m_spacecraft_trajectory[m_steps_best].CState::Cart().Slice(3,5));

    m_best_deflection = calcDeflection(
      t0                         // s/c impact time (on asteroid)
      ,t2                        // original time of asteroid (NEA) impact on target
      ,nea_s0.Slice(0,2)         // asteroid position (t=t0)
      ,nea_s0.Slice(3,5)         // asteroid velocity (t=t0)
      ,sc_v0                     // s/c velocity (t=t0)
      ,m_impact_efficiency       // how much linear momentum will be transferred to NEA
      );
  }
  else m_best_deflection = 0.0;

   // set normal fitness function
   setFitnessValues(&m_approx_fitness,0.0);

  } // if (m_clock_angle_history_approximation && m_cone_angle_history_approximation)
  return 0;
#undef NEW_CODE
}; // template<class T> void evaluateSteeringAngleApprox(...)


// -----------------------------------------------------------------------------
void CSimulation::evaluate(const Vector&lsl_vector)
{
  initEOM(&m_spacecraft);
  double h = getIntegrationInterval()/getMaxStepsCount();

  if ( getMinStartTime() != getMaxStartTime() )
    FatalError("CSimulation::Evaluate(): minimum and maximum start time are not equal");

  double t = getMinStartTime();

  if (m_spacecraft.isSolarSail())
    m_spacecraft.Init(0.0);
  else if (m_spacecraft.isNEP()||m_spacecraft.isSEP())
    m_spacecraft.Init();
  else
    FatalError("CSimulation::Evaluate(): unknown spacecraft type");

  // set spacecraft steering strategy

  m_spacecraft.setSteeringStrategy(lsl_vector);

  // set initial spacecraft values

  if (getMinLaunchVelocity() != getMaxLaunchVelocity())
    FatalError("CSimulation::Evaluate(): minimum and maximum launch velocity are not equal");
  if (getMinLaunchAzimuthAngle()!=getMaxLaunchAzimuthAngle())
    FatalError("CSimulation::Evaluate(): minimum and maximum launch azimuth are not equal");
  if (getMinLaunchElevationAngle()!=getMaxLaunchElevationAngle())
    FatalError("CSimulation::Evaluate(): minimum and maximum launch elevation are not equal");
  setLaunchState(t, getMinLaunchVelocity(), getMinLaunchAzimuthAngle(),
   getMinLaunchElevationAngle());

  m_spacecraft.ResetDeltaV();

  m_spacecraft.performSteering(lsl_vector, m_solsys->State(getTargetBodyId(), t));

  m_spacecraft.UpdateAccDeltaVec(0.0);

  // spacecraft trajectory
  resetCurrentStepCounter();
  addToSpacecraftTrajectory(m_spacecraft);

  // transformation of CPolarState into Vector

  Vector ys(8), x(3), v(3);

  initDESSolver(ys, t);

  // integration
  while (getCurrentStepCount()<getMaxStepsCount()&&
   m_spacecraft.CheckSolarDistance(getMinSolarDistance())&&
   m_spacecraft.CheckTemperature()) {
     incCurrentStepCounter();

     integrate(ys, t, h);

     m_spacecraft.performSteering(lsl_vector,
      m_solsys->State(getTargetBodyId(), t));
     m_spacecraft.CheckRemainingPropellant(h*DAY2SEC);
     m_spacecraft.UpdateAccDeltaVec(h);

     addToSpacecraftTrajectory(m_spacecraft);
  }

  m_steps_best=getCurrentStepCount();

  freeEOM();
};  // CSimulation::evaluate(const Vector&lsl_vector)


// -----------------------------------------------------------------------------
void CSimulation::evaluate(
  const CSCAttitudeHistory& att_history)
{

  resetCurrentStepCounter();
  initEOM(&m_spacecraft);
  setMaxStepsCount(att_history.getDimension());
  m_spacecraft_trajectory.resize(getMaxStepsCount()+1);

  double h=att_history.getMJD(getCurrentStepCount()+1)-att_history.getMJD(getCurrentStepCount());

  double t=att_history.getMJD(getCurrentStepCount());

  if (m_spacecraft.isSolarSail())
    m_spacecraft.Init(t);
  else if (m_spacecraft.isNEP()||m_spacecraft.isSEP())
    m_spacecraft.Init();
  else
    FatalError("CSimulation::Evaluate(): unknown spacecraft type");

  // set initial spacecraft values

  if (getMinLaunchVelocity() != getMaxLaunchVelocity())
    FatalError("CSimulation::Evaluate(): minimum and maximum launch velocity are not equal");
  if (getMinLaunchAzimuthAngle()!=getMaxLaunchAzimuthAngle())
    FatalError("CSimulation::Evaluate(): minimum and maximum launch azimuth are not equal");
  if (getMinLaunchElevationAngle()!=getMaxLaunchElevationAngle())
    FatalError("CSimulation::Evaluate(): minimum and maximum launch elevation are not equal");
  setLaunchState(t, getMinLaunchVelocity(), getMinLaunchAzimuthAngle(),
   getMinLaunchElevationAngle());

  m_spacecraft.ResetDeltaV();

  double act_clock = att_history.getClockAngle(getCurrentStepCount());
  double act_cone  = att_history.getConeAngle(getCurrentStepCount());
  m_spacecraft.performSteering(act_clock, act_cone);
  double next_clock(0.0), next_cone(0.0);

  m_spacecraft.UpdateAccDeltaVec(0.0);

  // spacecraft trajectory
  addToSpacecraftTrajectory(m_spacecraft);

  // transformation of CPolarState into Vector
  Vector ys(8), x(3), v(3);

  initDESSolver(ys, t);

  // integration
  int file_steps=getCurrentStepCount();

  double t_end=att_history.getMJD(att_history.getDimension()-1);

  double some_small_error=1e-6; // due to problems with GESOP

  while (t<t_end-some_small_error&&
   m_spacecraft.CheckSolarDistance(getMinSolarDistance())&&
   m_spacecraft.CheckTemperature())
  {
    double act_mjd  =att_history.getMJD(file_steps);
    act_clock=att_history.getClockAngle(file_steps);
    act_cone =att_history.getConeAngle(file_steps);

    file_steps++;

    double next_mjd=att_history.getMJD(file_steps);
    if (next_mjd==act_mjd&&next_mjd!=t_end)
      next_mjd=att_history.getMJD(++file_steps);

    next_clock =att_history.getClockAngle(file_steps);
    next_cone= att_history.getConeAngle(file_steps);

    incCurrentStepCounter();

    h=next_mjd-act_mjd;

    if (h==0.0) break;

    integrate(ys, t, h);

    double clock_set=att_history.getClockAngle(file_steps);
    double cone_set= std::min(att_history.getConeAngle(file_steps), HALFPI);
    m_spacecraft.performSteering(clock_set, cone_set);
    m_spacecraft.CheckRemainingPropellant(h*DAY2SEC);
    m_spacecraft.UpdateAccDeltaVec(h);

    addToSpacecraftTrajectory(m_spacecraft);
  }

  m_steps_best=getCurrentStepCount();

  freeEOM();

}; // end CSimulation::Evaluate

// -----------------------------------------------------------------------------
void CSimulation::writeGTOC2Csv(const std::string& file_name)
{
  using namespace std;

  if (!file_name.empty())
  {
    std::string s(file_name);    
    const unsigned int p = s.find('.');
    if (p!=string::npos)
    {
      s.replace(s.begin()+p,s.end(),".gtoc2");
    } else s.append(".gtoc2");

    ofstream os(s.c_str());
    Assure(os, s.c_str());

    unsigned int idx = 0;

    string old_name = "";
    while (m_gtoc2_best_bodies[idx]!=0)
    {
      const BodyID body_id = m_gtoc2_best_bodies[idx];
      string _name = m_solsys->getBodyName(m_gtoc2_best_bodies[idx]);
      string name(_name.begin()+1,_name.end());
      // data output
      os << "# Heliocentric trajectory of spacecraft\n";
      os << "#\n";
      os << "#PHASE NUMBER: " << noshowpos << dec << (idx+1) << endl;
      os << "#DESCRIPTION: From";
      switch (idx)
      {
      case 0: os << " the Earth";break;
      default:
        os << " Asteroid " << old_name;
        break;
      } // switch (idx)
      os << " to Asteroid " << name << endl;
      os << "#";
      const unsigned int w = 22;
      os << setw(w) << "Time (MJD)" << setw(0);
      os << setw(w) << "x (km)" << setw(0);
      os << setw(w) << "y (km)" << setw(0);
      os << setw(w) << "z (km)" << setw(0);
      os << setw(w) << "vx (km/s)" << setw(0);
      os << setw(w) << "vy (km/s)" << setw(0);
      os << setw(w) << "vz (km/s)" << setw(0);
      os << setw(w) << "mass (kg)" << setw(0);
      os << setw(w) << "Thrust_x (N)" << setw(0);
      os << setw(w) << "Thrust_y (N)" << setw(0);
      os << setw(w) << "Thrust_z (N)" << setw(0);
      os << endl;

      int j = 0;
      os << setprecision(13);os << scientific;
      for (j=0;j<m_steps_best;++j)
      {
        os << setw(0) << " ";
        os << setw(w) << (m_time_best-m_spacecraft_trajectory[j].MJD());
        os << setw(w) << m_spacecraft_trajectory[j].CState::x()*AU2KM;
        os << setw(w) << m_spacecraft_trajectory[j].CState::y()*AU2KM;
        os << setw(w) << m_spacecraft_trajectory[j].CState::z()*AU2KM;
        os << setw(w) << m_spacecraft_trajectory[j].CState::dx()*AUD2KMS;
        os << setw(w) << m_spacecraft_trajectory[j].CState::dy()*AUD2KMS;
        os << setw(w) << m_spacecraft_trajectory[j].CState::dz()*AUD2KMS;
        os << setw(w) << m_spacecraft_trajectory[j].getMass();
        const MathLib::Vector t(m_spacecraft_trajectory[j].ThrustAcc_xyz()*
          m_spacecraft_trajectory[j].getMass());
        for (int k=0;k<3;++k) os << setw(w) << t[k];
        os << endl;
      } // for (j=0;j<m_steps_best;++j)

      os << endl;
      old_name = name;
      ++idx;
    } // while (m_gtoc2_best_bodies[idx]!=0)
    
  } // if (!file_name.empty())
} // void CSimulation::writeGTOC2Csv(const std::string file_name)

// -----------------------------------------------------------------------------
void CSimulation::writeACTCsv(const std::string& file_name)
{

  using namespace std;

  if (!file_name.empty())
  {
    std::string s(file_name);    
    const unsigned int p = s.find('.');
    if (p!=string::npos)
    {
      s.replace(s.begin()+p,s.end(),".act");
    } else s.append(".gtoc2");

    ofstream act(s.c_str());
    Assure(act, s.c_str());


    const int fw(20),fp(14);
    Vector acc_xyz;

    act << "# PHASE NUMBER: 1" << endl;
    act << "# DESCRIPTION: From the Earth to " << m_solsys->Name(getTargetBodyId()) << endl;
    act << uppercase << scientific << showpoint << showpos << left << setw(22) << setprecision(fp) << "#    Time (MJD2000)" << " ";
    act << uppercase << scientific << showpoint << showpos << left << setw(22) << setprecision(fp) << "        x (km)" << " ";
    act << uppercase << scientific << showpoint << showpos << left << setw(22) << setprecision(fp) << "        y (km)" << " ";
    act << uppercase << scientific << showpoint << showpos << left << setw(22) << setprecision(fp) << "        z (km)" << " ";
    act << uppercase << scientific << showpoint << showpos << left << setw(22) << setprecision(fp) << "       vx (km/s)" << " ";
    act << uppercase << scientific << showpoint << showpos << left << setw(22) << setprecision(fp) << "       vy (km/s)" << " ";
    act << uppercase << scientific << showpoint << showpos << left << setw(22) << setprecision(fp) << "       vz (km/s)" << " ";
    act << uppercase << scientific << showpoint << showpos << left << setw(22) << setprecision(fp) << "       mass (kg)" << " ";
    act << uppercase << scientific << showpoint << showpos << left << setw(22) << setprecision(fp) << "         ux (N)" << " ";
    act << uppercase << scientific << showpoint << showpos << left << setw(22) << setprecision(fp) << "         uy (N)" << " ";
    act << uppercase << scientific << showpoint << showpos << left << setw(22) << setprecision(fp) << "         uz (N)" << " ";
    act << endl;

    for (int i=0; i<=m_steps_best; ++i)
    {
      act << uppercase << scientific << showpoint << showpos << right << setw(fw) << setprecision(fp) << (m_spacecraft_trajectory[i].MJD()-MJD_J2000) << " ";
      act << uppercase << scientific << showpoint << showpos << right << setw(fw) << setprecision(fp) << (m_spacecraft_trajectory[i].CState::x()*AU2KM) << " ";
      act << uppercase << scientific << showpoint << showpos << right << setw(fw) << setprecision(fp) << (m_spacecraft_trajectory[i].CState::y()*AU2KM) << " ";
      act << uppercase << scientific << showpoint << showpos << right << setw(fw) << setprecision(fp) << (m_spacecraft_trajectory[i].CState::z()*AU2KM) << " ";
      act << uppercase << scientific << showpoint << showpos << right << setw(fw) << setprecision(fp) << (m_spacecraft_trajectory[i].CState::dx()*AUD2KMS) << " ";
      act << uppercase << scientific << showpoint << showpos << right << setw(fw) << setprecision(fp) << (m_spacecraft_trajectory[i].CState::dy()*AUD2KMS) << " ";
      act << uppercase << scientific << showpoint << showpos << right << setw(fw) << setprecision(fp) << (m_spacecraft_trajectory[i].CState::dz()*AUD2KMS) << " ";
      act << uppercase << scientific << showpoint << showpos << right << setw(fw) << setprecision(fp) <<  m_spacecraft_trajectory[i].getMass() << " ";
      acc_xyz.Set(m_spacecraft_trajectory[i].ThrustAcc_xyz()*m_spacecraft_trajectory[i].getMass());
      act << uppercase << scientific << showpoint << showpos << right << setw(fw) << setprecision(fp) <<  acc_xyz(0) << " ";
      act << uppercase << scientific << showpoint << showpos << right << setw(fw) << setprecision(fp) <<  acc_xyz(1) << " ";
      act << uppercase << scientific << showpoint << showpos << right << setw(fw) << setprecision(fp) <<  acc_xyz(2) << " ";
      act << endl;
    } // for (int i=0; i<=m_steps_best; ++i)
  } // if (!file_name.empty())
} // void writeACTCsv(const std::string filename)


// ----------------------------------------------------------------------------
void CSimulation::PrintCSV(const char* filename)
{

  ofstream csv(filename);
  Assure(csv, filename);

  csv<<"TOF [d], ";

  // spacecraft
  m_spacecraft.PrintHeader(csv); csv<<", ";

  // initial body
  if (getInitialStateType()==BODY_STATE)
  {
    m_solsys->PrintHeader(
      getInitialBody(), csv, m_solsys->Name(getInitialBody()));
    csv<<", ";
    m_spacecraft.PrintDeltaHeader(csv, m_solsys->Name(getInitialBody()));
    if (getTargetStateType()!=FREE_ORBIT_STATE) csv<<", ";
    // the line above is a workaround and must be changed
  }

  // target body
  if (getTargetStateType()==BODY_STATE||getTargetStateType()==BODY_ORBIT_STATE) {
    m_solsys->PrintHeader(getTargetBodyId(), csv, m_solsys->Name(getTargetBodyId()));
    csv<<", ";
    m_spacecraft.PrintDeltaHeader(csv, m_solsys->Name(getTargetBodyId()));
    csv<<", ";
    m_spacecraft.PrintProximityHeader(csv, m_solsys->Name(getTargetBodyId()));
  }

  // NC output, if applicable
  if (m_nc!=NULL) {
    m_spacecraft.PrintNCInputHeader(csv);
    m_spacecraft.PrintNCOutputHeader(csv);
  }

  csv<<endl;

  for (int i=0; i<=m_steps_best; ++i) 
  {
    csv<<m_spacecraft_trajectory[i].MJD()-
     m_spacecraft_trajectory[0].MJD()<<", ";

    // spacecraft
    m_spacecraft_trajectory[i].Print(csv); csv<<", ";

    // initial body
    if (getInitialStateType()==BODY_STATE) 
    {
      m_solsys->Print(getInitialBody(), m_spacecraft_trajectory[i].MJD(), csv);
      csv<<", ";
      m_spacecraft_trajectory[i].PrintDelta(csv, getInitialBody(), true);
      if (getTargetStateType()!=FREE_ORBIT_STATE) csv<<", ";
      // the line above is a workaround and must be changed
    }

    // target body
    if ((getTargetStateType() == BODY_STATE) || 
      (getTargetStateType()==BODY_ORBIT_STATE)) 
    {
      m_solsys->Print(getTargetBodyId(), m_spacecraft_trajectory[i].MJD(), csv);
      csv<<", ";
      m_spacecraft_trajectory[i].PrintDelta(csv, getTargetBodyId()); csv<<", ";
      m_spacecraft_trajectory[i].PrintProximity(csv, getTargetBodyId());
    }

    // NC output, if applicable
    if (m_nc!=NULL)
    {
      m_spacecraft_trajectory[i].PrintNCInput(csv);
      m_spacecraft_trajectory[i].PrintNCOutput(csv);
    }

    csv<<endl;
  }
}; // end CSimulation::PrintCSV


// -----------------------------------------------------------------------------
void CSimulation::PrintDAT(const char* filename)
{
  ofstream dat(filename);
  Assure(dat, filename);

  const double start_time = m_spacecraft_trajectory[0].MJD();
  const double final_time = m_spacecraft_trajectory[m_steps_best].MJD();
  const double total_flight_time = (final_time - start_time);

  const bool draw_start_body_orbit = 
    (m_solsys->e(getInitialBody(), start_time)<1.0) && (total_flight_time < m_solsys->OrbitPeriod(getInitialBody(), start_time));

  const bool draw_target_body_orbit = 
    (m_solsys->e(getTargetBodyId(), start_time)<1.0) && (total_flight_time < m_solsys->OrbitPeriod(getTargetBodyId(), start_time));

  // print header

  dat<<"TITLE = \"InTrance 1.7 Plot\"" <<endl;

  dat << 
    "VARIABLES = \"TOF\", \"MJD\", \"x\", \"y\", \"z\", \"ax\", \"ay\", \"az\"";

  dat<<", \"alpha-+180\", \"alpha0-360\", \"beta\"";

  if (getCurrentSC().isSolarSail() == true)
  {
    if (getCurrentSC().isRealSail() == true) dat<<", \"T\"";
  } else // no solar sail
  {
    dat<<", \"throttle\"";
    if (getCurrentSC().isRIT())
    {
      dat << ", \"n\"";
    } // if (getCurrentSC().isRIT())
  } // else branch of if (getCurrentSC().isSolarSail() == true)

  if (getInitialStateType()==BODY_STATE)
  {
    dat<<", \"xI\", \"yI\", \"zI\"";
    if (draw_start_body_orbit == true)
    {
      dat << ", \"xI2\", \"yI2\", \"zI2\"";
    } // if (draw_target_body_orbit == true)
  } // if (getInitialStateType()==BODY_STATE)

  switch (getTargetStateType())
  {
  case BODY_STATE:
  case BODY_ORBIT_STATE:
    dat << ", \"xT\", \"yT\", \"zT\"";
    if (draw_target_body_orbit == true)
    {
      dat << ", \"xT2\", \"yT2\", \"zT2\"";
    } // if (draw_target_body_orbit == true)
    break;
  default:break;
  } // switch (getTargetStateType())

  dat << endl; // end of header

  // values
  Vector v_temp(3),s_temp(6);
  int j;
  for (int i=0; i<=m_steps_best; ++i)
  {
    const CSC* sc = &(m_spacecraft_trajectory[i]);
    const double time = sc->MJD();
    const double flight_time = (time - start_time);

    // flight time & date
    dat << flight_time << '\t' << time << '\t';

    // spacecraft position
    dat<<sc->CState::x()<<'\t'<< sc->CState::y()<<'\t'<< sc->CState::z()<<'\t';

    /// cartesian thrust acceleration of s/c
    v_temp = sc->ThrustAcc_xyz();
    for (j=0;j<3;++j) dat << v_temp(j) << '\t';

    // steering angle 1 (clock angle)
    double alpha1=m_spacecraft_trajectory[i].getClockAngle();
    dat<<alpha1<<'\t';

    double alpha2=(alpha1>0.0)?alpha1:alpha1+TWOPI;
    dat<<alpha2<<'\t';

    // steering angle 2 (cone angle)
    dat<<m_spacecraft_trajectory[i].getConeAngle()<<'\t';

    if (sc->isSolarSail())
    {
      if (sc->isRealSail()) dat << sc->getTemperature()<<'\t';
    } // if (sc->isSolarSail())
    else
    {
      dat << sc->getThrottle()<<'\t';
      if (sc->isRIT())
      {
        dat << sc->getNumberOfUsedThrusters() << '\t';
      } // if (sc->isRIT())
    } // else branch of "if (sc->isSolarSail())"

    // initial body
    if (getInitialStateType()==BODY_STATE)
    {
      const BodyID body = getInitialBody();
      dat<<m_solsys->x(body, time)<<'\t';
      dat<<m_solsys->y(body, time)<<'\t';
      dat<<m_solsys->z(body, time)<<'\t';
      if (draw_start_body_orbit == true)
      {
        const double date = 
          (flight_time/total_flight_time)*
          m_solsys->OrbitPeriod(body, start_time)+start_time;
        dat << m_solsys->x(body, date) << '\t';
        dat << m_solsys->y(body, date) << '\t';
        dat << m_solsys->z(body, date) << '\t';
      } // if (draw_target_body_orbit == true)
    } // if (getInitialStateType()==BODY_STATE)

    // target body
    if ( (getTargetStateType() == BODY_STATE) || 
      (getTargetStateType() == BODY_ORBIT_STATE))
    {
      const BodyID target = getTargetBodyId();
      dat << m_solsys->x(target, time)<<'\t';
      dat << m_solsys->y(target, time)<<'\t';
      dat << m_solsys->z(target, time)<<'\t';
      if (draw_target_body_orbit == true)
      {
        const double date = 
          (flight_time/total_flight_time)*
          m_solsys->OrbitPeriod(target, start_time)+start_time;
        dat << m_solsys->x(target, date) << '\t';
        dat << m_solsys->y(target, date) << '\t';
        dat << m_solsys->z(target, date) << '\t';
      } // if (draw_target_body_orbit == true)
    } // if ( (getTargetStateType() == BODY_STATE) || ...
    dat << endl;
  } // for (int i=0; i<=m_steps_best; ++i)
} // end CSimulation::PrintDAT


// ----------------------------------------------------------------------------
void CSimulation::PrintGESOP(const char*filename) {

  ofstream dat(filename);
  Assure(dat, filename);

  // print header

  dat<<"Time" <<'\t';
  dat<<"r" <<'\t';
  dat<<"phi" <<'\t';
  dat<<"theta" <<'\t';
  dat<<"dr" <<'\t';
  dat<<"dphi" <<'\t';
  dat<<"dtheta" <<'\t';
  dat<<"Clock Angle" <<'\t';
  if (m_spacecraft_trajectory[0].isSolarSail())
  {
    dat<<"Cone Angle" << endl;
  }
  else
  {
    dat<<"Cone Angle" <<'\t';
    dat<<"Throttle" <<endl;
  }

  dat<<"t [days]" <<'\t';
  dat<<"r [AU]" <<'\t';
  dat<<"phi [rad]" <<'\t';
  dat<<"theta [rad]" <<'\t';
  dat<<"dr [AU/day]" <<'\t';
  dat<<"dphi [rad/day]" <<'\t';
  dat<<"dtheta [rad/day]" <<'\t';
  dat<<"alpha [rad]" <<'\t';
  if (m_spacecraft_trajectory[0].isSolarSail())
  {
    dat<<"beta [rad]" << endl;
  }
  else
  {
    dat<<"beta [rad]" <<'\t';
    dat<<"throttle" <<endl;
  }

  for (int i=0; i<=m_steps_best; ++i) // print values
  {
    dat<<m_spacecraft_trajectory[i].MJD()<<'\t';
    dat<<m_spacecraft_trajectory[i].CState::r()<<'\t';
    dat<<m_spacecraft_trajectory[i].CState::phi()<<'\t';
    dat<<m_spacecraft_trajectory[i].CState::theta()<<'\t';
    dat<<m_spacecraft_trajectory[i].CState::dr()<<'\t';
    dat<<m_spacecraft_trajectory[i].CState::dphi()<<'\t';
    dat<<m_spacecraft_trajectory[i].CState::dtheta()<<'\t';
    double clockangle = m_spacecraft_trajectory[i].getClockAngle();
    MathLib::AssertAngularRange2pi(clockangle);
    dat<<clockangle<<'\t';
    if (m_spacecraft_trajectory[0].isSolarSail())
    {
      dat<<m_spacecraft_trajectory[i].getConeAngle()<<endl;
    }
    else
    {
      dat<<m_spacecraft_trajectory[i].getConeAngle()<<'\t';
      dat<<m_spacecraft_trajectory[i].getThrottle()<<endl;
    }

  }

}; // end CSimulation::PrintGESOP


// ----------------------------------------------------------------------------
// quick and dirty
void MultiplyQuaternion_(const double q1[4], const double q2[4]
, double result[4])
{
  result[0] = q2[3] * q1[0] + q2[0] * q1[3] + q2[1] * q1[2] - q2[2] * q1[1];
  result[1] = q2[3] * q1[1] + q2[1] * q1[3] + q2[2] * q1[0] - q2[0] * q1[2];
  result[2] = q2[3] * q1[2] + q2[2] * q1[3] + q2[0] * q1[1] - q2[1] * q1[0];
  result[3] = q2[3] * q1[3] - q2[0] * q1[0] - q2[1] * q1[1] - q2[2] * q1[2];
}


// ----------------------------------------------------------------------------
void CSimulation::PrintVRML(const char* filename)
{
  /*! \todo [AO]: implementation of exceptions in case VRML output file is not
  accessible*/
  ofstream wrl(filename);
  Assure(wrl, filename);

  int pind, pstop, i;
  double factor = 1.0;
  const double t_start = m_spacecraft_trajectory[0].MJD();
  const double t_end =  m_spacecraft_trajectory[m_steps_best].MJD();
  const double dt_flight_days = (t_end-t_start);
  const double dt_flight_years = dt_flight_days*DAY2JYR;

  // header
  wrl << "#VRML V2.0 utf8\n";

  // navigation info
  wrl << "\nNavigationInfo\n";
  wrl << "{\n";
  wrl << " type \"EXAMINE\"\n";
  wrl << "}\n";

  // light source
  wrl << "\nPointLight\n";
  wrl << "{\n";
  wrl << " ambientIntensity 1\n";
  wrl << " intensity 1\n";
  wrl << " location 0 0 0\n";
  wrl << "}\n";

  // viewpoints
  wrl << "\nViewpoint\n";
  wrl << "{\n";
  wrl << " position 0.0 0.0 3.0\n";
  wrl << " orientation 0.0 0.0 1.0 0.0\n";
  wrl << " description \"polar view\"\n";
  wrl << "}\n";
  wrl << "\nViewpoint\n";
  wrl << "{\n";
  wrl << " position 0.0 -3.0 0.0\n";
  wrl << " orientation 1.0 0.0 0.0 1.57\n";
  wrl << " description \"ecliptic view\"\n";
  wrl << "}\n";

  // clock
  wrl << "\nDEF time TimeSensor\n";
  wrl << "{\n";
  // every year of flight time is represented by one minute of visualization
  wrl << " cycleInterval " << (dt_flight_years * 60.0) << endl;
  wrl << " loop TRUE\n";
  wrl << " startTime -1\n";
  wrl << "}\n";
  wrl << "\nDEF StartSensor VisibilitySensor {}\n";
  wrl << "ROUTE StartSensor.enterTime TO time.set_startTime\n";

  // prototype of planet's appearance script code
  wrl << "\nPROTO Planet [\n";
  wrl << "field SFColor pl_emissiveColor 0.0 0.0 0.0\n";
  wrl << "field SFColor pl_diffuseColor 0.0 0.0 0.0\n";
  wrl << "field SFFloat pl_radius 0.001\n";
  wrl << "field SFFloat pl_transparency 0.4]\n";
  wrl << "{\n";
  wrl << " Shape\n";
  wrl << " {\n";
  wrl << "  appearance Appearance\n";
  wrl << "  {\n";
  wrl << "   material Material\n";
  wrl << "   {\n";
  wrl << "    emissiveColor IS pl_emissiveColor\n";
  wrl << "    diffuseColor IS pl_diffuseColor\n";
  wrl << "    transparency IS pl_transparency\n";
  wrl << "   }\n";
  wrl << "  }\n";
  wrl << "  geometry Sphere\n";
  wrl << "  {\n";
  wrl << "   radius IS pl_radius\n";
  wrl << "  }\n";
  wrl << " }\n";
  wrl << "}\n";

  wrl << "\nDEF Sun Planet\n";
  wrl << "{\n";
  wrl << " pl_emissiveColor 1.0 0.9 0.0\n";
  wrl << " pl_radius 0.00465\n";
  wrl << "}\n";

  wrl << "\nDEF SunSphere Planet\n";
  wrl << "{\n";
  wrl << " pl_diffuseColor 1.0 0.8 0.0\n";
  wrl << " pl_radius " << getMinSolarDistance() << endl;
  wrl << " pl_transparency 0.8\n";
  wrl << "}\n";

  // axes
  wrl << "\nDEF Axes Shape\n";
  wrl << "{" << endl;
  wrl << " appearance Appearance\n";
  wrl << " {\n";
  wrl << "  material Material\n";
  wrl << "  {" << endl;
  wrl << "   emissiveColor 1.0 1.0 1.0\n";
  wrl << "  }\n";
  wrl << " } # appearance\n";
  wrl << " geometry IndexedLineSet\n";
  wrl << " {\n";
  wrl << "  coord Coordinate\n";
  wrl << "  {\n";
  wrl << "   point\n";
  wrl << "   [\n";
  wrl << "    -1.0  0.0  0.0," << endl;
  wrl << "     1.0  0.0  0.0," << endl;
  wrl << "     0.0 -1.0  0.0," << endl;
  wrl << "     0.0  1.0  0.0," << endl;
  wrl << "     0.0  0.0 -1.0," << endl;
  wrl << "     0.0  0.0  1.0\n";
  wrl << "   ]\n";
  wrl << "  } # coord\n";
  wrl << "  coordIndex\n";
  wrl <<"   [\n";
  wrl << "   0, 1, -1" << endl;
  wrl << "   2, 3, -1" << endl;
  wrl << "   4, 5\n";
  wrl << "  ]\n";
  wrl << " } # geometry\n";
  wrl << "} # DEF AXES\n\n";

  wrl << "Transform\n";
  wrl << "{\n";
  wrl << " children\n";
  wrl << " [\n";
  wrl << "  Shape\n";
  wrl << "  {\n";
  wrl << "   appearance Appearance\n";
  wrl << "   {\n";
  wrl << "    material Material\n";
  wrl << "    {\n";
  wrl << "     emissiveColor 1.0 1.0 1.0\n";
  wrl << "    } # material\n";
  wrl << "   } # appearance" << endl;
  wrl << "   geometry Text\n";
  wrl << "   {\n";
  wrl << "    string \"x\"\n";
  wrl << "    fontStyle FontStyle\n";
  wrl << "    {\n";
  wrl << "     family \"SANS\"\n";
  wrl << "     size 0.2\n";
  wrl << "    } # fontStyle\n";
  wrl << "   } # geometry\n";
  wrl << "  } # shape\n";
  wrl << " ] # children\n";
  wrl << " translation 0.9 -0.15 0.0\n";
  wrl << "} # Transform\n";

  wrl << "Transform\n";
  wrl << "{\n";
  wrl << " children\n";
  wrl << " [\n";
  wrl << "  Shape\n";
  wrl << "  {\n";
  wrl << "   appearance Appearance\n";
  wrl << "   {\n";
  wrl << "    material Material\n";
  wrl << "    {\n";
  wrl << "     emissiveColor 1.0 1.0 1.0\n";
  wrl << "    } # material\n";
  wrl << "   } # appearance" << endl;
  wrl << "   geometry Text\n";
  wrl << "   {\n";
  wrl << "    string \"y\"\n";
  wrl << "    fontStyle FontStyle\n";
  wrl << "    {\n";
  wrl << "     family \"SANS\"\n";
  wrl << "     size 0.2\n";
  wrl << "    } # fontStyle\n";
  wrl << "   } # geometry\n";
  wrl << "  } # shape\n";
  wrl << " ] # children\n";
  wrl << " translation -0.15 0.9 0.0\n";
  wrl << "} # Transform\n";

  wrl << "Transform\n";
  wrl << "{\n";
  wrl << " children\n";
  wrl << " [\n";
  wrl << "  Shape\n";
  wrl << "  {\n";
  wrl << "   appearance Appearance\n";
  wrl << "   {\n";
  wrl << "    material Material\n";
  wrl << "    {\n";
  wrl << "     emissiveColor 1.0 1.0 1.0\n";
  wrl << "    } # material\n";
  wrl << "   } # appearance" << endl;
  wrl << "   geometry Text\n";
  wrl << "   {\n";
  wrl << "    string \"z\"\n";
  wrl << "    fontStyle FontStyle\n";
  wrl << "    {\n";
  wrl << "     family \"SANS\"\n";
  wrl << "     size 0.2\n";
  wrl << "    } # fontStyle\n";
  wrl << "   } # geometry\n";
  wrl << "  } # shape\n";
  wrl << " ] # children\n";
  wrl << " rotation 1.0 0.0 0.0 1.5708\n";
  wrl << " translation -0.15 0.0 0.9\n";
  wrl << "} # Transform\n";

  //wrl << "Transform {" << endl;
  //wrl << "children [" << endl;
  //wrl << "Shape {" << endl;
  //wrl << "appearance Appearance {" << endl;
  //wrl << "material Material {" << endl;
  //wrl << "emissiveColor 1.0 1.0 1.0 }}" << endl;
  //wrl << "geometry Text {" << endl;
  //wrl << "string \"y\"" << endl;
  //wrl << "fontStyle FontStyle {" << endl;
  //wrl << "family \"SANS\"" << endl;
  //wrl << "size 0.2 }}}]" << endl;
  //wrl << "translation -0.15 0.9 0.0 }" << endl;

  //wrl << "Transform {" << endl;
  //wrl << "children [" << endl;
  //wrl << "Shape {" << endl;
  //wrl << "appearance Appearance {" << endl;
  //wrl << "material Material {" << endl;
  //wrl << "emissiveColor 1.0 1.0 1.0 }}" << endl;
  //wrl << "geometry Text {" << endl;
  //wrl << "string \"z\"" << endl;
  //wrl << "fontStyle FontStyle {" << endl;
  //wrl << "family \"SANS\"" << endl;
  //wrl << "size 0.2 }}}]" << endl;
  //wrl << "rotation 1.0 0.0 0.0 1.5708" << endl;
  //wrl << "translation -0.15 0.0 0.9 }" << endl;

  // table of time keys
  vector<double> time_keys(m_steps_best+1, 0.0);
  for (i=0;i<=m_steps_best;++i)
  {
    time_keys[i] = (m_spacecraft_trajectory[i].MJD()-t_start)/dt_flight_days;
  } // for (i=0;i<=m_steps_best;++i)

  // planet trajectories
  if (isHeliocentric() == false) factor = 1.0e4;
  bool init_traj_drawn = false;
  bool tgt_traj_drawn = false;
  if (isHeliocentric())
  {
    const double a_init = m_solsys->a(getInitialBody(), MJD_J2000);
    const double a_tgt = m_solsys->a(getTargetBodyId(), MJD_J2000);
    const double a_stop = std::max(a_init, a_tgt);

    for (pind = MERCURY; pind <= PLUTO; ++pind)
    {
      if (m_solsys->a(pind, MJD_J2000) <= a_stop)
      {
        CBody* body = m_solsys->getBody(pind);
        double red,green,blue;
        if (pind == getTargetBodyId())
        {
          red = 0.8; green = 1.0; blue = 0.2;
          tgt_traj_drawn = true;
        }
        else if (pind == getInitialBody())
        {
          red = 0.4; green = 0.4; blue = 1.0;
          init_traj_drawn = true;
        }
        else
        {
          red = (1.0  - 0.001*double(pind*pind*pind));
          green = (0.7  - 0.05 *double(pind));
          blue = (-0.12 + 0.12 *double(pind));
        }
        body->printVrmlTrajectory(wrl, 9, 4, t_start, 100, red, green, blue);
      } // if (m_solsys->a(pind, MJD_J2000) <= a_stop)
    } // for (pind = MERCURY; pind <= pstop; ++pind)
  } // if (isHeliocentric())


  // trajectory of Initial Body
  if (init_traj_drawn==false)
  {
    m_solsys->getBody(getInitialBody())->printVrmlTrajectory(wrl, 9, 4,
      t_start, 100, 0.4, 0.4, 1.0);
    init_traj_drawn = true;
  } // if (init_traj_drawn==false)


  // 2nd Global Trajectory Optimization Competition
  if (getTrajectoryType() != GTOC2TRAJ)
  {
    // target body trajectory
    if (tgt_traj_drawn==false)
    {
      m_solsys->getBody(getTargetBodyId())->printVrmlTrajectory(wrl, 9, 4,
        t_end, 100, 0.8, 1.0, 0.2);
      tgt_traj_drawn = true;
    } // if (init_traj_drawn==false)
  } //   if (getTrajectoryType() != GTOC2TRAJ)
  else
  {
    for (unsigned int i=0;i<m_gtoc2_best_bodies.size();++i)
    {
      if (m_gtoc2_best_bodies[i]!=0)
      {
        BodyID id = m_gtoc2_best_bodies[i];
        m_solsys->getBody(id)->printVrmlTrajectory(wrl, 9, 4, t_end
          , 100, 0.8,1.0,0.2);
      } // if (m_gtoc2_best_bodies[i]!=0)
    } // for (int i=0;i<m_gtoc2_best_bodies.size();++i)
  }


  //wrl << "DEF InitialBodyTrajectory Shape\n";
  //wrl << "{\n";
  //wrl << " appearance Appearance\n";
  //wrl << " {\n";
  //wrl << "  material Material\n";
  //wrl << "  {\n";
  //wrl << "   emissiveColor 0.4 0.4 1.0\n";
  //wrl << "  } # material\n";
  //wrl << " } # appearance\n";
  //wrl << " geometry IndexedLineSet\n";
  //wrl << " {\n";
  //wrl << "  coord Coordinate\n";
  //wrl << "  {\n";
  //wrl << "   point\n";
  //wrl << "   [\n";
  //double t = MJD_J2000;
  //double period = m_solsys->OrbitPeriod(getInitialBody(), t);
  //const int steps = 100;
  //const double t_inc = (period*1.01)/double(steps);
  //wrl << scientific << setprecision(4) << showpos;
  //for (i=0; i<(steps+1); ++i,t+=t_inc)
  //{
  //  if ( (i%2) == 0 ) wrl << "\n    ";
  //  else wrl << "  ";
  //  wrl << setw(9) << m_solsys->x(getInitialBody(), t) << " ";
  //  wrl << setw(9) << m_solsys->y(getInitialBody(), t) << " ";
  //  wrl << setw(9) << m_solsys->z(getInitialBody(), t) << " ";
  //} //       for (i=0; i<(steps+1); ++i,t+=t_inc)
  //wrl << "\n   ] # point\n";
  //wrl << "  } # coord\n";
  //wrl << "  coordIndex\n";
  //wrl << "  [\n";
  //wrl << noshowpos;
  //for (i=0;i<(steps);++i)
  //{
  //  if (i%10==0) wrl << "\n   ";
  //  wrl << i << ",";          
  //} // for (i=0;i<steps;++i)
  //wrl << steps << endl;
  //wrl << "\n  ] # coordIndex\n";
  //wrl << " } # geometry\n";
  //wrl << "} # DEF InitialBodyTrajectory\n\n";
  //// Trajectory of Target Body
  //wrl << "DEF TargetBodyTrajectory Shape"                                            << endl;
  //wrl << "{"                                                << endl;
  //wrl << "  appearance Appearance"                          << endl;
  //wrl << "  {"                                              << endl;
  //wrl << "    material Material"                            << endl;
  //wrl << "    {"                                            << endl;
  //wrl << "      emissiveColor ";
  //wrl << "0.0 1.0 0.0" << endl;
  //wrl << "    }"                                            << endl;
  //wrl << "  }"                                              << endl;
  //wrl << "  geometry IndexedLineSet"                        << endl;
  //wrl << "  {"                                              << endl;
  //wrl << "    coord Coordinate"                             << endl;
  //wrl << "    { point ["                                    << endl;
  //period = m_solsys->OrbitPeriod(getTargetBodyId(), MJD_J2000);
  //for (t = MJD_J2000; t <= MJD_J2000+1.01*period; t += 1.01*period/100.0)
  //{
  //  wrl << factor*m_solsys->x(getTargetBodyId(), t) << " ";
  //  wrl << factor*m_solsys->y(getTargetBodyId(), t) << " ";
  //  wrl << factor*m_solsys->z(getTargetBodyId(), t)                          << endl;
  //}
  //wrl << "            ]"                                    << endl;
  //wrl << "    }"                                            << endl;
  //wrl << "    coordIndex ["                                 << endl;
  //i = 0;
  //for (t = MJD_J2000; t <= MJD_J2000+1.01*period; t += 1.01*period/100.0, ++i)
  //{
  //  if (t != MJD_J2000+1.01*period)
  //    wrl << i << ", ";
  //  else
  //    wrl << i;
  //}
  //wrl << "               ]"                                 << endl;
  //wrl << "  }"                                              << endl;
  //wrl << "}"   <<endl;


  // spacecraft trajectory
  wrl << "DEF SCTrajectory Shape\n";
  wrl << "{"                                                << endl;
  wrl << "  appearance Appearance"                          << endl;
  wrl << "  {"                                              << endl;
  wrl << "    material Material"                            << endl;
  wrl << "    {"                                            << endl;
  wrl << "      emissiveColor 1.0 1.0 0.0"                  << endl;
  wrl << "    }"                                            << endl;
  wrl << "  }"                                              << endl;
  wrl << "  geometry IndexedLineSet"                        << endl;
  wrl << "  {"                                              << endl;
  wrl << "    coord Coordinate"                             << endl;
  wrl << "    { point ["                                    << endl;
  for ( i=0; i <= m_steps_best; ++i)
  {
    wrl << factor*m_spacecraft_trajectory[i].CState::x() << " ";
    wrl << factor*m_spacecraft_trajectory[i].CState::y() << " ";
    wrl << factor*m_spacecraft_trajectory[i].CState::z() << endl;
  }
  wrl << "            ]"                                    << endl;
  wrl << "    }"                                            << endl;
  wrl << "    coordIndex ["                                 << endl;
  for (i=0; i <= m_steps_best; ++i)
  {
    if (i != m_steps_best)
      wrl << i << ", ";
    else
      wrl << i;
  }
  wrl << "               ]"                                 << endl;
  wrl << "  }"                                              << endl;
  wrl << "}"                                                << endl;

  // spacecraft body interpolator
  wrl << "\nDEF SCInterpolator PositionInterpolator"          << endl;
  wrl << "{"                                                << endl;
  wrl << "  key ["                                          << endl;
  for (i=0; i <= m_steps_best; ++i)
  {
    wrl << time_keys[i];
    if (i!=m_steps_best) wrl << ',';
  }
  wrl << "      ]"                                          << endl;
  wrl << "  keyValue ["                                     << endl;
  for (i=0; i <= m_steps_best; ++i)
  {
    wrl << factor*m_spacecraft_trajectory[i].CState::x() << " ";
    wrl << factor*m_spacecraft_trajectory[i].CState::y() << " ";
    wrl << factor*m_spacecraft_trajectory[i].CState::z();
    if (i != m_steps_best)
      wrl << "," << endl;
  }
  wrl << "           ]"                                     << endl;
  wrl << "}"                                                << endl;
  // spacecraft rotation interpolator
  wrl << "DEF SCRotInterpolator OrientationInterpolator"    << endl;
  wrl << "{"                                                << endl;
  wrl << "  key ["                                          << endl;
  for (i=0; i <= m_steps_best; ++i)
  {
    wrl << time_keys[i];
    if (i!=m_steps_best) wrl << ',';
  }
  wrl << "      ]"                                          << endl;
  wrl << "  keyValue ["                                     << endl;
  for (i=0; i <= m_steps_best; ++i)
  {
    const double vomega   = m_spacecraft_trajectory[i].CState::omega();
    const double vi       = m_spacecraft_trajectory[i].CState::i();
    const double vOmega   = m_spacecraft_trajectory[i].CState::Omega();
    const double vf       = m_spacecraft_trajectory[i].CState::f();
    const double valpha   = m_spacecraft_trajectory[i].getClockAngle();
    const double vbeta    = m_spacecraft_trajectory[i].getConeAngle();
    double x_axis[3];
    x_axis[0] = 1.0;
    x_axis[1] = 0.0;
    x_axis[2] = 0.0;
    double y_axis[3];
    y_axis[0] = 0.0;
    y_axis[1] = 1.0;
    y_axis[2] = 0.0;
    double z_axis[3];
    z_axis[0] = 0.0;
    z_axis[1] = 0.0;
    z_axis[2] = 1.0;
    double q_cylinder[4];
    ConvertToQuaternion(z_axis, -HALFPI, q_cylinder);
    double q_omega_[4];
    ConvertToQuaternion(z_axis, vomega, q_omega_);
    double q_omega[4];
    MultiplyQuaternion_(q_cylinder, q_omega_, q_omega);
    double q_i[4];
    ConvertToQuaternion(x_axis, vi, q_i);
    double q_omega_i[4];
    MultiplyQuaternion_(q_omega, q_i, q_omega_i);
    double q_Omega[4];
    ConvertToQuaternion(z_axis, vOmega, q_Omega);
    double q_omega_i_Omega[4];
    MultiplyQuaternion_(q_omega_i, q_Omega, q_omega_i_Omega);
    double h_axis[3];
    h_axis[0] = m_spacecraft_trajectory[i].CState::hVector()(0);
    h_axis[1] = m_spacecraft_trajectory[i].CState::hVector()(1);
    h_axis[2] = m_spacecraft_trajectory[i].CState::hVector()(2);
    double q_f_beta[4];
    ConvertToQuaternion(h_axis, vf+vbeta, q_f_beta);
    double q_omega_i_Omega_f_beta[4];
    MultiplyQuaternion_(q_omega_i_Omega,q_f_beta,q_omega_i_Omega_f_beta);
    double r_axis[3];
    r_axis[0] = m_spacecraft_trajectory[i].CState::rVector()(0);
    r_axis[1] = m_spacecraft_trajectory[i].CState::rVector()(1);
    r_axis[2] = m_spacecraft_trajectory[i].CState::rVector()(2);
    double q_alpha[4];
    ConvertToQuaternion(r_axis, valpha, q_alpha);
    double q_omega_i_Omega_f_beta_alpha[4];
    MultiplyQuaternion_(q_omega_i_Omega_f_beta, q_alpha,q_omega_i_Omega_f_beta_alpha);
    double axis[3];
    double angle;
    ConvertFromQuaternion(q_omega_i_Omega_f_beta_alpha, axis, &angle);
    wrl << axis[0] << " " << axis[1] << " " << axis[2] << " " << angle;
    if (i != m_steps_best)
      wrl << ", " << endl;
    else
      wrl << endl;
  }
  wrl << "           ]"                                     << endl;
  wrl << "}"                                                << endl;

  // spacecraft
  m_spacecraft.printVrmlAppearance(wrl);

  if (!m_spacecraft.isSolarSail())
  {
    // fuel scaling interpolator
    wrl << "DEF FuelInterpolator PositionInterpolator {" << endl;
    wrl << "key [" << endl;
    for (i=0; i <= m_steps_best; ++i)
    {
      wrl << time_keys[i];
      if (i!=m_steps_best) wrl << ',';
    }
    wrl << "]";
    wrl << "keyValue [" << endl;
    for (i=0; i <= m_steps_best; ++i)
    {
      wrl << 1.0 << " ";
      wrl << m_spacecraft_trajectory[i].getPropellantMass() / m_spacecraft_trajectory[i].getInitialPropellantMass() << " ";
      wrl << 1.0;
      if (i != m_steps_best)
        wrl << ", " << endl;
    }
    wrl << " ]}" << endl;
    wrl << "ROUTE time.fraction_changed TO FuelInterpolator.set_fraction" << endl;
    wrl << "ROUTE FuelInterpolator.value_changed TO Fuel.set_scale" << endl;
  }
  // arrow scaling interpolator
  wrl << "DEF ArrowInterpolator PositionInterpolator {" << endl;
  wrl << "key [" << endl;
  for (i=0; i <= m_steps_best; ++i)
  {
    wrl << time_keys[i];
    if (i!=m_steps_best) wrl << ',';
  }
  wrl << "]";
  wrl << "keyValue [" << endl;
  for (i=0; i <= m_steps_best; ++i)
  {
    wrl << 1.0 << " ";
    if (m_spacecraft_trajectory[i].MaxThrustAcceleration() != 0.0) // new due to error in VRML file (BD 03.11.05)
      wrl << m_spacecraft_trajectory[i].ThrustAcc()/m_spacecraft_trajectory[i].MaxThrustAcceleration() << " ";
    else
      wrl << 0.0 << " ";
    wrl << 1.0;
    if (i != m_steps_best)
      wrl << ", " << endl;
  }
  wrl << " ]}" << endl;
  // spacecraft routes
  wrl << "ROUTE time.fraction_changed TO ";
  wrl << getCurrentSC().CBody::getName();
  wrl << "Interpolator.set_fraction" << endl;
  wrl << "ROUTE SCInterpolator.value_changed TO SC.set_translation" << endl;
  wrl << "ROUTE time.fraction_changed TO SCRotInterpolator.set_fraction" << endl;
  wrl << "ROUTE SCRotInterpolator.value_changed TO SC.set_rotation" << endl;
  wrl << "ROUTE time.fraction_changed TO ArrowInterpolator.set_fraction" << endl;
  wrl << "ROUTE ArrowInterpolator.value_changed TO Arrow.set_scale" << endl;
  // initial body interpolator
  if (getInitialStateType() == BODY_STATE)
  {
    wrl << "DEF IBodyInterpolator PositionInterpolator"       << endl;
    wrl << "{"                                                << endl;
    wrl << "  key ["                                          << endl;
    for (i=0; i <= m_steps_best; ++i)
    {
      wrl << time_keys[i];
      if (i!=m_steps_best) wrl << ',';
    }
    wrl << "      ]"                                          << endl;
    wrl << "  keyValue ["                                     << endl;
    for (i=0; i <= m_steps_best; ++i)
    {
      wrl << m_solsys->x(getInitialBody(), 
        m_spacecraft_trajectory[i].CState::MJD()) << " ";
      wrl << m_solsys->y(getInitialBody(), 
        m_spacecraft_trajectory[i].CState::MJD()) << " ";
      wrl << m_solsys->z(getInitialBody(), 
        m_spacecraft_trajectory[i].CState::MJD());
      if (i != m_steps_best)
        wrl << "," << endl;
    }
    wrl << "           ]"                                     << endl;
    wrl << "}"                                                << endl;
    // initial body
    wrl << "DEF IBody Transform"                              << endl;
    wrl << "{"                                                << endl;
    wrl << "  children ["                                     << endl;
    wrl << "Shape"                                            << endl;
    wrl << "{"                                                << endl;
    wrl << "  appearance Appearance"                          << endl;
    wrl << "  {"                                              << endl;
    wrl << "    material Material"                            << endl;
    wrl << "    {"                                            << endl;
    wrl << "      emissiveColor 0.4 0.4 0.1"                  << endl;
    wrl << "      transparency 0.4"                           << endl;
    wrl << "    }"                                            << endl;
    wrl << "  }"                                              << endl;
    wrl << "  geometry Sphere"                                << endl;
    wrl << "  {"                                              << endl;
    wrl << "    radius ";
    double initial_hillradius = m_solsys->RHill(getInitialBody(), m_spacecraft_trajectory[0].CState::MJD());
    if (initial_hillradius < 0.01) initial_hillradius = 0.01;
    wrl << initial_hillradius << endl;
    wrl << "  }"                                              << endl;
    wrl << "}"                                                << endl;
    wrl << "           ]"                                     << endl;
    wrl << "}"                                                << endl;
    // initial body routes
    wrl << "ROUTE time.fraction_changed TO IBodyInterpolator.set_fraction" << endl;
    wrl << "ROUTE IBodyInterpolator.value_changed TO IBody.set_translation" << endl;
  }
  // target body interpolator
  if (getTargetStateType() == BODY_STATE ||
    getTargetStateType() == BODY_ORBIT_STATE ||
    getTargetStateType() == ORBITAL_ELEMENTS_STATE)
  {
    BodyID tgt_body = 0;
    if (getTrajectoryType() != GTOC2TRAJ)
    {
      tgt_body = getTargetBodyId();
    } // if (getTrajectoryType() != GTOC2TRAJ)
    else
    {
      tgt_body = m_gtoc2_best_bodies[0];
    }
    wrl << "DEF TBodyInterpolator PositionInterpolator"       << endl;
    wrl << "{"                                                << endl;
    wrl << "  key ["                                          << endl;
    for (i=0; i <= m_steps_best; ++i)
    {
      wrl << time_keys[i];
      if (i!=m_steps_best) wrl << ',';
    }
    wrl << "      ]"                                          << endl;
    wrl << "  keyValue ["                                     << endl;
    for (i=0; i <= m_steps_best; ++i)
    {
      wrl << factor*m_solsys->x(tgt_body,m_spacecraft_trajectory[i].CState::MJD()) << " ";
      wrl << factor*m_solsys->y(tgt_body,m_spacecraft_trajectory[i].CState::MJD()) << " ";
      wrl << factor*m_solsys->z(tgt_body,m_spacecraft_trajectory[i].CState::MJD());
      if (i != m_steps_best)
        wrl << "," << endl;
    }
    wrl << "           ]"                                     << endl;
    wrl << "}"                                                << endl;
    // target body
    wrl << "DEF TBody Transform"                              << endl;
    wrl << "{"                                                << endl;
    wrl << "  children ["                                     << endl;
    wrl << "Shape"                                            << endl;
    wrl << "{"                                                << endl;
    wrl << "  appearance Appearance"                          << endl;
    wrl << "  {"                                              << endl;
    wrl << "    material Material"                            << endl;
    wrl << "    {"                                            << endl;
    wrl << "      emissiveColor 0.8 1.0 0.1"                  << endl;
    wrl << "      transparency 0.4"                           << endl;
    wrl << "    }"                                            << endl;
    wrl << "  }"                                              << endl;
    wrl << "  geometry Sphere"                                << endl;
    wrl << "  {"                                              << endl;
    wrl << "    radius ";
    double target_hillradius = std::max(0.01, m_solsys->RHill(tgt_body, t_start));
    wrl << target_hillradius << endl;
    wrl << "  }"                                              << endl;
    wrl << "}"                                                << endl;
    wrl << "           ]"                                     << endl;
    wrl << "}"                                                << endl;
    // target body routes
    wrl << "ROUTE time.fraction_changed TO TBodyInterpolator.set_fraction" << endl;
    wrl << "ROUTE TBodyInterpolator.value_changed TO TBody.set_translation" << endl;
  }
  // planets body interpolator
  if (isHeliocentric())
  {
    if (m_solsys->a(getTargetBodyId(),MJD_J2000)<m_solsys->a(MARS,MJD_J2000))
      pstop=MARS;
    else
      pstop=PLUTO;
    for (pind = MERCURY; pind <= pstop; ++pind)
    {
      if (pind != static_cast<int>(getInitialBody()) && pind != static_cast<int>(getTargetBodyId()))
      {
        wrl << "DEF ";
        wrl << m_solsys->getBodyName(pind);
        wrl << "BodyInterpolator PositionInterpolator\n";
        wrl << "{"                                                << endl;
        wrl << "  key ["                                          << endl;
        for (i=0; i <= m_steps_best; ++i)
        {
          wrl << time_keys[i];
          if (i!=m_steps_best) wrl << ',';
        }
        wrl << "      ]"                                          << endl;
        wrl << "  keyValue ["                                     << endl;
        for (i=0; i <= m_steps_best; ++i)
        {
          wrl << m_solsys->x(pind, m_spacecraft_trajectory[i].CState::MJD()) << " ";
          wrl << m_solsys->y(pind, m_spacecraft_trajectory[i].CState::MJD()) << " ";
          wrl << m_solsys->z(pind, m_spacecraft_trajectory[i].CState::MJD());
          if (i != m_steps_best)
            wrl << "," << endl;
        }
        wrl << "           ]"                                     << endl;
        wrl << "}"                                                << endl;
        // planets body
        wrl << "DEF ";
        wrl << m_solsys->getBodyName(pind);
        wrl << "Body Transform" << endl;
        wrl << "{"                                                << endl;
        wrl << "  children ["                                     << endl;
        wrl << "Shape"                                            << endl;
        wrl << "{"                                                << endl;
        wrl << "  appearance Appearance"                          << endl;
        wrl << "  {"                                              << endl;
        wrl << "    material Material"                            << endl;
        wrl << "    {"                                            << endl;
        wrl << "      emissiveColor ";
        wrl <<  1.0  - 0.001*double(pind*pind*pind)               << " "; // red
        wrl <<  0.7  - 0.05 *double(pind)                         << " "; // green
        wrl << -0.12 + 0.12 *double(pind)                         << endl; // blue
        wrl << "      transparency 0.4"                           << endl;
        wrl << "    }"                                            << endl;
        wrl << "  }"                                              << endl;
        wrl << "  geometry Sphere"                                << endl;
        wrl << "  {"                                              << endl;
        wrl << "    radius ";
        double initial_hillradius = m_solsys->RHill(pind, m_spacecraft_trajectory[0].MJD());
        if (initial_hillradius < 0.01) initial_hillradius = 0.01;
        wrl << initial_hillradius << endl;
        wrl << "  }"                                              << endl;
        wrl << "}"                                                << endl;
        wrl << "           ]"                                     << endl;
        wrl << "}"                                                << endl;
        // planets body routes
        wrl << "ROUTE time.fraction_changed TO ";
        wrl << m_solsys->getBodyName(pind);
        wrl << "BodyInterpolator.set_fraction" << endl;
      } // end if
    } // end for
  } // if (isHeliocentric())

} // end CSimulation::PrintVRML


// helper function to improve readability of code
inline void VrmlAxisAngleOutput(Vector& q, int i, std::ofstream& file)
{
  Vector axis(3);
  double angle(0.0);

  Quaternion2AxisAngle(q[i],q[i+1],q[i+2],q[i+3], axis, angle);
  for (int j=0;j<3;++j) file << axis[j] << " ";
  file << angle;
}


// ----------------------------------------------------------------------------
void CSimulation::PrintANG(const char*filename)
// ----------------------------------------------------------------------------
{
  ofstream ang(filename);
  Assure(ang, filename);

  ang<<m_steps_best+1<<endl;

  for (int i=0; i<=m_steps_best; i++) {
    // time
    ang<<setprecision(12)<<m_spacecraft_trajectory[i].MJD()<<'\t';
    // clock angle
    ang<<setprecision(12)<<m_spacecraft_trajectory[i].getClockAngle()<<'\t';
    // cone angle
    ang<<setprecision(12)<<m_spacecraft_trajectory[i].getConeAngle()<<endl;
  }
}; // end CSimulation::PrintANG


//------------------------------------------------------------------------------
void CSimulation::saveApproximationClockTime(const std::string& filename) const
{
  if (filename.empty()!=true)
  {
    std::ofstream o(filename.c_str());
    Global::FileOutputEndl(m_best_clock_knot_time, o);
  }
} // void CSimulation::saveApproximationClockTime(const char* filename) const


//------------------------------------------------------------------------------
void CSimulation::saveApproximationClockValues(const std::string& filename) const
{
  if (filename.empty()!=true)
  {
    std::ofstream o(filename.c_str());
    Global::FileOutputEndl(m_best_clock_knot_values, o);
  }
} // void CSimulation::saveApproximationClockValues(const char* filename) const


//------------------------------------------------------------------------------
void CSimulation::saveApproxParams(
                           const std::vector<double>& kt ///< knot time
                           , const std::vector<double>& kv ///< knot value
                           , const std::string& filename     ///< filename
                           , unsigned int sz                ///< # of values
                           ) const
{
  if (filename.empty() != true)
  {
    if ( (sz<=kt.size()) && (sz<=kv.size()) )
    {
      std::ofstream ofs(filename.c_str(), ios_base::binary);
      if (ofs.is_open() == true)
      {
        ofs.write((const char*)&sz, sizeof(unsigned int));// write # of data sets
        double d(0.0);
        for (unsigned int i=0;i<sz;++i) // write data
        {
          d = m_best_clock_knot_time[i];
          ofs.write((const char*)&d, sizeof(double));
          d = m_best_clock_knot_values[i];
          ofs.write((const char*)&d, sizeof(double));
        }
      }
    }
  }
}//void CSimulation::saveClockApproxParams(const string&, const valarray....)


//------------------------------------------------------------------------------
void CSimulation::saveApproximationConeTime(const std::string& filename) const
{
  if (filename.empty()!=true)
  {
    std::ofstream o(filename.c_str());
    Global::FileOutputEndl(m_best_cone_knot_time, o);
  }
} // void CSimulation::saveApproximationClockTime(const char* filename) const


//------------------------------------------------------------------------------
void CSimulation::saveApproximationConeValues(const std::string& filename) const
{
  if (filename.empty()!=true)
  {
    std::ofstream o(filename.c_str());
    Global::FileOutputEndl(m_best_cone_knot_values, o);
  }
} // void CSimulation::saveApproximationClockValues(const char* filename) const


//------------------------------------------------------------------------------
void CSimulation::swapReferenceTrajectory(SCTraj& t)
{
  m_ref_sc_traj.swap(t);
  const unsigned int n_ref_steps = m_ref_sc_traj.size();

  const unsigned int n_max_approx_steps =
    unsigned int(ceil(n_ref_steps*m_approx_res_factor));

  m_last_saved_trajectory_approx.resize(n_max_approx_steps);
  m_last_saved_clock_angle_approx_values.resize(n_max_approx_steps);
  m_last_saved_clock_angle_approx_values_dot.resize(n_max_approx_steps);
  m_last_saved_cone_angle_approx_values.resize(n_max_approx_steps);
  m_last_saved_cone_angle_approx_values_dot.resize(n_max_approx_steps);
} // void CSimulation::swapReferenceTrajectory(SCTraj& t)


//------------------------------------------------------------------------------
void CSimulation::setReferenceTrajectory(const SCTraj& t)
{
#if defined COMMONCPP
  SimpleGuard g(m_ref_sc_traj_lock);
#endif //COMMONCPP
  m_ref_sc_traj=t;
  m_last_saved_trajectory_approx.resize(m_ref_sc_traj.size());
  m_last_saved_clock_angle_approx_values.resize(m_ref_sc_traj.size());
  m_last_saved_clock_angle_approx_values_dot.resize(m_ref_sc_traj.size());
  m_last_saved_cone_angle_approx_values.resize(m_ref_sc_traj.size());
  m_last_saved_cone_angle_approx_values_dot.resize(m_ref_sc_traj.size());
} // void CSimulation::setReferenceTrajectory(const SCTraj& t)


//------------------------------------------------------------------------------
Trajectory CSimulation::getTrajectoryData(const SCTraj& source) const
{
  if (source.empty() == true) return 0;   // check for empty source data
  Trajectory traj = 0;
  try
  {
    traj = new CTrajectory();
  }
  catch (std::bad_alloc&) { traj = 0;}
  if (traj != 0)
  {
    getTrajectoryData(source, *traj, source.size());
    if (traj->empty() == true) 
    {
      delete traj;
      traj = 0;
    }
  }
  else Output("Trajectory CSimulation::getTrajectoryData(const SCTraj& source)\
              : trajectory data couldn't be obtained!\n");
  return traj;
} // Trajectory CSimulation::getTrajectoryData(const SCTraj& source) const
