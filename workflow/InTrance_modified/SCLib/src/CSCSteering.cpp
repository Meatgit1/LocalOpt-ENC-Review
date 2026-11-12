/*!
\file CSCSteering.cpp
\brief Two classes that provide the functionality for steering a spacecraft.

\version

  1.7

\author

  Bernd Dachwald\n
  German Aerospace Center (DLR)\n
  Institute of Space Simulation\n
  Linder Hoehe\n
  51147 Cologne\n
  Germany\n
  bernd.dachwald@dlr.de

\date

  17.09.2006 (AO)

\b Compiler:

  Microsoft Visual C++ .NET

\b Copyright:

  (c) 2001-2006 by Bernd Dachwald,Andreas Ohndorf

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/

//SCLib
#include <configSC.h>
#include <algorithm>

/*
#ifdef BORLAND
 #include <stdlib.h>
#else
 #include <minmax.h>
#endif
*/
#include <Kepler.h>
#include <CSC.h>
#include <CNC.h>


using namespace std;
using namespace Global;
using namespace SpaceLib;

// ----------------------------------------------------------------------------
int CSCSteering::getNumberOfLSL(void) const
// ----------------------------------------------------------------------------
// private
{
  return m_number_of_lsl;

}; // end CSCSteering::getNumberOfLSL


// ----------------------------------------------------------------------------
int CSCSteering::getNumberOfLStSt(void) const
// ----------------------------------------------------------------------------
// private
{
  return m_number_of_stst;

}; // end CSCSteering::getNumberOfLStSt


// ----------------------------------------------------------------------------
StStIndex CSCSteering::getSteeringStrategy(void) const
// ----------------------------------------------------------------------------
// private
{
  return m_stst_index;

}; // end CSCSteering::getSteeringStrategy


// ----------------------------------------------------------------------------
void CSCSteering::setSteeringStrategy(const StStIndex stst)
// ----------------------------------------------------------------------------
// private
{
  m_stst_index = stst;

}; // end CSCSteering::setSteeringStrategy


// ----------------------------------------------------------------------------
// private
void CSCSteering::setNCInputCartesian(int           first_neuron,
                                      int           last_neuron,
                                      const CState& ps,
                                      const CState& ts,
                                      CSC*          sc,
                                      const CState* ga) {
//AO

  const double ps_dis_factor = (ps.Unit() == SI ? M2AU : 1.0);
  const double ts_dis_factor = (ts.Unit() == SI ? M2AU : 1.0);
  const double ga_dis_factor = ((ga!=0?ga->Unit():AUD) == SI ? M2AU : 1.0);
  const double ps_vel_factor = (ps.Unit() == SI ? MS2AUD : 1.0);
  const double ts_vel_factor = (ts.Unit() == SI ? MS2AUD : 1.0);
  const double ga_vel_factor = ((ga!=0?ga->Unit():AUD) == SI ? MS2AUD : 1.0);

//AO

  int d = last_neuron - first_neuron + 1;

  bool delta_steering;
  if (d == 6 || (d == 12 && m_ga_steering)) delta_steering = true;
   else if ((d == 12 && !m_ga_steering) || d == 18) delta_steering = false;
   else FatalError("CSCSteering::setNCInputCartesian(): error in function input");

  // d == 6
  if (delta_steering && !m_ga_steering) {

    sc->m_input_ncp(first_neuron     ) =  ps.x()*ps_dis_factor  - ts.x()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  1) =  ps.y()*ps_dis_factor  - ts.y()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  2) =  ps.z()*ps_dis_factor  - ts.z()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  3) = (ps.dx()*ps_vel_factor - ts.dx()*ts_vel_factor) * 100.0;
    sc->m_input_ncp(first_neuron +  4) = (ps.dy()*ps_vel_factor - ts.dy()*ts_vel_factor) * 100.0;
    sc->m_input_ncp(first_neuron +  5) = (ps.dz()*ps_vel_factor - ts.dz()*ts_vel_factor) * 100.0;

  } else if (delta_steering && m_ga_steering && ga)  { // d == 12

    sc->m_input_ncp(first_neuron     ) =  ps.x()*ps_dis_factor  - ts.x()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  1) =  ps.y()*ps_dis_factor  - ts.y()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  2) =  ps.z()*ps_dis_factor  - ts.z()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  3) = (ps.dx()*ps_vel_factor - ts.dx()*ts_vel_factor) * 100.0;
    sc->m_input_ncp(first_neuron +  4) = (ps.dy()*ps_vel_factor - ts.dy()*ts_vel_factor) * 100.0;
    sc->m_input_ncp(first_neuron +  5) = (ps.dz()*ps_vel_factor - ts.dz()*ts_vel_factor) * 100.0;
    sc->m_input_ncp(first_neuron +  6) =  ts.x()*ts_dis_factor  - ga->x()*ga_dis_factor;
    sc->m_input_ncp(first_neuron +  7) =  ts.y()*ts_dis_factor  - ga->y()*ga_dis_factor;
    sc->m_input_ncp(first_neuron +  8) =  ts.z()*ts_dis_factor  - ga->z()*ga_dis_factor;
    sc->m_input_ncp(first_neuron +  9) = (ts.dx()*ts_vel_factor - ga->dx()*ga_vel_factor) * 100.0;
    sc->m_input_ncp(first_neuron + 10) = (ts.dy()*ts_vel_factor - ga->dy()*ga_vel_factor) * 100.0;
    sc->m_input_ncp(first_neuron + 11) = (ts.dz()*ts_vel_factor - ga->dz()*ga_vel_factor) * 100.0;

  } else if (delta_steering && m_ga_steering && !ga) { // d == 12

    sc->m_input_ncp(first_neuron     ) =  ps.x()  - ts.x()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  1) =  ps.y()  - ts.y()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  2) =  ps.z()  - ts.z()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  3) = (ps.dx()*ps_vel_factor - ts.dx()*ts_vel_factor) * 100.0;
    sc->m_input_ncp(first_neuron +  4) = (ps.dy()*ps_vel_factor - ts.dy()*ts_vel_factor) * 100.0;
    sc->m_input_ncp(first_neuron +  5) = (ps.dz()*ps_vel_factor - ts.dz()*ts_vel_factor) * 100.0;
    sc->m_input_ncp(first_neuron +  6) = 0.0;
    sc->m_input_ncp(first_neuron +  7) = 0.0;
    sc->m_input_ncp(first_neuron +  8) = 0.0;
    sc->m_input_ncp(first_neuron +  9) = 0.0;
    sc->m_input_ncp(first_neuron + 10) = 0.0;
    sc->m_input_ncp(first_neuron + 11) = 0.0;

  } else if (!delta_steering && !m_ga_steering) { // d == 12

    sc->m_input_ncp(first_neuron     ) = ps.x()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  1) = ps.y()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  2) = ps.z()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  3) = ps.dx()*ps_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron +  4) = ps.dy()*ps_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron +  5) = ps.dz()*ps_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron +  6) = ts.x()*ts_vel_factor;
    sc->m_input_ncp(first_neuron +  7) = ts.y()*ts_vel_factor;
    sc->m_input_ncp(first_neuron +  8) = ts.z()*ts_vel_factor;
    sc->m_input_ncp(first_neuron +  9) = ts.dx()*ts_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron + 10) = ts.dy()*ts_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron + 11) = ts.dz()*ts_vel_factor * 100.0;

  }  else if (!delta_steering && m_ga_steering && ga) { // d == 18

    sc->m_input_ncp(first_neuron     ) = ps.x()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  1) = ps.y()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  2) = ps.z()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  3) = ps.dx()*ps_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron +  4) = ps.dy()*ps_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron +  5) = ps.dz()*ps_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron +  6) = ts.x()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  7) = ts.y()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  8) = ts.z()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  9) = ts.dx()*ts_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron + 10) = ts.dy()*ts_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron + 11) = ts.dz()*ts_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron + 12) = ga->x()*ga_dis_factor;
    sc->m_input_ncp(first_neuron + 13) = ga->y()*ga_dis_factor;
    sc->m_input_ncp(first_neuron + 14) = ga->z()*ga_dis_factor;
    sc->m_input_ncp(first_neuron + 15) = ga->dx()*ga_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron + 16) = ga->dy()*ga_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron + 17) = ga->dz()*ga_vel_factor * 100.0;

  } else if (!delta_steering && m_ga_steering && !ga) { // d == 18

    sc->m_input_ncp(first_neuron     ) = ps.x()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  1) = ps.y()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  2) = ps.z()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  3) = ps.dx()*ps_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron +  4) = ps.dy()*ps_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron +  5) = ps.dz()*ps_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron +  6) = ts.x()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  7) = ts.y()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  8) = ts.z()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  9) = ts.dx()*ts_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron + 10) = ts.dy()*ts_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron + 11) = ts.dz()*ts_vel_factor * 100.0;
    sc->m_input_ncp(first_neuron + 12) = 0.0;
    sc->m_input_ncp(first_neuron + 13) = 0.0;
    sc->m_input_ncp(first_neuron + 14) = 0.0;
    sc->m_input_ncp(first_neuron + 15) = 0.0;
    sc->m_input_ncp(first_neuron + 16) = 0.0;
    sc->m_input_ncp(first_neuron + 17) = 0.0;

  } else FatalError("CSCSteering::setNCInputCartesian(): error in function input");

}; // end CSCSteering::setNCInputCartesian


// ----------------------------------------------------------------------------
// private
void CSCSteering::setNCInputPolar(int           first_neuron,
                                  int           last_neuron,
                                  const CState& ps,
                                  const CState& ts,
                                  CSC*          sc,
                                  const CState* ga) {

//AO
  const double ps_dis_factor = (ps.Unit() == SI ? M2AU : 1.0);
  const double ts_dis_factor = (ts.Unit() == SI ? M2AU : 1.0);
  const double ga_dis_factor = ((ga!=0?ga->Unit():AUD) == SI ? M2AU : 1.0);
  const double ps_vel_factor = (ps.Unit() == SI ? MS2AUD : 1.0);
  const double ts_vel_factor = (ts.Unit() == SI ? MS2AUD : 1.0);
  const double ga_vel_factor = ((ga!=0?ga->Unit():AUD) == SI ? MS2AUD : 1.0);
//AO

  double dr_factor = 100.0;
  double dphi_factor = 100.0;
  double dtheta_factor = 100.0;
  if (sc->CState::isBodycentric())
  {
    dr_factor = 1.0;
    dphi_factor = 1.0;
    dtheta_factor = 1.0;
  }

  int d = last_neuron - first_neuron + 1;

  bool delta_steering;

  if (d == 6 || (d == 12 && m_ga_steering)) delta_steering = true;
   else if ((d == 12 && !m_ga_steering) || d == 18) delta_steering = false;
   else FatalError("CSCSteering::setNCInputPolar(): error in function input");

  if (delta_steering && !m_ga_steering) { //d == 6

    sc->m_input_ncp(first_neuron     ) =  ps.r()*ps_dis_factor - ts.r()*ts_dis_factor;

    sc->m_input_ncp(first_neuron +  1) =  ps.phi()    - ts.phi();
    AssertAngularRange(sc->m_input_ncp(first_neuron+1));

    sc->m_input_ncp(first_neuron +  2) =  ps.theta()  - ts.theta();
    AssertAngularRange(sc->m_input_ncp(first_neuron+2));

    sc->m_input_ncp(first_neuron +  3) = (ps.dr()*ps_vel_factor - ts.dr()*ts_vel_factor)*dr_factor;
    sc->m_input_ncp(first_neuron +  4) = (ps.dphi()   - ts.dphi())   * dphi_factor;
    sc->m_input_ncp(first_neuron +  5) = (ps.dtheta() - ts.dtheta()) * dtheta_factor;

  } else if (delta_steering && m_ga_steering && ga) { // d == 12

    sc->m_input_ncp(first_neuron     ) =  ps.r()*ps_dis_factor - ts.r()*ts_dis_factor;

    sc->m_input_ncp(first_neuron +  1) =  ps.phi()    - ts.phi();
    AssertAngularRange(sc->m_input_ncp(first_neuron+1));

    sc->m_input_ncp(first_neuron +  2) =  ps.theta()  - ts.theta();
    AssertAngularRange(sc->m_input_ncp(first_neuron+2));

    sc->m_input_ncp(first_neuron +  3) = (ps.dr()*ps_vel_factor - ts.dr()*ts_vel_factor)*dr_factor;
    sc->m_input_ncp(first_neuron +  4) = (ps.dphi()   - ts.dphi())   * dphi_factor;
    sc->m_input_ncp(first_neuron +  5) = (ps.dtheta() - ts.dtheta()) * dtheta_factor;

    sc->m_input_ncp(first_neuron +  6) =  ts.r()*ps_dis_factor - ga->r()*ga_dis_factor;
    sc->m_input_ncp(first_neuron +  7) =  ts.phi()    - ga->phi();
    AssertAngularRange(sc->m_input_ncp(first_neuron+7));

    sc->m_input_ncp(first_neuron +  8) =  ts.theta()  - ga->theta();
    AssertAngularRange(sc->m_input_ncp(first_neuron+8));

    sc->m_input_ncp(first_neuron +  9) = (ts.dr()*ts_dis_factor - ga->dr()*ts_vel_factor)*dr_factor;
    sc->m_input_ncp(first_neuron + 10) = (ts.dphi()   - ga->dphi())   * dphi_factor;
    sc->m_input_ncp(first_neuron + 11) = (ts.dtheta() - ga->dtheta()) * dtheta_factor;

  } else if (delta_steering && m_ga_steering && !ga) {// d == 12

    sc->m_input_ncp(first_neuron     ) =  ps.r()*ps_dis_factor - ts.r()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  1) =  ps.phi() - ts.phi();
    AssertAngularRange(sc->m_input_ncp(first_neuron+1));

    sc->m_input_ncp(first_neuron +  2) =  ps.theta()  - ts.theta();
    AssertAngularRange(sc->m_input_ncp(first_neuron+2));

    sc->m_input_ncp(first_neuron +  3) = (ps.dr()*ps_vel_factor - ts.dr()*ts_vel_factor)* dr_factor;
    sc->m_input_ncp(first_neuron +  4) = (ps.dphi()   - ts.dphi())   * dphi_factor;
    sc->m_input_ncp(first_neuron +  5) = (ps.dtheta() - ts.dtheta()) * dtheta_factor;

    sc->m_input_ncp(first_neuron +  6) =  0.0;
    sc->m_input_ncp(first_neuron +  7) =  0.0;
    sc->m_input_ncp(first_neuron +  8) =  0.0;
    sc->m_input_ncp(first_neuron +  9) =  0.0;
    sc->m_input_ncp(first_neuron + 10) =  0.0;
    sc->m_input_ncp(first_neuron + 11) =  0.0;

  } else if (!delta_steering && !m_ga_steering) { // d == 12

    sc->m_input_ncp(first_neuron     ) = ps.r()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  1) = ps.phi();
    sc->m_input_ncp(first_neuron +  2) = ps.theta();
    sc->m_input_ncp(first_neuron +  3) = ps.dr()*ps_vel_factor * dr_factor;
    sc->m_input_ncp(first_neuron +  4) = ps.dphi()   * dphi_factor;
    sc->m_input_ncp(first_neuron +  5) = ps.dtheta() * dtheta_factor;
    sc->m_input_ncp(first_neuron +  6) = ts.r()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  7) = ts.phi();
    sc->m_input_ncp(first_neuron +  8) = ts.theta();
    sc->m_input_ncp(first_neuron +  9) = ts.dr()*ts_vel_factor * dr_factor;
    sc->m_input_ncp(first_neuron + 10) = ts.dphi()   * dphi_factor;
    sc->m_input_ncp(first_neuron + 11) = ts.dtheta() * dtheta_factor;

  } else if (!delta_steering && m_ga_steering && ga) { // d == 18

    sc->m_input_ncp(first_neuron     ) = ps.r()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  1) = ps.phi();
    sc->m_input_ncp(first_neuron +  2) = ps.theta();
    sc->m_input_ncp(first_neuron +  3) = ps.dr()*ps_vel_factor * dr_factor;
    sc->m_input_ncp(first_neuron +  4) = ps.dphi()   * dphi_factor;
    sc->m_input_ncp(first_neuron +  5) = ps.dtheta() * dtheta_factor;
    sc->m_input_ncp(first_neuron +  6) = ts.r()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  7) = ts.phi();
    sc->m_input_ncp(first_neuron +  8) = ts.theta();
    sc->m_input_ncp(first_neuron +  9) = ts.dr()*ts_vel_factor*dr_factor;
    sc->m_input_ncp(first_neuron + 10) = ts.dphi()*dphi_factor;
    sc->m_input_ncp(first_neuron + 11) = ts.dtheta()*dtheta_factor;
    sc->m_input_ncp(first_neuron + 12) = ga->r()*ga_dis_factor;
    sc->m_input_ncp(first_neuron + 13) = ga->phi();
    sc->m_input_ncp(first_neuron + 14) = ga->theta();
    sc->m_input_ncp(first_neuron + 15) = ga->dr()*ga_vel_factor*dr_factor;
    sc->m_input_ncp(first_neuron + 16) = ga->dphi()*dphi_factor;
    sc->m_input_ncp(first_neuron + 17) = ga->dtheta()*dtheta_factor;

  } else if (!delta_steering && m_ga_steering && !ga) {// d == 18

    sc->m_input_ncp(first_neuron     ) = ps.r()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  1) = ps.phi();
    sc->m_input_ncp(first_neuron +  2) = ps.theta();
    sc->m_input_ncp(first_neuron +  3) = ps.dr()*ps_vel_factor * dr_factor;
    sc->m_input_ncp(first_neuron +  4) = ps.dphi()* dphi_factor;
    sc->m_input_ncp(first_neuron +  5) = ps.dtheta()*dtheta_factor;
    sc->m_input_ncp(first_neuron +  6) = ts.r()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  7) = ts.phi();
    sc->m_input_ncp(first_neuron +  8) = ts.theta();
    sc->m_input_ncp(first_neuron +  9) = ts.dr()*ts_vel_factor*dr_factor;
    sc->m_input_ncp(first_neuron + 10) = ts.dphi()*dphi_factor;
    sc->m_input_ncp(first_neuron + 11) = ts.dtheta()*dtheta_factor;
    sc->m_input_ncp(first_neuron + 12) = 0.0;
    sc->m_input_ncp(first_neuron + 13) = 0.0;
    sc->m_input_ncp(first_neuron + 14) = 0.0;
    sc->m_input_ncp(first_neuron + 15) = 0.0;
    sc->m_input_ncp(first_neuron + 16) = 0.0;
    sc->m_input_ncp(first_neuron + 17) = 0.0;
  }
  else
    FatalError("CSCSteering::setNCInputPolar(): error in function input");

}; // end CSCSteering::setNCInputPolar


// ----------------------------------------------------------------------------
// private
void CSCSteering::setNCInputKepler(int           first_neuron,
                                   int           last_neuron,
                                   const CState& ps,
                                   const CState& ts,
                                   CSC*          sc,
                                   const CState* ga)
{

//AO
  const double ps_dis_factor = (ps.Unit() == SI ? M2AU : 1.0);
  const double ts_dis_factor = (ts.Unit() == SI ? M2AU : 1.0);
  const double ga_dis_factor = ((ga!=0?ga->Unit():AUD) == SI ? M2AU : 1.0);
  const double ps_vel_factor = (ps.Unit() == SI ? MS2AUD : 1.0);
  const double ts_vel_factor = (ts.Unit() == SI ? MS2AUD : 1.0);
  const double ga_vel_factor = ((ga!=0?ga->Unit():AUD) == SI ? MS2AUD : 1.0);
//AO

  int d = last_neuron - first_neuron + 1;

  bool delta_steering;
  if (d == 6 || (d == 12 && m_ga_steering)) delta_steering = true;
  else if ((d == 12 && !m_ga_steering) || d == 18) delta_steering = false;
  else FatalError("CSCSteering::setNCInputKepler(): error in function input");

  if (delta_steering && !m_ga_steering) {// d == 6

    sc->m_input_ncp(first_neuron     ) = ps.a()*ps_dis_factor - ts.a()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  1) = ps.e()     - ts.e();
    sc->m_input_ncp(first_neuron +  2) = ps.i()     - ts.i();
    AssertAngularRange(sc->m_input_ncp(first_neuron+2));

    sc->m_input_ncp(first_neuron +  3) = ps.Omega() - ts.Omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron+3));

    sc->m_input_ncp(first_neuron +  4) = ps.omega() - ts.omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron+4));

    sc->m_input_ncp(first_neuron +  5) = ps.f()     - ts.f();
    AssertAngularRange(sc->m_input_ncp(first_neuron+5));

  } else if (delta_steering && m_ga_steering && ga) {// d == 12

    sc->m_input_ncp(first_neuron     ) = ps.a()*ps_dis_factor - ts.a()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  1) = ps.e()     - ts.e();
    sc->m_input_ncp(first_neuron +  2) = ps.i()     - ts.i();
    AssertAngularRange(sc->m_input_ncp(first_neuron+2));

    sc->m_input_ncp(first_neuron +  3) = ps.Omega() - ts.Omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron+3));

    sc->m_input_ncp(first_neuron +  4) = ps.omega() - ts.omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron+4));

    sc->m_input_ncp(first_neuron +  5) = ps.f()     - ts.f();
    AssertAngularRange(sc->m_input_ncp(first_neuron+5));

    sc->m_input_ncp(first_neuron +  6) = ts.a()*ts_dis_factor - ga->a()*ga_dis_factor;
    sc->m_input_ncp(first_neuron +  7) = ts.e()     - ga->e();
    sc->m_input_ncp(first_neuron +  8) = ts.i()     - ga->i();
    AssertAngularRange(sc->m_input_ncp(first_neuron+8));

    sc->m_input_ncp(first_neuron +  9) = ts.Omega() - ga->Omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron+9));

    sc->m_input_ncp(first_neuron + 10) = ts.omega() - ga->omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron+10));

    sc->m_input_ncp(first_neuron + 11) = ts.f()     - ga->f();
    AssertAngularRange(sc->m_input_ncp(first_neuron+11));

  } else if (delta_steering && m_ga_steering && !ga) { // d == 12

    sc->m_input_ncp(first_neuron     ) = ps.a()*ps_dis_factor - ts.a()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  1) = ps.e()     - ts.e();
    sc->m_input_ncp(first_neuron +  2) = ps.i()     - ts.i();
    AssertAngularRange(sc->m_input_ncp(first_neuron+2));

    sc->m_input_ncp(first_neuron +  3) = ps.Omega() - ts.Omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron+3));

    sc->m_input_ncp(first_neuron +  4) = ps.omega() - ts.omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron+4));

    sc->m_input_ncp(first_neuron +  5) = ps.f()     - ts.f();
    AssertAngularRange(sc->m_input_ncp(first_neuron+5));

    sc->m_input_ncp(first_neuron +  6) = 0.0;
    sc->m_input_ncp(first_neuron +  7) = 0.0;
    sc->m_input_ncp(first_neuron +  8) = 0.0;
    sc->m_input_ncp(first_neuron +  9) = 0.0;
    sc->m_input_ncp(first_neuron + 10) = 0.0;
    sc->m_input_ncp(first_neuron + 11) = 0.0;

  } else if (!delta_steering && !m_ga_steering) {// d == 12

    sc->m_input_ncp(first_neuron     ) = ps.a()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  1) = ps.e();
    sc->m_input_ncp(first_neuron +  2) = ps.i();
    sc->m_input_ncp(first_neuron +  3) = ps.Omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron +  3));

    sc->m_input_ncp(first_neuron +  4) = ps.omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron +  4));

    sc->m_input_ncp(first_neuron +  5) = ps.f();
    AssertAngularRange(sc->m_input_ncp(first_neuron +  5));

    sc->m_input_ncp(first_neuron +  6) = ts.a()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  7) = ts.e();
    sc->m_input_ncp(first_neuron +  8) = ts.i();
    sc->m_input_ncp(first_neuron +  9) = ts.Omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron +  9));

    sc->m_input_ncp(first_neuron + 10) = ts.omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron + 10));

    sc->m_input_ncp(first_neuron + 11) = ts.f();
    AssertAngularRange(sc->m_input_ncp(first_neuron + 11));

  } else if (!delta_steering && m_ga_steering && ga) {// d == 18

    sc->m_input_ncp(first_neuron     ) = ps.a()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  1) = ps.e();
    sc->m_input_ncp(first_neuron +  2) = ps.i();
    sc->m_input_ncp(first_neuron +  3) = ps.Omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron +  3));

    sc->m_input_ncp(first_neuron +  4) = ps.omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron +  4));

    sc->m_input_ncp(first_neuron +  5) = ps.f();
    AssertAngularRange(sc->m_input_ncp(first_neuron +  5));

    sc->m_input_ncp(first_neuron +  6) = ts.a()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  7) = ts.e();
    sc->m_input_ncp(first_neuron +  8) = ts.i();
    sc->m_input_ncp(first_neuron +  9) = ts.Omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron +  9));

    sc->m_input_ncp(first_neuron + 10) = ts.omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron + 10));

    sc->m_input_ncp(first_neuron + 11) = ts.f();
    AssertAngularRange(sc->m_input_ncp(first_neuron + 11));

    sc->m_input_ncp(first_neuron + 12) = ga->a()*ga_dis_factor;
    sc->m_input_ncp(first_neuron + 13) = ga->e();
    sc->m_input_ncp(first_neuron + 14) = ga->i();
    sc->m_input_ncp(first_neuron + 15) = ga->Omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron + 15));

    sc->m_input_ncp(first_neuron + 16) = ga->omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron + 16));

    sc->m_input_ncp(first_neuron + 17) = ga->f();
    AssertAngularRange(sc->m_input_ncp(first_neuron + 17));

  } else if (!delta_steering && m_ga_steering && !ga) {// d == 18

    sc->m_input_ncp(first_neuron     ) = ps.a()*ps_dis_factor;
    sc->m_input_ncp(first_neuron +  1) = ps.e();
    sc->m_input_ncp(first_neuron +  2) = ps.i();
    sc->m_input_ncp(first_neuron +  3) = ps.Omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron +  3));

    sc->m_input_ncp(first_neuron +  4) = ps.omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron +  4));

    sc->m_input_ncp(first_neuron +  5) = ps.f();
    AssertAngularRange(sc->m_input_ncp(first_neuron +  5));

    sc->m_input_ncp(first_neuron +  6) = ts.a()*ts_dis_factor;
    sc->m_input_ncp(first_neuron +  7) = ts.e();
    sc->m_input_ncp(first_neuron +  8) = ts.i();
    sc->m_input_ncp(first_neuron +  9) = ts.Omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron +  9));

    sc->m_input_ncp(first_neuron + 10) = ts.omega();
    AssertAngularRange(sc->m_input_ncp(first_neuron + 10));

    sc->m_input_ncp(first_neuron + 11) = ts.f();
    AssertAngularRange(sc->m_input_ncp(first_neuron + 11));

    sc->m_input_ncp(first_neuron + 12) = 0.0;
    sc->m_input_ncp(first_neuron + 13) = 0.0;
    sc->m_input_ncp(first_neuron + 14) = 0.0;
    sc->m_input_ncp(first_neuron + 15) = 0.0;
    sc->m_input_ncp(first_neuron + 16) = 0.0;
    sc->m_input_ncp(first_neuron + 17) = 0.0;
  } else  FatalError("CSCSteering::setNCInputKepler(): error in function input");

}; // end CSCSteering::setNCInputKepler


// ----------------------------------------------------------------------------
void CSCSteering::setNCInputPropellant(int neuron, CSC* sc)
// ----------------------------------------------------------------------------
// private
{
  sc->m_input_ncp(neuron) = sc->getPropellantMass()/sc->getInitialPropellantMass();

}; // end CSCSteering::setNCInputPropellant


// ----------------------------------------------------------------------------
void CSCSteering::setNCInput(CNC*          nc,
                             const CState& ps,
                             const CState& ts,
                             CSC*          sc,
                             const CState* ga)
// ----------------------------------------------------------------------------
// private
{
  const int input_dimension = nc->getNumberOfInputNeurons();
  const NCInputType nctype = nc->getInputType();

  int steering_input_dimension = (sc->isSolarSail()) ? (input_dimension) : (input_dimension - 1);

  int input_blocksize;

  // neurocontroller type is cartesian OR polar OR Keplerian elements
  // ----------------------------------------------------------------
  if (nctype == C_INPUT || nctype == P_INPUT || nctype == E_INPUT) {
    if      (steering_input_dimension == 6 && !m_ga_steering)
      input_blocksize = 6;
    else if (steering_input_dimension == 12 && !m_ga_steering)
      input_blocksize = 12;
    else if (steering_input_dimension == 12 &&  m_ga_steering)
      input_blocksize = 9;
    else if (steering_input_dimension == 18 &&  m_ga_steering)
      input_blocksize = 18;
    else
      FatalError("CSCSteering::performSteering(): internal error");
  }
  // neurocontroller type is cartesian/polar OR cartesian/Keplerian elements
  // OR polar/Keplerian elements
  // ----------------------------------------------------------------
  else if (nctype == CP_INPUT || nctype == CE_INPUT || nctype == PE_INPUT) {
    if      (steering_input_dimension == 12 && !m_ga_steering)
      input_blocksize = 6;
    else if (steering_input_dimension == 24 && !m_ga_steering)
      input_blocksize = 12;
    else if (steering_input_dimension == 24 &&  m_ga_steering)
      input_blocksize = 9;
    else if (steering_input_dimension == 36 &&  m_ga_steering)
      input_blocksize = 18;
    else
      FatalError("CSCSteering::performSteering(): internal error");
  }
  // neurocontroller includes all tpyes of coordinates
  // cartesian, polar and Keplerian elements
  // -----------------------------------------------------------------
  else if (nctype == CPE_INPUT) {
    if      (steering_input_dimension == 18 && !m_ga_steering)
      input_blocksize = 6;
    else if (steering_input_dimension == 36 && !m_ga_steering)
      input_blocksize = 12;
    else if (steering_input_dimension == 36 &&  m_ga_steering)
      input_blocksize = 9;
    else if (steering_input_dimension == 54 &&  m_ga_steering)
      input_blocksize = 18;
    else
      FatalError("CSCSteering::performSteering(): internal error");
  } else FatalError("CSCSteering::performSteering(): internal error");

  if (nctype == C_INPUT)
   // =================
  {
    setNCInputCartesian(0, steering_input_dimension - 1, ps, ts, sc);
  }
  else if (nctype == P_INPUT)
        // =================
  {
    setNCInputPolar(0, steering_input_dimension - 1, ps, ts, sc);
  }
  else if (nctype == E_INPUT)
        // =================
  {
    setNCInputKepler(0, steering_input_dimension - 1, ps, ts, sc);
  }
  else if (nctype == CP_INPUT)
        // ==================
  {
    setNCInputCartesian(0, input_blocksize - 1, ps, ts, sc);
    setNCInputPolar(input_blocksize, 2*input_blocksize - 1, ps, ts, sc);
  }
  else if (nctype == CE_INPUT)
        // ==================
  {
    setNCInputCartesian(0, input_blocksize - 1, ps, ts, sc);
    setNCInputKepler(input_blocksize, 2*input_blocksize - 1, ps, ts, sc);
  }
  else if (nctype == PE_INPUT)
        // ==================
  {
    setNCInputPolar(0, input_blocksize - 1, ps, ts, sc);
    setNCInputKepler(input_blocksize, 2*input_blocksize - 1, ps, ts, sc);
  }
  else if (nctype == CPE_INPUT)
        // ===================
  {
    setNCInputCartesian(0, input_blocksize - 1, ps, ts, sc);
    setNCInputPolar(input_blocksize, 2*input_blocksize - 1, ps, ts, sc);
    setNCInputKepler(2*input_blocksize, 3*input_blocksize - 1, ps, ts, sc);
  }
  else
    FatalError("CSCSteering::performSteering(): internal error");

  if (!sc->isSolarSail())
    setNCInputPropellant(input_dimension - 1, sc);

}; // end CSCSteering::setNCInput






// ----------------------------------------------------------------------------
void CSCSteering::performSteering(CSC*          sc,
                                  CNC*          nc,
                                  const CState& ps,
                                  const CState& ts,
                                  const CState* ga)
// ----------------------------------------------------------------------------
// private
{
  int first_lsl, last_lsl;

  double clock = 0.0, cone  = 0.0, xi_r  = 0.0, xi_t  = 0.0, xi_h  = 0.0;

  const int  input_dim = nc->getNumberOfInputNeurons();
  const int output_dim = nc->getNumberOfOutputNeurons();

  sc->setCtrlInputDim(input_dim);
  sc->setCtrlOutputDim(output_dim);

  const NCInputType nctype = nc->getInputType();

  int steering_output_dimension(output_dim);
  if (sc->isSolarSail() == false) steering_output_dimension -= 1;
  if (m_ga_steering     == true) steering_output_dimension  -= 2;

  setNCInput(nc, ps, ts, sc, ga);          // set neurocontroller input
  nc->propagate(sc->getCtrlInputVec(), sc->getCtrlOutputVec());

  // interpret neurocontroller output
  // ---------------------------------------------------
  if (m_stst_index == DEOPT || m_stst_index == DEMINMAX)
  // ---------------------------------------------------
  {
    if (m_stst_index == DEOPT)
    {
      first_lsl = static_cast<int>(aopt);
      last_lsl  = static_cast<int>(omegaopt);
    }
    else // if (m_stst_index == DEMINMAX)
    {
      first_lsl = static_cast<int>(amax);
      last_lsl  = static_cast<int>(omegamin);
    }

    // get requested direction over all local steering laws

    Vector nc_steeringpart(steering_output_dimension);
    int i;
    for (i = 0; i < steering_output_dimension; ++i)
    {
      nc_steeringpart(i) = sc->getCtrlOutput(i);//nc_output(i);
    }

    double sum_nc_steeringpart = Sum(nc_steeringpart);
    Vector nc_lslpart_norm = (sum_nc_steeringpart != 0.0) ? (nc_steeringpart / sum_nc_steeringpart) : Vector(steering_output_dimension);

    for (i = first_lsl; i <= last_lsl; ++i)
    {
      m_lslvec(i) = nc_lslpart_norm(i-first_lsl);

      (this->*lsl_ptr[i])(ps, ts);

      double xi_norm = sqrt(m_xi_r*m_xi_r*+m_xi_t*m_xi_t+m_xi_h*m_xi_h);
      if (xi_norm != 0.0)
      {
        xi_r += m_lslvec(i)*m_xi_r/xi_norm;
        xi_t += m_lslvec(i)*m_xi_t/xi_norm;
        xi_h += m_lslvec(i)*m_xi_h/xi_norm;
      }
    }

  }

  // --------------------------------
  else if (m_stst_index == DIRECTION)
  // --------------------------------
  {
    xi_r = 2.0*(sc->getCtrlOutput(0)-0.5); // scaling [0,1] => [-1,1]
    xi_t = 2.0*(sc->getCtrlOutput(1)-0.5); // scaling [0,1] => [-1,1]
    xi_h = 2.0*(sc->getCtrlOutput(2)-0.5); // scaling [0,1] => [-1,1]

    double xi_norm = sqrt(xi_r*xi_r+xi_t*xi_t+xi_h*xi_h);
    if (xi_norm != 0.0)
    {
      xi_r /= xi_norm;
      xi_t /= xi_norm;
      xi_h /= xi_norm;
    }
  }
  else
    FatalError("CSCSteering::performSteering(): wrong steering strategy");

  if (sc->isEP())   // set throttle
  {
    double throttle(0.0);
    if (m_ga_steering==true)
    {
      throttle = (std::max(std::min(sc->getCtrlOutput(output_dim-3),1.0),0.0));
    }
    else 
    {
      double t1(sc->getCtrlOutput(output_dim-1));
      double t2(sc->getCtrlOutput(output_dim-2));
      const double off1(0.83);
      const double fac1(0.6);
      throttle= std::max(std::min(fac1*(t1-t2+off1),1.0),0.0);
    }     
    sc->setThrottle(throttle);
  }
  sc->setAttitude(xi_r, xi_t, xi_h);  // set attitude
  if (!SpaceLib::AssertElevationRange(sc->getConeAngle()))
  {
    //cout << "cone: " << sc->getConeAngle() << endl;
    //cout << "xi_r: " << xi_r << "  xi_t: " << xi_t << "  xi_h: " << xi_h << endl;
  }
};

// ----------------------------------------------------------------------------
void CSCSteering::performSteering(CSC* sc, CNC* nc)
{
  double clock = 0.0, cone  = 0.0, xi_r  = 0.0, xi_t  = 0.0, xi_h  = 0.0;

  sc->setCtrlInputDim(nc->getNumberOfInputNeurons());
  sc->setCtrlOutputDim(nc->getNumberOfOutputNeurons());

  const NCInputType nctype = nc->getInputType();

  int steering_output_dimension(sc->getCtrlOutputDim());
  if (sc->isSolarSail()==false) --steering_output_dimension;
  if (m_ga_steering==true) steering_output_dimension -= 2;

  Vector ctrl_output(sc->getCtrlOutputDim()), ctrl_input(sc->getCtrlInputDim());

  for (int k = 0; k < sc->getCtrlInputDim(); ++k) 
  {
    ctrl_input(k) = sc->getCtrlInput(k);
  }

  nc->propagate(ctrl_input, ctrl_output); // apply controller

  for (int k = 0; k < sc->getCtrlOutputDim(); ++k)
  {
    sc->setCtrlOutput(k, ctrl_output(k));
  }

  double xi_norm = 0.0;
  // interpret neurocontroller output
  switch (m_stst_index)
  {
  case DIRECTION:
    xi_r = 2.0*(ctrl_output(0)-0.5); // scaling [0,1] => [-1,1]
    xi_t = 2.0*(ctrl_output(1)-0.5); // scaling [0,1] => [-1,1]
    xi_h = 2.0*(ctrl_output(2)-0.5); // scaling [0,1] => [-1,1]

    xi_norm = sqrt(xi_r*xi_r+xi_t*xi_t+xi_h*xi_h);
    if (xi_norm != 0.0)
    {
      xi_r /= xi_norm;
      xi_t /= xi_norm;
      xi_h /= xi_norm;
    }
    break;
  case DEOPT:
  case DEMINMAX:
  default:
    Warning("CSCSteering::performSteering(): wrong steering strategy");
    break;
  }

  if (sc->isEP()) // set throttle
  {
    double throttle(0.0);
    if (m_ga_steering==true)
    {
      throttle = (std::max(std::min(ctrl_output(sc->getCtrlOutputDim()-3),1.0),0.0));
    }
    else 
    {
      double t1(ctrl_output(sc->getCtrlOutputDim()-1)),t2(ctrl_output(sc->getCtrlOutputDim()-2));
      const double off1(0.83);
      const double fac1(0.6);
      throttle= std::max(std::min(fac1*(t1-t2+off1),1.0),0.0);
    }     
    sc->setThrottle(throttle);
  }

  sc->setAttitude(xi_r, xi_t, xi_h); // set s/c attitude
}

// ----------------------------------------------------------------------------
void CSCSteering::performSteering(CSC*  sc,
                                  const Vector& lsl_vector,
                                  const CState& ps,
                                  const CState& ts)
// ----------------------------------------------------------------------------
// private
{
  int first_lsl, last_lsl;

  double clock = 0.0,
         cone  = 0.0,
         xi_r  = 0.0,
         xi_t  = 0.0,
         xi_h  = 0.0;

  // ---------------------------------------------------
  if (m_stst_index == DEOPT || m_stst_index == DEMINMAX)
  // ---------------------------------------------------
  {
    if (m_stst_index == DEOPT)
    {
      first_lsl = static_cast<int>(aopt);
      last_lsl  = static_cast<int>(omegaopt);
    }
    else // if (m_stst_index == DEMINMAX)
    {
      first_lsl = static_cast<int>(amax);
      last_lsl  = static_cast<int>(omegamin);
    }

    // get requested direction over all local steering laws

    Vector lsl_vector_norm = lsl_vector / Sum(lsl_vector);

    for (int i = first_lsl; i <= last_lsl; i++)
    {
      m_lslvec(i) = lsl_vector_norm(i-first_lsl);

      (this->*lsl_ptr[i])(ps, ts);

      double xi_norm = sqrt(m_xi_r*m_xi_r+m_xi_t*m_xi_t+m_xi_h*m_xi_h);
      if (xi_norm != 0.0)
      {
        xi_r += m_lslvec(i)*m_xi_r/xi_norm;
        xi_t += m_lslvec(i)*m_xi_t/xi_norm;
        xi_h += m_lslvec(i)*m_xi_h/xi_norm;
      }
    }

    sc->setAttitude(xi_r, xi_t, xi_h);
  }

  // --------------------------------
  else if (m_stst_index == ALPHABETA)
  // --------------------------------
  {
    clock = lsl_vector(0);
    cone  = lsl_vector(1);

    sc->setAttitude(clock, cone);
  }

  else
    FatalError("CSCSteering::performSteering(): wrong steering strategy");

};

// ----------------------------------------------------------------------------
void CSCSteering::performSteering(CSC* sc,
                                  const double clock_angle,
                                  const double cone_angle)
// ----------------------------------------------------------------------------
// private
{
  sc->setAttitude(clock_angle, cone_angle);

}; // end CSCSteering::performSteering
