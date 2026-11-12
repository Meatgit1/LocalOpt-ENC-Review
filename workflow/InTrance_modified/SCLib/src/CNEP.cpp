/*!
\file CNEP.cpp
\brief A class that provides the functionality for a nuclear electric propulsion (NEP) system.

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

  23 Feb 2005 (AO)

\b Compiler:

  Microsoft Visual C++ .NET

\b Copyright:

  (c) 2001-2005 by Bernd Dachwald

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
#include <CSC.h>
#include <CEP.h>

using namespace std;

void CNEP::CheckRemainingPropellant(const double dt)
{
  if (getMassFlow() * dt > getPropellantMass()) setEmpty(true);
};

double CNEP::getThrust(void)
{
  return (isEmpty()) ? 0.0 : getThrottle() * getMaxThrust(); // [N]
};

double CNEP::getMassFlow(void)
{
  return (isEmpty()) ? 0.0 : getThrottle() * getMaxThrust() / (G0 * getIsp()); // [kg/s]
};

//------------------------------------------------------------------------------
void CNEP::EOM_RK_POLAR(double t,	const Vector& y, Vector& yp, void* pAux)
{
  CSC* pSpacecraft = static_cast<CSC*>(pAux);
  // Schedule-only fail-closed domain/path checks before native RHS arithmetic.
  // All force expressions below remain unchanged.
  const double mjd = pSpacecraft->MJD();

  const double r         = y(0);
  const double phi       = y(1);
  const double theta     = y(2);
  const double sintheta  = sin(theta);
  const double costheta  = cos(theta);
  const double tantheta  = tan(theta);
  const double dr        = y(3);
  const double dphi      = y(4);
  const double dtheta    = y(5);
  const double mp        = y(6);
  const double zeta      = atan2(dtheta,dphi*costheta);
  const double alpha     = pSpacecraft->getClockAngle();
  const double beta      = pSpacecraft->getConeAngle();
  const double sinalpha  = sin(alpha+zeta);
  const double cosalpha  = cos(alpha+zeta);
  const double sinbeta   = sin(beta);
  const double cosbeta   = cos(beta);
  const double m         = pSpacecraft->getSCMass()+mp;
  const double thrust    = pSpacecraft->getThrust() * THRUST_SI2AUD;
  const double flow      = pSpacecraft->getMassFlow() / SEC2DAY;

  const double Psi_r     =          cosbeta;
  const double Psi_phi   = cosalpha*sinbeta;
  const double Psi_theta = sinalpha*sinbeta;

  const double acc       = thrust/m;

  const double ddr     = r*dtheta*dtheta
                         + r*dphi*dphi*costheta*costheta
                         - GM_SUN_AUD/(r*r)
                         + acc*Psi_r;

  const double ddphi   = - 2.0*dr*dphi/r
                         + 2.0*dphi*dtheta*tantheta
                         + acc*Psi_phi/(r*costheta);

  const double ddtheta = - 2.0*dr*dtheta/r
                         - dphi*dphi*sintheta*costheta
                         + acc*Psi_theta/r;

  const double dmp     = - flow;

  yp = Vector(dr, dphi, dtheta, ddr, ddphi, ddtheta, dmp, 0.0);

  pSpacecraft->IncrementDECounter(); // Increment function call count
}

//------------------------------------------------------------------------------
void CNEP::EOM_RK_CARTESIAN(
                                    double t
                                  , const MathLib::Vector& y
                                  , MathLib::Vector & yp
                                  , void * pAux)
{
  CSC* sc = static_cast<CSC *>(pAux);
  CSolarSystem* ss;
  ss = sc->CBody::getRefSolSys();

  const double gm_sun   = ss->GM_AUd(SUN);
  const double gm       = (sc->CState::getCenterState())->getGM(AUD);
  const double alpha    = sc->getClockAngle();
  const double beta     = sc->getConeAngle();
  const double sinalpha = sin(alpha);
  const double cosalpha = cos(alpha);
  const double sinbeta  = sin(beta);
  const double cosbeta  = cos(beta);
  const double mp       = y(6);
  const double m        = sc->getSCMass()+mp;
  const double thrust   = sc->getThrust() * THRUST_SI2AUD;
  const double flow     = sc->getMassFlow() / SEC2DAY;

  const double r_sc = sqrt(pow(y(0),int(2))+pow(y(1),int(2))+pow(y(2),int(2)));
  // orbit frame unit vector in radial direction
  const MathLib::Vector e_vec_r(y(0)/r_sc, y(1)/r_sc, y(2)/r_sc);
  
  // specific orbit angular momentum vector
  const MathLib::Vector h(Cross(y.Slice(0,2), y.Slice(3,5)));

  // orbit frame unit vector in angular momentum vector direction
  const MathLib::Vector e_vec_h(h/Norm(h));

  /// orbit frame unit vector in tranversal direction
  const MathLib::Vector e_vec_t(Cross(e_vec_h, e_vec_r));

  const double a_grav = (-gm/pow(r_sc,int(2)));       // central body grav. acc.

  double a_thrust = thrust/m;
  MathLib::Vector a_vec_thrust(cosbeta,sinbeta*cosalpha,sinbeta*sinalpha);
  a_vec_thrust*=a_thrust;
  a_vec_thrust = 
      Dot(a_vec_thrust, e_vec_r)*e_vec_r
    + Dot(a_vec_thrust, e_vec_t)*e_vec_t
    + Dot(a_vec_thrust, e_vec_h)*e_vec_h;

  // now we can set the state derivative vector
  unsigned int i;                        
  for (i = 0;i<3;++i) 
  {
    yp(i)   = y(i+3);                    // first, derivative of position vector
    yp(i+3) = (a_grav*e_vec_r(i)) + a_vec_thrust(i);  // second, velocity vector
  }
  yp(6) = - flow;                                                    // massflow
  yp(7) = 0.0;
  
  sc->IncrementDECounter(); // Increment function call count
} // void EOM_RK_CARTESIAN(...)
