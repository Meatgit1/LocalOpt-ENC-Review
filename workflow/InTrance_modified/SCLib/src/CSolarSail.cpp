/*!
\file CSolarSail.cpp
\brief A class that provides the functionality for a solar sail.

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

 30.10.2006 (AO)

\b Copyright:

  (c) 2001-2006 by Bernd Dachwald

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/
#include <iostream>
#include <iomanip>

//SCLib
#include <configSC.h>
#include <SpaceConst.h>
#include <CSC.h>
#include <NativeBridgeProbe.h>

#include <CBody.h>

#include <errno.h>

using namespace std;
using namespace MathLib;

void CSolarSail::setOpticalModel(const bool   ideal_or_real,
                                 const double srd /* = 0.0 */,
                                 const double max_temperature /* = HUGE_VAL */,
                                 const double rho_0 /* = 0.0 */,
                                 const double s_0 /* = 0.0 */,
                                 const double epsilonf_0 /* = 0.0 */,
                                 const double epsilonb_0 /* = 0.0 */,
                                 const double rho_infty /* = 0.0 */,
                                 const double s_infty /* = 0.0 */,
                                 const double epsilonf_infty /* = 0.0 */,
                                 const double epsilonb_infty /* = 0.0 */,
                                 const double rho_half /* = 0.0 */,
                                 const double s_half /* = 0.0 */,
                                 const double epsilonf_half /* = 0.0 */,
                                 const double epsilonb_half /* = 0.0 */,
                                 const double Bf /* = 0.0 */,
                                 const double Bb /* = 0.0 */)
{
  m_ideal_or_real_sail = ideal_or_real;

  if(isRealSail())
  {
    m_srd             = srd;
    m_max_temperature = max_temperature;
    m_rho_0           = rho_0;
    m_s_0             = s_0;
    m_epsilonf_0      = epsilonf_0;
    m_epsilonb_0      = epsilonb_0;
          m_rho_infty       = rho_infty;
    m_s_infty         = s_infty;
    m_epsilonf_infty  = epsilonf_infty;
    m_epsilonb_infty  = epsilonb_infty;
          m_rho_half        = rho_half;
    m_s_half          = s_half;
    m_epsilonf_half   = epsilonf_half;
    m_epsilonb_half   = epsilonb_half;
    m_rho_lambda      = (rho_half      == 0.0) ? 0.0 : log(2.0)/rho_half;
    m_s_lambda        = (s_half        == 0.0) ? 0.0 : log(2.0)/s_half;
    m_epsilonf_lambda = (epsilonf_half == 0.0) ? 0.0 : log(2.0)/epsilonf_half;
    m_epsilonb_lambda = (epsilonb_half == 0.0) ? 0.0 : log(2.0)/epsilonb_half;
    m_Bf              = Bf;
    m_Bb              = Bb;
  }
};

double CSolarSail::getTemperature(const double r, const double coneangle) const
/// input must be in AUD units, result is also in AUD units!!!
{
  if(isRealSail())
  {
    const double term1 = (1.0-getrho())/(getepsilonf()+getepsilonb());
    const double term2 = S0_AUD/(SIGMA_AUD*r*r)*cos(coneangle);
    return pow(term1*term2, 0.25);
  }
  else
    return 0.0;
};

double CSolarSail::getSolarEnergyFlux(const double r, const double coneangle) const
/// input must be in AUD units, result is also in AUD units!!!
{
  return cos(coneangle)/(r*r);
};

void CSolarSail::EOM_RK_POLAR(double t
                        , const MathLib::Vector& y
                        , MathLib::Vector& yp
                        , void * pAux)
{
  
  CSC * pSpacecraft = static_cast<CSC *>(pAux);

  // Opt-in diagnostics only; the original force and derivative expressions below are unchanged.
  NativeBridgeProbe::check_schedule_stage(t,y);

  const double r         = y(0);
  const double phi       = y(1);
  const double theta     = y(2);
  const double sintheta  = sin(theta);
  const double costheta  = cos(theta);
  const double tantheta  = tan(theta);
  const double dr        = y(3);
  const double dphi      = y(4);
  const double dtheta    = y(5);
  const double lambda    = getLightness();
  const double zeta      = atan2(dtheta,dphi*costheta);
  const double alpha     = pSpacecraft->getClockAngle();
  const double beta      = pSpacecraft->getConeAngle();
  const double sinalpha  = sin(alpha+zeta);
  const double cosalpha  = cos(alpha+zeta);
  const double sinbeta   = sin(beta);
  const double cosbeta   = cos(beta);

  double factor    = lambda*GM_SUN_AUD/(r*r);

  double Psi_r, Psi_phi, Psi_theta;

  if (isIdealSail())
  {
    Psi_r     =          cosbeta*cosbeta*cosbeta;
    Psi_phi   = cosalpha*sinbeta*cosbeta*cosbeta;
    Psi_theta = sinalpha*sinbeta*cosbeta*cosbeta;
  }
  else
  {
    const double Hsinbeta  = getH()*sinbeta;
    const double GcosbetaK = getG()*cosbeta + getK();

    const double beta_hat = beta - atan2(Hsinbeta, GcosbetaK);
    const double Q        = sqrt(GcosbetaK*GcosbetaK + Hsinbeta*Hsinbeta) * cosbeta;
    const double Q_prime  = Q / (getG0()+getK0());

    Psi_r     = Q_prime         *cos(beta_hat);
    Psi_phi   = Q_prime*cosalpha*sin(beta_hat);
    Psi_theta = Q_prime*sinalpha*sin(beta_hat);
  }
#define SAIL_JETTISON
#ifdef SAIL_JETTISON
  if (r > 5.0) factor = 0.0;
#endif SAIL_JETTISON

  const double ddr     =   r*dtheta*dtheta
                         + r*dphi*dphi*costheta*costheta
                         - GM_SUN_AUD/(r*r) + factor*Psi_r;

  const double ddphi   = - 2.0*dr*dphi/r
                         + 2.0*dphi*dtheta*tantheta
                         + factor*Psi_phi/(r*costheta);

  const double ddtheta = - 2.0*dr*dtheta/r
                         - dphi*dphi*sintheta*costheta
                         + factor*Psi_theta/r;

  const double dsrd    = getSolarEnergyFlux(r, beta);

  yp = Vector(dr, dphi, dtheta, ddr, ddphi, ddtheta, 0.0, dsrd);

  if (NativeBridgeProbe::enabled())
    NativeBridgeProbe::record_sail_rhs(t,y,yp,*pSpacecraft,
      factor*Psi_r,factor*Psi_phi,factor*Psi_theta,zeta);

  pSpacecraft->IncrementDECounter(); // Increment function call count
};

#ifndef SLOW_CODE
Vector e_vec_r(3);
Vector e_vec_h(3);
Vector e_vec_t(3);
Vector e_vec_l(3);
Vector e_vec_n(3);
Vector h(3);
Vector r_vec_sun(3);
Vector s_vec_c(6);
Vector a_vec_thrust(3);
#endif // SLOW_CODE
 
//------------------------------------------------------------------------------
void CSolarSail::EOM_RK_CARTESIAN(double t,const Vector& y,Vector& yp,void* pAux)
{
  // if all works as designed all conditional compilation statements and comments
  // shall be deleted
  using namespace MathLib;
  CSC* sc = static_cast<CSC *>(pAux);
  CSolarSystem* ss;
  ss = sc->CBody::getRefSolSys();

  const double gm_sun   = ss->GM_AUd(SUN);
  const double gm       = (sc->CState::getCenterState())->getGM(AUD);
  const double lambda   = getLightness();
  const double alpha    = sc->getClockAngle();
  const double beta     = sc->getConeAngle();
  const double sinalpha = sin(alpha);
  const double cosalpha = cos(alpha);
  const double sinbeta  = sin(beta);
  const double cosbeta  = cos(beta);  

  const double r_sc = sqrt(pow(y(0),int(2))+pow(y(1),int(2))+pow(y(2),int(2)));
  // orbit frame unit vector in radial direction
#ifndef SLOW_CODE
  int i;
  for (i=0;i<3;++i) e_vec_r(i) = y(i)/r_sc;
#else
  const Vector e_vec_r(y(0)/r_sc, y(1)/r_sc, y(2)/r_sc);
#endif
  //const Vector e_vec_r_test(y(0)/r_sc, y(1)/r_sc, y(2)/r_sc);
  //assert(Norm(e_vec_r-e_vec_r_test)==0);

  // specific orbit angular momentum vector h = r x v
#ifndef SLOW_CODE
  h(0) = y(1)*y(5)-y(2)*y(4);
  h(1) = y(2)*y(3)-y(0)*y(5);
  h(2) = y(0)*y(4)-y(1)*y(3);
#else
  const Vector h(Cross(y.Slice(0,2), y.Slice(3,5)));
#endif
  //const Vector h_test(Cross(y.Slice(0,2), y.Slice(3,5)));
  //assert(Norm(h-h_test)==0);

  // orbit frame unit vector in angular momentum vector direction
#ifndef SLOW_CODE
  const double h_norm(Norm(h));
  for (i=0;i<3;++i) e_vec_h(i) = (h(i)/h_norm);
#else
  const Vector e_vec_h(h/Norm(h));
#endif
  //const Vector e_vec_h_test(h/Norm(h));
  //assert(Norm(e_vec_h-e_vec_h_test)==0);

  /// orbit frame unit vector in tranversal direction
#ifndef SLOW_CODE
  e_vec_t(0) = e_vec_h(1)*e_vec_r(2)-e_vec_h(2)*e_vec_r(1);
  e_vec_t(1) = e_vec_h(2)*e_vec_r(0)-e_vec_h(0)*e_vec_r(2);
  e_vec_t(2) = e_vec_h(0)*e_vec_r(1)-e_vec_h(1)*e_vec_r(0);
#else
  const Vector e_vec_t(Cross(e_vec_h, e_vec_r));
#endif
  //const Vector e_vec_t_test(Cross(e_vec_h, e_vec_r));
  //assert(Norm(e_vec_t-e_vec_t_test)==0);

  const double a_grav = (-gm/pow(r_sc,int(2)));       // central body grav. acc.

#ifndef SLOW_CODE
  s_vec_c = ss->Cart(sc->getCenterBodyId(),SUN,t,AUD);
  for (i=0;i<3;++i) r_vec_sun(i) = s_vec_c(i)+y(i);
#else
  const Vector s_vec_c(ss->Cart(sc->getCenterBodyId(),SUN,t,AUD));
  const Vector r_vec_sun(s_vec_c(0)+y(0),s_vec_c(1)+y(1),s_vec_c(2)+y(2));
#endif
  //const Vector s_vec_c_test(ss->Cart(sc->getCenterBodyId(),SUN,t,AUD));
  //const Vector r_vec_sun_test(s_vec_c_test(0)+y(0),s_vec_c_test(1)+y(1),s_vec_c_test(2)+y(2));
  //assert(Norm(r_vec_sun-r_vec_sun_test)==0);

  const double r_sun_square = 
    pow(r_vec_sun(0),int(2))+pow(r_vec_sun(1),int(2))+pow(r_vec_sun(2),int(2));
  const double r_sun = sqrt(r_sun_square);                  // distance from sun

  // unit vector of sunlight direction
#ifndef SLOW_CODE
  for (i=0;i<3;++i) e_vec_l(i) = (r_vec_sun(i)/r_sun);
#else
  const Vector e_vec_l(r_vec_sun(0)/r_sun, r_vec_sun(1)/r_sun, r_vec_sun(2)/r_sun);
#endif
  //const Vector e_vec_l_test(r_vec_sun(0)/r_sun, r_vec_sun(1)/r_sun, r_vec_sun(2)/r_sun);
  //assert(Norm(e_vec_l-e_vec_l_test) == 0);

#ifndef SLOW_CODE
  double fac[3];
  fac[0] = cosbeta; fac[1] = sinbeta*cosalpha; fac[2] = sinbeta*sinalpha; 
  for (i=0;i<3;++i) 
    e_vec_n(i) = fac[0]*e_vec_r(i) + fac[1]*e_vec_t(i) + fac[2]*e_vec_h(i);
#else
  const Vector e_vec_n(
    cosbeta*e_vec_r + sinbeta*cosalpha*e_vec_t + sinbeta*sinalpha*e_vec_h);
#endif
  //const Vector e_vec_n_test(
  //  cosbeta*e_vec_r + sinbeta*cosalpha*e_vec_t + sinbeta*sinalpha*e_vec_h);
  //assert(Norm(e_vec_n-e_vec_n_test)==0);

  const double cos_incident_angle = bounded(Dot(e_vec_n, e_vec_l),0.0,1.0);
  const double incident_angle = acos(cos_incident_angle);

  double a_thrust = 0;
  a_vec_thrust.Set(0.0);
  if (isIdealSail())
  {
    a_thrust = lambda*(gm_sun/r_sun_square)*pow(cos_incident_angle,int(2));
    if (a_thrust>0)
    {
      a_vec_thrust = e_vec_n;
      a_vec_thrust *= a_thrust;
    }
  }
  else
  {
    assert(0==1); // throw error to signal incomplete implementation
  }

  // now we can set the state derivative vector
  unsigned int j;                        
  for (j = 0;j<3;++j) 
  {
    yp(j)   = y(j+3);                    // first, derivative of position vector
    yp(j+3) = (a_grav*e_vec_r(j)) + a_vec_thrust(j);  // second, velocity vector
  }

  yp(6) = 0.0;                          // massflow (none because of solar sail)
  yp(7) = getSolarEnergyFlux(r_sun, incident_angle);

  sc->IncrementDECounter(); // Increment function call count
} // void EOM_RK_CARTESIAN(...)


//------------------------------------------------------------------------------
void CSolarSail::printVrmlAppearance(std::ostream& os,const string& name) 
const throw ()
{
  os << "DEF " << name << " Transform\n";
  os << "{" << endl;
  os << " children\n";
  os << " [" << endl;
  os << "  Shape\n";
  os << "  {" << endl;
  os << "   appearance Appearance\n";
  os << "   {" << endl;
  os << "    material Material\n";
  os << "    {" << endl;
  os << "     ambientIntensity 0.0667" << endl;
  os << "     diffuseColor 0.48 0.53 0.6" << endl;
  os << "     specularColor 0.32 0.4 0.4" << endl;
  os << "     shininess 0.54" << endl;
  os << "     emissiveColor 0.25 0.27 0.31" << endl;
  os << "     transparency 0.4\n";
  os << "    }\n";
  os << "   }" << endl;
  os << "   geometry Cylinder\n";
  os << "   {" << endl;
  os << "    radius 0.1" << endl;
  os << "    height 0.001\n";
  os << "   }\n";
  os << "  }," << endl;


  os << "  DEF Arrow Transform\n";
  os << "  {" << endl;
  os << "   children\n";
  os << "   [" << endl;
  os << "    Shape\n";
  os << "    {" << endl;
  os << "     appearance Appearance\n";
  os << "     {" << endl;
  os << "      material Material\n";
  os << "      {" << endl;
  os << "       emissiveColor 1.0 0.0 0.0\n";
  os << "      }\n";
  os << "     }" << endl;
  os << "     geometry Cylinder\n";
  os << "     {" << endl;
  os << "      radius 0.01" << endl;
  os << "      height 0.1\n";
  os << "     }\n";
  os << "    }," << endl;
  os << "    Transform\n";
  os << "    {" << endl;
  os << "     children\n";
  os << "     [" << endl;
  os << "      Shape\n";
  os << "      {" << endl;
  os << "       appearance Appearance\n";
  os << "       {" << endl;
  os << "        material Material\n";
  os << "        {" << endl;
  os << "         emissiveColor 1.0 0.0 0.0\n";
  os << "        }\n";
  os << "       }" << endl;
  os << "       geometry Cone\n";
  os << "       {" << endl;
  os << "        bottomRadius 0.02" << endl;
  os << "        height 0.05\n";
  os << "       }\n";
  os << "      }\n";
  os << "     ]" << endl;
  os << "     translation 0.0 0.075 0.0\n";
  os << "    }\n";
  os << "   ]" << endl;
  os << "   translation 0.0 0.055 0.0" << endl;
  os << "   center 0.0 -0.05 0.0\n";
  os << "  }\n";
  os << " ]\n";
  os << "}" << endl << endl;
} // virtual void printVrmlAppearance(std::ostream& os) const throw ();
