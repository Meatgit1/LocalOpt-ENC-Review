/*!
\file CSEP.cpp
\brief A class that provides the functionality for a solar electric propulsion (SEP) system.

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

#include <CSC.h>

using namespace std;

#ifdef BORLAND
  template<class T> inline T max(const T lhs,T rhs) {return lhs<rhs?rhs:lhs;};
  template<class T> inline T min(const T lhs,T rhs) {return lhs>rhs?rhs:lhs;};
#endif BORLAND

//------------------------------------------------------------------------------
double CSEP::getPower(void) const
{
  const double rp = pow(m_state->r(), getPowerExponent()); // [AU^2]

  const double available_power = getCharacteristicPower() / rp; // [kW]

  const double useful_power = std::min(available_power, getNumberOfThrusters() * getMaxThrusterPower()); // [kW]

  return getThrottle() * useful_power; // [kW]
} // double getPower(void) const


//------------------------------------------------------------------------------
double CSEP::getMaxPower(void) const
{
  const double rp = pow(m_state->r(), getPowerExponent()); // [AU^2]

  const double available_power = getCharacteristicPower() / rp; // [kW]

  const double useful_power = std::min(available_power, getNumberOfThrusters() * getMaxThrusterPower()); // [kW]

  return useful_power; // [kW]
} // double getMaxPower(void) const


//------------------------------------------------------------------------------
void CSEP::CheckRemainingPropellant(const double dt)
{
  if (getMassFlow() * dt > getPropellantMass()) setEmpty(true);
} // void CheckRemainingPropellant(const double dt)

//------------------------------------------------------------------------------
void CSEP::ThrusterAllocation(const double p)
{
  m_nt = getNumberOfThrusters();        // start with all thrusters
  m_pt = 0.0;                           // start with zero power

  const double pmax = getMaxThrusterPower();
  const double popt = getOptThrusterPower();
  const double pmin = getMinThrusterPower();

  while (m_nt > 1.0 && m_pt < popt)     // true at the beginning
  {
    m_pt = p / m_nt;                    // calculate power for one thruster
    if (m_pt < popt) m_nt -= 1.0;       // if power is less than the power for optimum
                                        // Isp, then use one thruster less
  }                                     // now m_nt holds the number of running
                                        // thrusters (at least one)
  m_pt = p / m_nt;                      // and m_pt is the power PPU input for each
                                        // thruster
  if (m_pt > pmax)                      // however, if m_pt > pmax
  {
    m_nt += 1.0;                        // we have to use one thruster more
    m_pt = p / m_nt;                    // and m_pt is thus less
  }

  if (m_nt == 1.0)                      // if there is only one thruster left
  {
    if (p < pmin)                       // and we have less than pmin available
    {
      m_nt = 0.0;                       // then all thrusters are switched off
      m_pt = 0.0;
    }
  }
} // void ThrusterAllocation(const double p)


//------------------------------------------------------------------------------
double CSEP::getIsp(void) throw()
{
  if (!isEmpty() && (getMassFlow()>0.0))
  {
    return getThrust()/(getMassFlow()*G0);
  }
  else return 0.0;
} // double getIsp(void)


//------------------------------------------------------------------------------
void CSEP::EOM_RK_POLAR(double t, const Vector& y, Vector& yp, void* pAux)
{
  CSC* pSpacecraft = static_cast<CSC*>(pAux);

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

  yp(0) = dr;
  yp(1) = dphi;
  yp(2) = dtheta;
  yp(3) = ddr;
  yp(4) = ddphi;
  yp(5) = ddtheta;
  yp(6) = dmp;
  yp(7) = 0.0;

  pSpacecraft->IncrementDECounter(); // Increment function call count
} // void EOM_RK_POLAR(double t, const Vector& y, Vector& yp, void* pAux)

//------------------------------------------------------------------------------
void CSEP::printVrmlAppearance(ostream& os, const string& name) const throw ()
{

  os << "DEF " << name << " Transform\n";
  os << "{" << endl;
  os << " children\n";
  os << " [" << endl;
  os << "  DEF Capsule Transform\n";
  os << "  {" << endl;
  os << "   children\n";
  os << "   [" << endl;
  os << "    Shape\n";
  os << "    {" << endl;
  os << "     appearance Appearance\n";
  os << "     {" << endl;
  os << "      material Material\n";
  os << "      {" << endl;
  os << "       ambientIntensity 0.0667" << endl;
  os << "       diffuseColor 0.48 0.53 0.6" << endl;
  os << "       specularColor 0.32 0.4 0.4" << endl;
  os << "       shininess 0.54" << endl;
  os << "       emissiveColor 0.25 0.27 0.31\n";
  os << "      }\n";
  os << "     }" << endl;
  os << "     geometry Cone\n";
  os << "     {\n";    
  os << "      bottomRadius 0.05" << endl;
  os << "      height 0.1\n";
  os << "     }\n";
  os << "    }\n";
  os << "   ]" << endl;
  os << "   translation 0.0 0.05 0.0\n";
  os << "  }," << endl;
  os << "  DEF Tank Transform\n";
  os << "  {" << endl;
  os << "   children\n";
  os << "   [" << endl;
  os << "    Shape\n";
  os << "    {" << endl;
  os << "     appearance Appearance\n";
  os << "     {" << endl;
  os << "      material Material\n";
  os << "      {" << endl;
  os << "       ambientIntensity 0.0667" << endl;
  os << "       diffuseColor 0.48 0.53 0.6" << endl;
  os << "       specularColor 0.32 0.4 0.4" << endl;
  os << "       shininess 0.54" << endl;
  os << "       emissiveColor 0.25 0.27 0.31" << endl;
  os << "       transparency 0.4\n";
  os << "      }\n";
  os << "     }" << endl;
  os << "     geometry Cylinder\n";
  os << "     {" << endl;
  os << "      radius 0.05" << endl;
  os << "      height 0.1\n";
  os << "     }\n";
  os << "    }\n";
  os << "   ]" << endl;
  os << "   translation 0.0 -0.05 0.0" << endl;
  os << "   center 0.0 -0.1 0.0\n";
  os << "  }," << endl;
  os << "  DEF Fuel Transform\n";
  os << "  {" << endl;
  os << "   children\n";
  os << "   [" << endl;
  os << "    Shape\n";
  os << "    {" << endl;
  os << "     appearance Appearance\n";
  os << "     {" << endl;
  os << "      material Material\n";
  os << "      {" << endl;
  os << "       ambientIntensity .0467" << endl;
  os << "       diffuseColor .68 .68 .22" << endl;
  os << "       specularColor .39 .39 .27" << endl;
  os << "       shininess .19" << endl;
  os << "       emissiveColor .18 .18 .06\n";
  os << "      }\n";
  os << "     }" << endl;
  os << "     geometry Cylinder\n";
  os << "     {" << endl;
  os << "      radius 0.04" << endl;
  os << "      height 0.09\n";
  os << "     }\n";
  os << "    }\n";
  os << "   ]" << endl;
  os << "   translation 0.0 -0.05 0.0" << endl;
  os << "   center 0.0 0.045 0.0\n";
  os << "  }," << endl;


  os << "  DEF SAPanel Transform\n";
  os << "  {" << endl;
  os << "   children\n";
  os << "   [" << endl;
  os << "    DEF Stange Transform\n";
  os << "    {" << endl;
  os << "     children\n";
  os << "     [" << endl;
  os << "      Shape\n";
  os << "      {" << endl;
  os << "       appearance Appearance\n";
  os << "       {" << endl;
  os << "        material Material\n";
  os << "        {" << endl;
  os << "         ambientIntensity 0.0667" << endl;
  os << "         diffuseColor 0.48 0.53 0.6" << endl;
  os << "         specularColor 0.32 0.4 0.4" << endl;
  os << "         shininess 0.54" << endl;
  os << "         emissiveColor 0.25 0.27 0.31" << endl;
  os << "         transparency 0.4\n";
  os << "        }\n";
  os << "       }" << endl;
  os << "       geometry Cylinder\n";
  os << "       {" << endl;
  os << "        radius 0.005" << endl;
  os << "        height 0.15\n";
  os << "       }\n";
  os << "      }\n";
  os << "     ]\n";
  os << "    }" << endl;
  os << "    DEF SA1 Transform\n";
  os << "    {" << endl;
  os << "     children\n";
  os << "     [" << endl;
  os << "      Shape\n";
  os << "      {" << endl;
  os << "       appearance Appearance\n";
  os << "       {" << endl;
  os << "        material Material\n";
  os << "        {" << endl;
  os << "         ambientIntensity 0" << endl;
  os << "         diffuseColor 0.1 0.11 0.3" << endl;
  os << "         shininess 0.1" << endl;
  os << "         emissiveColor 0.1 0.11 0.3" << endl;
  os << "         specularColor 0.19 0.19 0.19\n";
  os << "        }\n";
  os << "       }" << endl;
  os << "       geometry Box\n";
  os << "       {" << endl;
  os << "        size 0.3 0.001 0.1\n";
  os << "       }\n";
  os << "      }\n";
  os << "     ]" << endl;
  os << "     translation 0.0 0.21 0.0" << endl;
  os << "     rotation 0.0 0.0 1.0 1.57\n";
  os << "    }" << endl;
  os << "    DEF SA2 Transform\n";
  os << "    {" << endl;
  os << "     children\n";
  os << "     [" << endl;
  os << "      Shape\n";
  os << "      {" << endl;
  os << "       appearance Appearance\n";
  os << "       {" << endl;
  os << "        material Material\n";
  os << "        {" << endl;
  os << "         ambientIntensity 0" << endl;
  os << "         diffuseColor 0.1 0.11 0.3" << endl;
  os << "         shininess 0.1" << endl;
  os << "         emissiveColor 0.1 0.11 0.3" << endl;
  os << "         specularColor 0.19 0.19 0.19\n";
  os << "        }\n";
  os << "       }" << endl;
  os << "       geometry Box\n";
  os << "       {" << endl;
  os << "        size 0.3 0.001 0.1\n";
  os << "       }\n";
  os << "      }\n";
  os << "     ]" << endl;
  os << "     translation 0.0 -0.21 0.0" << endl;
  os << "     rotation 0.0 0.0 1.0 1.57\n";
  os << "    }\n";
  os << "   ]" << endl;
  os << "   translation 0.0 -0.05 0.0" << endl;
  os << "   rotation 1.0 0.0 0.0 1.57\n";
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
  os << "   translation 0.0 0.15 0.0" << endl;
  os << "   center 0.0 -0.05 0.0\n";
  os << "  }\n";
  os << " ]\n";
  os << "}" << endl << endl;

} // void printVrmlAppearance(std::ostream& os) const throw ();