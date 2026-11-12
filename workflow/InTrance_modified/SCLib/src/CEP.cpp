/*!
\file CEP.cpp
\brief A class that provides the functionality for an electric propulsion (EP) system.

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

  07 Jan 2005

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

using namespace std;

//------------------------------------------------------------------------------
void CEP::setThrottle(const double throttle)
{
#ifdef development_version
  if(throttle < 0.0 || throttle > 1.0) FatalError ("CSEP::setThrottle(): throttle must be in the range (0.0,1.0)");
#endif
  switch (getThrottleType())
  {
  case VARIABLE:
    m_throttle = throttle;
    break;
  case BANGBANG:
    m_throttle = (throttle > 0.5) ? 1.0 : 0.0;
    break;
  default:
    FatalError("CEP::setThrottle(): undefined throttle type");
    break;
  } // switch (getThrottleType())
} // void CEP::setThrottle(const double throttle)

//------------------------------------------------------------------------------
void CEP::printVrmlAppearance(ostream& os, const string& sc_name) const throw ()
{
    os << "DEF " << sc_name << " Transform\n";
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

} // void printVrmlAppearance(ostream& os) const throw ()