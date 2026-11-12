/*!
\file CRIT.cpp
\brief A class that provides the functionality for a RIT system.

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

//SCLib
#include <configSC.h>
#include <CSC.h>

using namespace std;

//------------------------------------------------------------------------------
double CRIT::getMassFlow(void)
{
  ThrusterAllocation(getPower());

//  return (isEmpty()) ? 0.0 : (getNumberOfUsedThrusters() 
//    * (getCoefficient(M_FLOW_0) * getThrust()/(getNumberOfUsedThrusters()==0.0?1.0:getNumberOfUsedThrusters()) + getCoefficient(M_FLOW_1))); // [kg/s]
  //return (isEmpty()) ? 0.0 : (getNumberOfUsedThrusters() 
  //  * (getCoefficient(M_FLOW_0) * getThrust()/(getNumberOfUsedThrusters()==0.0?1.0:getNumberOfUsedThrusters()) + getCoefficient(M_FLOW_1)));

  if (isEmpty() == false)
  {
    return getNumberOfUsedThrusters()*(getCoefficient(M_FLOW_0)*getPowerOfUsedThruster()+getCoefficient(M_FLOW_1));  // [kg/s]
  } else return 0.0;

} // double CRIT::getMassFlow(void)

//------------------------------------------------------------------------------
double CRIT::getMaxThrust(void)
{
  ThrusterAllocation(getMaxPower());

  //return (isEmpty()) ? 0.0 : (getNumberOfUsedThrusters() * ((getPowerOfUsedThruster()*1e3 - getCoefficient(THRUST_0)) / getCoefficient(THRUST_1))); // [N]
  return (isEmpty()) ? 0.0 : (getNumberOfUsedThrusters() * ((getPowerOfUsedThruster() - getCoefficient(THRUST_0)) * getCoefficient(THRUST_1))); // [N]
} // double CRIT::getMaxThrust(void)

//------------------------------------------------------------------------------
double CRIT::getThrust(void)
{
  ThrusterAllocation(getPower());

//  return (isEmpty()) ? 0.0 : (getNumberOfUsedThrusters() * ((getPowerOfUsedThruster()*1e3 - getCoefficient(THRUST_0)) / getCoefficient(THRUST_1))); // [N]
  return isEmpty() ? 0.0 : (getNumberOfUsedThrusters() * ((getPowerOfUsedThruster() - getCoefficient(THRUST_0)) * getCoefficient(THRUST_1))); // [N]
} // double CRIT::getThrust(void)
