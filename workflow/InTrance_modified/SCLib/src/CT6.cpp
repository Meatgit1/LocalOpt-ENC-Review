/*!
\file CT6.cpp
\brief A class that provides the functionality for a QinetiQ T6 system.

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

#define ISP_OPT                         4600.0
#define REQ_KILOWATTS_PER_NEWTON_THRUST   30.1

using namespace std;

double CT6::getMassFlow(void)
{
  ThrusterAllocation(getPower());

  return (isEmpty()) ? 0.0 : (getNumberOfUsedThrusters() * getThrust() / ISP_OPT / G0); // [kg/s]
};

double CT6::getMaxThrust(void)
{
  ThrusterAllocation(getMaxPower());

  return (isEmpty()) ? 0.0 : (getNumberOfUsedThrusters() * (getPowerOfUsedThruster() - getMinThrusterPower()) / REQ_KILOWATTS_PER_NEWTON_THRUST); // [N]
};

double CT6::getThrust(void)
{
  ThrusterAllocation(getPower());

  return (isEmpty()) ? 0.0 : (getNumberOfUsedThrusters() * (getPowerOfUsedThruster() - getMinThrusterPower()) / REQ_KILOWATTS_PER_NEWTON_THRUST); // [N]
};