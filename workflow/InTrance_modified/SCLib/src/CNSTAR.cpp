/*!
\file CNSTAR.cpp
\brief A class that provides the functionality for a NSTAR system.

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

#define TCOEFF0 -3.4318
#define TCOEFF1 37.365
#define MCOEFF0  0.74343
#define MCOEFF1  0.20951
#define MCOEFF2  0.25205

using namespace std;

double CNSTAR::getMassFlow(void)
{
  ThrusterAllocation(getPower());

  return (isEmpty()) ? 0.0 : (getNumberOfUsedThrusters() * (MCOEFF0 + MCOEFF1 * getPowerOfUsedThruster() + MCOEFF2 * getPowerOfUsedThruster() * getPowerOfUsedThruster()) / 1e6); // [kg/s]
};

double CNSTAR::getMaxThrust(void)
{
  ThrusterAllocation(getMaxPower());

  return (isEmpty()) ? 0.0 : (getNumberOfUsedThrusters() * (TCOEFF0 + TCOEFF1 * getPowerOfUsedThruster()) / 1e3); // [N]
};

double CNSTAR::getThrust(void)
{
  ThrusterAllocation(getPower());

  return (isEmpty()) ? 0.0 : (getNumberOfUsedThrusters() * (TCOEFF0 + TCOEFF1 * getPowerOfUsedThruster()) / 1e3); // [N]
};