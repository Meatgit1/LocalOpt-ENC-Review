/*!
\file SCConst.h
\brief A header that provides constants and type declarations for the SCLib.

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

  15.11.2004

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
#ifndef SCCONST_H_INCLUDED
#define SCCONST_H_INCLUDED

#include <SpaceConst.h>
#include <SLTime.h>

using namespace SpaceLib;

/// Spacecraft propulsion system type
/** Characterizes a spacecraft propulsion system. */
enum PropulsionSystemType
{
  UNDEFINED_PST, ///< undefined propulsion system type
  SOLARSAIL,     ///< solar sail
  NEP,           ///< nuclear electric propulsion system
  NSTAR,         ///< NASA's NSTAR ion propulsion system
  T6,            ///< QinetiQ's T6 ion propulsion system
  RIT,           ///< RIT ion propulsion system
};

/// Spacecraft propulsion system throttle type
/** Characterizes the behavior/control of the propulsion system. */
enum ThrottleType
{
  UNDEFINED_TT, ///< undefined throttle type
  BANGBANG,     ///< throttle is either on (1.0) or off (0.0)
  VARIABLE      ///< throttle can continuously vary in the range [0.0,1.0]
};

/**
* \brief converts thrust(acceleration) unit SI=>AUD
*
* Converts the unit of thrust (mass specific, thus an acceleration) from SI
* units to units in AUD system.
* SI units (m*s^-2)) => AUD system units (AU*d^-2); d = day
*/
const double THRUST_SI2AUD = DAY2SEC * DAY2SEC / AU;

/** Converts the unit of thrust (mass specific, thus an acceleration)from AUD
  *  units to units in AUD system. AUD system units (AU*d^-2) => SI units (m*s^-2))
  */
const double THRUST_AUD2SI = AU / DAY2SEC / DAY2SEC;  ///< convert thrust(acceleration) unit AUD=>SI

/** Converts the unit of massflow from SI system to AUD system.
  *   (kg*s^-1) => (kg*d^-1))
  */
const double MASSFLOW_SI2AUD = DAY2SEC; ///< convert massflow unit SI => AUD

/** Converts the unit of massflow from AUD system to SI system.
  *   (kg*d^-1) => (kg*s^-1))
  */
const double MASSFLOW_AUD2SI = SEC2DAY; ///< convert massflow unit AUD => SI

#endif // SCCONST_H_INCLUDED
