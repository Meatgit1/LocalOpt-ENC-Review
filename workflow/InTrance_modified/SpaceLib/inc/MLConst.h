/*!
\file MLConst.h
\version 1.7
\brief Mathematical constants and conversion factors.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 23.06.2005 (AO)

\b Purpose:
   This file provides different mathematical constants and conversion factors.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.

\b Copyright:
   (c) 2001-2004 Bernd Dachwald, Andreas Ohndorf
*/

#ifndef MATHCONST_H
#define MATHCONST_H

#include <limits>

namespace MathLib
{
  const double EPS_MACH = std::numeric_limits<double>::epsilon();
               ///< Difference between 1 and the smallest value greater than 1 that is representable.
  const double EPS = 10.0 * EPS_MACH;
               ///< Ten times the difference between 1 and the smallest value greater than 1 that is representable.

  // mathematical constants
  const double PI = 3.1415926535897932384626433832795; ///< The value of pi.
  const double HALFPI  = 0.5 * PI;                     ///< The value of pi/2.
  const double TWOPI   = 2.0 * PI;                     ///< The value of 2*pi.
  const double THREEPI = 3.0 * PI;                     ///< The value of 3*pi.

  // conversion factors for angular calculations
  const double DEG2RAD = PI / 180.0;          ///< Conversion factor: degrees -> radians.
  const double RAD2DEG = 1.0 / DEG2RAD;       ///< Conversion factor: degrees <- radians.
  const double HOUR2DEG = 15.0;               ///< Conversion factor: arc hours -> degrees.
  const double DEG2HOUR = 1.0 / HOUR2DEG;     ///< Conversion factor: arc hours <- degrees.
  const double RAD2HOUR = RAD2DEG * DEG2HOUR; ///< Conversion factor: radians -> arc hours.
  const double HOUR2RAD = 1.0 / RAD2HOUR;     ///< Conversion factor: radians <- arc hours.
  const double HOUR2MIN = 60.0;               ///< Conversion factor: hours -> minutes.
  const double MIN2HOUR = 1.0 / HOUR2MIN;     ///< Conversion factor: hours <- minutes.
  const double MIN2SEC = HOUR2MIN;            ///< Conversion factor: minutes -> seconds.
  const double SEC2MIN = MIN2HOUR;            ///< Conversion factor: minutes <- seconds.
  const double HOUR2SEC = HOUR2MIN * MIN2SEC; ///< Conversion factor: hours -> seconds.
  const double SEC2HOUR = 1.0 / HOUR2SEC;     ///< Conversion factor: hours <- seconds.
  const double ARCMIN2RAD = DEG2RAD / 60.0;   ///< Conversion factor: arc minutes -> radians.
  const double RAD2ARCMIN = 60.0 / DEG2RAD;   ///< Conversion factor: arc minutes <- radians.
  const double ARCSEC2RAD = DEG2RAD / 3600.0; ///< Conversion factor: arc seconds -> radians.
  const double RAD2ARCSEC = 3600.0 / DEG2RAD; ///< Conversion factor: arc seconds <- radians.

} // end of namespace MathLib

#endif // MATHCONST_H
