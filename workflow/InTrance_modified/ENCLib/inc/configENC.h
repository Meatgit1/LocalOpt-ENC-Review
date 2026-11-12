/*!
\file configENC.h
\brief ENCLib compilation configuration

\version

  1.7

\author

  Bernd Dachwald (AO)\n
  German Aerospace Center (DLR)\n
  Institute of Space Simulation\n
  Linder Hoehe\n
  51147 Cologne\n
  Germany\n
  bernd.dachwald@dlr.de

\author

  Andreas Ohndorf (AO)\n
  German Aerospace Center (DLR)\n
  Institute of Space Simulation\n
  Linder Hoehe\n
  51147 Cologne\n
  Germany\n
  andreas.ohndorf@gmx.de

\date

  09.05.2006 (AO)

\b Compiler:

  GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++

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
#ifndef CONFIGENC_H_INCLUDED
#define CONFIGENC_H_INCLUDED

#if defined _MSC_VER    // Microsoft C++ Compiler
 #pragma warning (disable: 4786) // "identifier was truncated to '255' characters in the browser information"
 #pragma warning (disable: 4250) // "inherits ??? via dominance"
#endif

#if defined _DEBUG
#define DEBUG 1
#endif

/*!
 * Tbd.
 */
typedef void (*KeyboardHandler) (const int first_key, const int second_key);

const double MIN_FITNESS(-2.0e10);  // minimal allowed fitness value

#endif