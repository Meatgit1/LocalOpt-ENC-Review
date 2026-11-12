/*!
\file SpaceLib.h
\version 1.7
\brief The header file that you can use to include the whole SpaceLib.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 09.05.2006 (AO)

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.

\b Copyright:
   (c) 2001-2006 Bernd Dachwald, Andreas Ohndorf
*/

#ifndef SPACELIB_H
#define SPACELIB_H

// standard C runtime library
#include <assert.h>

// C++ Standard Template Library (STL)
#include <string>

// own function declarations
#include <Kepler.h>
#include <Euler.h>
#include <CAttitude.h>
#include <CLagrange.h>
#include <CProximity.h>
#include <MathLib.h>
#include <CPhysUnit.h>


/** \namespace SpaceLib
    \brief The namespace of the SpaceLib
    The whole functionality of the SpaceLib library is encapsuled in the namespace SpaceLib.
    You don't have to use the 'using' statement, because this is done at the end of this
    header file. So, be aware of the possibilty of compiler errors saying for example
    'multiple declaration'. If this happens, use the name of the namespace with the scope
    operator. */

using namespace Global;
using namespace MathLib;
using namespace SpaceLib;

namespace SpaceLib
{

  ///*! \brief periapsis distance Rp = a(1-e)
  // *
  // * \param a semi-major axis of orbit
  // * \param e eccentricity
  // * \return periapsis distance in unit system of parameter a
  // */
  //template<class T, class U, class V> const T CalcPeriapsisDistance(const U a, const V e)
  //{
  //  assert( ((a<0.0) && (e>1.0)) || ((a>=0.0) && (e<=1.0)) );
  //  return T(a) * ( T(1) - T(e) );
  //}

  ///*! \brief apoapsis distance Rp = a(1+e)
  // *
  // * \param a semi-major axis of orbit (>=0.0)
  // * \param e eccentricity (<=1.0)
  // * \return apoapsis distance in unit system of parameter a
  // */
  //template<class T, class U, class V> const T CalcApoapsisDistance(const U a, const V e)
  //{
  //  assert( (a>=0.0) && (e<=1.0) );
  //  return T(a) * ( T(1) + T(e) );
  //}

  ///*! \brief periapsis velocity Vp = sqrt( GM*( (2/Rp) - (1/|a|) ) )
  // *
  // * \param a semi-major axis of orbit
  // * \param e eccentricity
  // * \param GM gravitational parameter of central body
  // * \return Periapsis velocity in units depending on input parameter units
  // */
  //template<class T, class U, class V, class W> 
  //  const T CalcPeriapsisVelocity(const U a, const V e, const W GM)
  //{
  //  return sqrt(T(GM)*(T(2.0)/CalcPeriapsisDistance<T>(a,e)-T(1.0)/T(abs(a))));
  //}

  ///*! \brief periapsis velocity Vp = sqrt( GM*( (2/Rp) - (1/|a|) ) )
  // *
  // * \param Rp periapsis distance
  // * \param Ra apoapsis distance
  // * \param GM gravitational parameter of central body
  // * \return Periapsis velocity in units depending on input parameter units
  // */
  //template<class T, class U, class W>
  //  const T CalcPeriapsisVelocity2(const U Rp, const U Ra, const W GM)
  //{
  //  return sqrt(T(2.0)*T(GM)*T(Ra)/(T(Rp)*(T(Ra)+T(Rp))));
  //}


  ///*! \brief mean motion = sqrt(GM/(a^3))
  // * \param a semimajor axis
  // * \param GM gravitational parameter of central body
  // * \return mean motion (unit depending on GM)
  // */
  //template<class T,  class U>
  //const T CalcMeanMotion(const T a, const U GM)
  //{
  //  return sqrt(T(GM)/pow(a,int(3)));
  //}
}

#endif // SPACELIB_H
