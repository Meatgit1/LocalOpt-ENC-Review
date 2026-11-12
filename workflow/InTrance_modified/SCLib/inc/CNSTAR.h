/*!
\file CNSTAR.h
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

\author

  Andreas Ohndorf\n
  German Aerospace Center (DLR)\n
  Institute of Space Simulation\n
  Linder Hoehe\n
  51147 Cologne\n
  Germany\n
  andreas.ohndorf@gmx.de

\date

  12.01.2006 (AO)

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
#ifndef CNSTAR_H_INCLUDED
#define CNSTAR_H_INCLUDED

#include <new>
#include "CSEP.h"

/// A class that provides the functionality for a NSTAR system.
class CNSTAR : public CSEP
{
  double getMassFlow(void);
  double getMaxThrust(void);
  double getMaxThrusterPower(void) const { return 2.0; };
  double getMinThrusterPower(void) const { return 0.5; };
  double getOptThrusterPower(void) const { return 1.83374; };
  double getThrust(void);

public:

  /// Default constructor
  CNSTAR(CState* state = 0) : CSEP(NSTAR, state) {};

  /// Copy constructor
  CNSTAR(const CNSTAR & nstar,CState* state = 0) : CSEP(nstar, state) {};

  /// Get propulsion system type
  PropulsionSystemType getPropulsionSystemType(void) const { return NSTAR; };

  /*!
   \brief make copy of object and set given host object pointer

    Because polymorphism is used in implementation of different propulsion
    systems a clone() - method is needed. Due to the fact that only the object
    itself knows its real type (except you try a dynamic_cast<..> with all
    possible object types and handle possibly thrown exceptions) it makes sense
    to have such a clone function which returns a base class pointer to the
    new object which is a copy of itself.

   \param state pointer to "host" object
   \return CPropSys Pointer to new CNSTAR object
  */
  virtual CPropSys* clone(CState* state)
  {
    try 
    {
      return (new CNSTAR(*this, state));
    } 
    catch (std::bad_alloc ) 
    {
      return 0;
    }
  };

  ~CNSTAR(void) {};

}; // end class CNSTAR

#endif // CNSTAR_H_INCLUDED
