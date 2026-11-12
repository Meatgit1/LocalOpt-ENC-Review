/*!
\file CProximityConstraint.h
\version 1.0
\brief Proximity constraint object for trajectory intergration.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 21.10.2006 (AO)

\b Purpose:
   Derived from class CConstraint this soft constraint is to signal whether
   a given maximum relative distance and maximum relative velocity w.r.t. a 
   given target body is not exceeded. If the used CProximity object
   tells, that those limits, which are stored within the CSimulation object, are
   exceeded the constraint's method "isMatched(void)" returns true and false
   if the limits are kept. In application InTrance this objects signal the
   occurence of a boundary criterion (target relative distance/velocity) match.

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

#ifndef CPROXIMITYCONSTRAINT_H_INCLUDED
#define CPROXIMITYCONSTRAINT_H_INCLUDED

#include <CConstraint.h>
#include <CSimulation.h>

namespace constraints
{

  // declaration

  /// proximity constraint
  class CProximityConstraint : virtual public CConstraint 
  { 
  public:
    CProximityConstraint(void) throw();
    CProximityConstraint(CSimulation* sim) throw();
    CProximityConstraint(const CProximityConstraint& constraint) throw();
    virtual bool isMatched(void) throw();
    virtual ~CProximityConstraint(void) throw();
  };

  // implementation

  //----------------------------------------------------------------------------
  inline CProximityConstraint::CProximityConstraint() throw()
    : CConstraint() {};

  //----------------------------------------------------------------------------
  inline CProximityConstraint::CProximityConstraint(CSimulation* sim) throw()
    : CConstraint(sim) {};

  //----------------------------------------------------------------------------
  inline CProximityConstraint::CProximityConstraint(
    const CProximityConstraint& constraint) throw() : CConstraint(constraint) {}

  //----------------------------------------------------------------------------
  inline bool CProximityConstraint::isMatched(void) throw()
  {
    return (m_sim->getCurrProximity().CheckProximity(
      m_sim->getMaxRelDistance(), m_sim->getMaxRelVelocity() ) == false);
  } // bool isMatched(void) throw()

  //----------------------------------------------------------------------------
  inline CProximityConstraint::~CProximityConstraint(void) throw() {};

} // namespace constraints

#endif