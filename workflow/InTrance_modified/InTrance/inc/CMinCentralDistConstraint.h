/*!
\file CMinCentralDistConstraint.h
\version 1.0
\brief Min. allowed distance from central body constraint.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 21.10.2006 (AO)

\b Purpose:
   Derived from CConstraint class, objects instantiated from
   CMinCentralDistConstraint class signal violation of a given allowed
   minimum distance from current central body. If the condition is fullfilled,
   e.g. the relative distance is lower than the allowed minimum distance, member
   method "isMatched(void)" returns false. As this constraint is of hard
   constraints type, usually the trajectory integration is terminated.

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

#ifndef CMINCENTRALDISTCONSTRAINT_H
#define CMINCENTRALDISTCONSTRAINT_H

#include <CConstraint.h>
#include <CSimulation.h>

namespace constraints
{

  /// constraint object of minimum allowed distance from central body
  class CMinCentralDistConstraint : virtual public CConstraint
  {
  public:
    CMinCentralDistConstraint(void) throw();
    CMinCentralDistConstraint(const CMinCentralDistConstraint& cs) throw();
    CMinCentralDistConstraint(CSimulation* sim) throw();
    ~CMinCentralDistConstraint(void) throw();
    virtual bool isMatched(void) throw();
  };

  //----------------------------------------------------------------------------
  inline CMinCentralDistConstraint::CMinCentralDistConstraint(void) throw()
    : CConstraint(true) {}

  //----------------------------------------------------------------------------
  inline CMinCentralDistConstraint::CMinCentralDistConstraint(const CMinCentralDistConstraint& cs) throw()
    : CConstraint(cs) {}

  //----------------------------------------------------------------------------
  inline CMinCentralDistConstraint::CMinCentralDistConstraint(CSimulation* sim)
    throw() : CConstraint(sim, true) {}

  //----------------------------------------------------------------------------
  inline CMinCentralDistConstraint::~CMinCentralDistConstraint(void) throw() {}

  //----------------------------------------------------------------------------
  inline bool CMinCentralDistConstraint::isMatched(void) throw()
  {
    if ( m_sim->getCurrentSC().CState::r() < m_sim->getMinSolarDistance() )
      return false;
    return true;
  } // bool isMatched(void) throw()
  
} // namespace constraints

#endif