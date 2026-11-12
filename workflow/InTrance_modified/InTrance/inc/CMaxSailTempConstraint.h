/*!
\file CMaxSailTempConstraint.h
\version 1.0
\brief Maximum sail temperature constraint.

\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 21.10.2006 (AO)

\b Purpose:
   This is a "hard" constraint derived from CConstraint class for Solar Sail
   purposes and breaks trajectory integration if the maximum sail temperature 
   is exceeded. Current sail temperature and allowed maximum sail temperature
   are obtained from current simulation spacecraft.

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

#ifndef CMAXSAILTEMPCONSTRAINT_H_INCLUDED
#define CMAXSAILTEMPCONSTRAINT_H_INCLUDED

#include <CConstraint.h>
#include <CSimulation.h>

namespace constraints
{

  // declaration

  /// constraint object for max. allowed solar sail temperature
  class CMaxSailTempConstraint : virtual public CConstraint
  {
  public:
    CMaxSailTempConstraint(void);
    CMaxSailTempConstraint(CSimulation* sim);
    virtual ~CMaxSailTempConstraint(void) throw();
    virtual bool isMatched(void) throw();
  };

  // implementation

  //----------------------------------------------------------------------------
  inline CMaxSailTempConstraint::CMaxSailTempConstraint(void)
    : CConstraint(true)  {}

  //----------------------------------------------------------------------------
  inline CMaxSailTempConstraint::CMaxSailTempConstraint(CSimulation* sim)
    : CConstraint(sim, true)  {}

  //----------------------------------------------------------------------------
  inline bool CMaxSailTempConstraint::isMatched(void) throw()
  {
    if ( m_sim->getCurrentSC().getTemperature() > 
      m_sim->getCurrentSC().getMaxTemperature() )
      return false;
    return true;
  } // bool isMatched(void) throw()

  //----------------------------------------------------------------------------
  inline CMaxSailTempConstraint::~CMaxSailTempConstraint(void) throw() {};

} // namespace constraints

#endif // CMAXSAILTEMPCONSTRAINT_H_INCLUDED