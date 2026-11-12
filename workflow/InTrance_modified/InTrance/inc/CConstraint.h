/*!
\file CConstraint.h
\version 1.0
\brief Simulation constraint for trajectory intergration.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 21.10.2006 (AO)

\b Purpose:
   The evaluation integration of CSimulation objects is usually forced to obey
   certain constraints. An obvious one is, for example, not to exceed a given
   maximum flight time. Another one may be not reach certain celestial
   bodies, e.g. the sun, because of thermal conditions. This can be assured by
   specifying an allowed minimum distance to the mentioned celestial body.
   In order to support such general constraints in a flexible, adaptable and
   extensible general design, the concept of constraint classes is introduced.
   Realised by class hierarchy having the virtual base class CConstraint
   it allows the treament of many kinds of simulation restrictions in a general
   way. The integration of new, not yet known, constrainst is possible through
   usage of the OOP-concepts of polymorphism, derivation and overriding. This
   means you only have to create a new constraint class by deriving it from 
   CConstraint at least and to overide the virtual member method 
   "isMatched(void)" which should return true if the simulation shall continue
   and false to break the integration loop.

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

#ifndef CCONSTRAINT_H_INCLUDED
#define CCONSTRAINT_H_INCLUDED

class CSimulation;

/// namespace of all classes and functions related to simulation integration
/// constraints
namespace constraints
{

  // declaration

  /// constraints base class (abstract class)
  class CConstraint 
  {
  private:
    bool m_hard_constraint;                    ///< constraint type member field

  protected:
    CSimulation* m_sim;            ///< pointer to simulation object of interest

  public: 
    /// standard constructor
    CConstraint(const bool ishard = false) throw();
    /// simulation object constructor
    CConstraint(CSimulation* sim, const bool hard = false) throw();
    /// copy constructor
    CConstraint(const CConstraint& constraint) throw();
    /// assign operator
    CConstraint& operator=(const CConstraint& constraint) throw ();
    /// equaloperator
    bool operator==(const CConstraint& constraint) const throw();
    /// "lower/less than" - operator
    const bool operator<(const CConstraint& constraint) const throw(); 
    /// destructor
    virtual ~CConstraint(void) throw();

    // modifier/access methods

    /// get currently used reference simulation object
    CSimulation* getSimPtr(void) const throw();

    /// set new reference simulation object
    void setSimPtr(CSimulation*) throw();

    /*! \brief is constraint still obeyed?
     * This method is not implemented and makes this class abstract. So
     * , in order to use constraints you have to derive classes from this
     * base constraint class and implement an "isMatched"-method.
    \return true, if constraint is fulfilled (no break), and false if not (break)
    */
    virtual bool isMatched(void) = 0;

    /*! \brief tells whether the constraint is a hard one or not
    * Misobeying of a hard criterion in context of individual evaluations leads 
    * to minimum fitness values. Individual not obeying soft constraints will
    * in general be given a fitness value while not obeying hard constraints
    * result in invalid solutions, e.g. a solution leading to a sun-hitting
    * trajectory is absolutely not desired and sorted out because of the minimum
    * fitness value it receives.
    * \return true, if it is a hard constraint, and false if not*/
    bool isHardConstraint(void) const;

    /*! \brief set type of constraint (hard or not)
     * \param hard true to set constraint to be a hard one, false otherwise
    */
    void setHardConstraint(const bool hard) throw();
  }; // class CConstraint


  // implementation

  //-----------------------------------------------------------------------------
  inline CConstraint::CConstraint(const bool ishard) throw()
  { 
    setSimPtr(0);
    setHardConstraint(ishard);
  } // CConstraint(void) throw()

  //-----------------------------------------------------------------------------
  inline CConstraint::CConstraint(CSimulation* sim, const bool hard) throw() 
  { 
    setSimPtr(sim); 
    setHardConstraint(hard);
  } // CConstraint(CSimulation* sim, const bool hard) throw() 

  //-----------------------------------------------------------------------------
  inline CConstraint::CConstraint(const CConstraint& constraint) throw()
  {
    setSimPtr(constraint.getSimPtr());
    setHardConstraint(constraint.isHardConstraint());
  } // CConstraint(const CConstraint& constraint) throw()

  //-----------------------------------------------------------------------------
  inline 
    CConstraint& CConstraint::operator=(const CConstraint& constraint) throw()
  {
    setSimPtr(constraint.getSimPtr());
    setHardConstraint(constraint.isHardConstraint());
    return *this;
  }

  //-----------------------------------------------------------------------------
  inline 
    bool CConstraint::operator==(const CConstraint& constraint) const throw()
  {
    if (getSimPtr() == constraint.getSimPtr())
      if (isHardConstraint() == constraint.isHardConstraint()) return true;
    else return false;
  }

  //-----------------------------------------------------------------------------
  inline 
   const bool CConstraint::operator<(const CConstraint& constraint) const throw()
  {
    return false;
  }

  //-----------------------------------------------------------------------------
  inline CConstraint::~CConstraint(void) throw()  {};

  //----------------------------------------------------------------------------
  inline CSimulation* CConstraint::getSimPtr(void) const throw() 
  { 
    return m_sim; 
  } // CSimulation* getSimPtr(void) const throw() 

  //----------------------------------------------------------------------------
  inline void CConstraint::setSimPtr(CSimulation* sim) throw() { m_sim = sim;};

  //----------------------------------------------------------------------------
  inline bool CConstraint::isHardConstraint(void) const throw()
  {
    return m_hard_constraint;
  } // bool isHardConstraint(void) const throw()

  //----------------------------------------------------------------------------
  inline void CConstraint::setHardConstraint(const bool hard) throw()
  {
    m_hard_constraint = hard;
  } // void setHardConstraint(const bool hard) throw()

} // namespace constraints
#endif CCONSTRAINED_H_INCLUDED
