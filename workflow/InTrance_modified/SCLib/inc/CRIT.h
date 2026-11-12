/*!
\file CRIT.h
\brief A class that provides the functionality for a RIT system.

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

  28.10.2006 (AO)

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
#ifndef RIT_H_INCLUDED
#define RIT_H_INCLUDED

#include "CSEP.h"
#include <new>

const double DFLT_MAX_THRUSTER_POWER = 6.209;
const double DFLT_OPT_THRUSTER_POWER = 6.209;
const double DFLT_MIN_THRUSTER_POWER = 0.142;
const double DFLT_RIT_THRUST_COEFF0 = 142.0;
const double DFLT_RIT_THRUST_COEFF1 = 34.67E3;
const double DFLT_RIT_M_FLOW_COEFF0 = 0.0207E-3;
const double DFLT_RIT_M_FLOW_COEFF1 = 0.0967E-6;

enum RIT_COEFF
{
   MAX_THRUSTER_POWER // 0
  ,OPT_THRUSTER_POWER // 1
  ,MIN_THRUSTER_POWER // 2
  ,THRUST_0           // 3
  ,THRUST_1           // 4
  ,M_FLOW_0           // 5
  ,M_FLOW_1           // 6
  ,N_RIT_COEFFICIENTS // 7
};

/// A class that provides the functionality for a RIT system.
class CRIT : public CSEP
{
private:

  double m_coefficients[N_RIT_COEFFICIENTS];

  void initDefaults(void) throw();
  double getMassFlow(void);
  double getMaxThrust(void);
  double getMaxThrusterPower(void) const throw();
  double getMinThrusterPower(void) const throw();
  double getOptThrusterPower(void) const throw();
  double getThrust(void);

public:

  /// Default constructor
  CRIT(CState* state = 0) throw();

  /// Copy constructor
  CRIT(const CRIT & rit,CState* state = 0) throw();

  /// assign operator
  const CRIT& operator=(const CRIT& rit) throw();

  /// Get propulsion system type
  PropulsionSystemType getPropulsionSystemType(void) const { return RIT; };

  /// set new propulsion system coeffient
  void setCoefficient(const RIT_COEFF which ,const double value) throw();

  /// current propulsion system setting
  double getCoefficient(const RIT_COEFF which) const throw();

  /// update members depending on parameter settings (should be called after
  /// changing the propulsion system parameters)
  void update(void) throw();

  /*!
   \brief make copy of object and set given host object pointer

    Because polymorphism is used in implementation of different propulsion
    systems a clone() - method is needed. Due to the fact that only the object
    itself knows its real type (except you try a dynamic_cast<..> with all
    possible object types and handle possibly thrown exceptions) it makes sense
    to have such a clone function which returns a base class pointer to the
    new object which is a copy of itself.

   \param state pointer to "host" object
   \return CPropSys Pointer to new CRIT object
  */
  virtual CPropSys* clone(CState* state) throw();

  /// Destructor
  ~CRIT(void) {};

}; // end class CRIT

// implementation

//------------------------------------------------------------------------------
inline CPropSys* CRIT::clone(CState* state) throw()
{
  try 
  {
    return (new CRIT((*this), state));
  } 
  catch (std::bad_alloc ) 
  {
    return 0;
  }
} // virtual CPropSys* clone(CState* state) throw();

//------------------------------------------------------------------------------
inline void CRIT::initDefaults(void) throw()
{
  setCoefficient(MAX_THRUSTER_POWER, DFLT_MAX_THRUSTER_POWER);
  setCoefficient(OPT_THRUSTER_POWER, DFLT_OPT_THRUSTER_POWER);
  setCoefficient(MIN_THRUSTER_POWER, DFLT_MIN_THRUSTER_POWER);
  setCoefficient(THRUST_0, DFLT_RIT_THRUST_COEFF0);
  setCoefficient(THRUST_1, DFLT_RIT_THRUST_COEFF1);
  setCoefficient(M_FLOW_0, DFLT_RIT_M_FLOW_COEFF0);
  setCoefficient(M_FLOW_1, DFLT_RIT_M_FLOW_COEFF1);
} // void CRIT::initDefaults(void) throw()

//------------------------------------------------------------------------------
inline CRIT::CRIT(CState* state) throw() : CSEP(RIT,state)
{
  initDefaults();
} // CRIT(CState*)

//------------------------------------------------------------------------------
inline CRIT::CRIT(const CRIT & rit,CState* state) throw() : CSEP(rit,state)
{
  for (unsigned int i=0;i<N_RIT_COEFFICIENTS;++i)
  {
    m_coefficients[i] = rit.m_coefficients[i];
  } // for (unsigned int i=0;i<N_RIT_COEFFICIENTS;++i)  
} // CRIT(const CRIT & rit,CState* state = 0)

//------------------------------------------------------------------------------
inline const CRIT& CRIT::operator=(const CRIT& rit) throw()
{
  CSEP::operator=(rit);
  for (unsigned int i=0;i<N_RIT_COEFFICIENTS;++i)
  {
    m_coefficients[i] = rit.m_coefficients[i];
  } // for (unsigned int i=0;i<N_RIT_COEFFICIENTS;++i)
  return (*this);
} // const CRIT& operator=(const CRIT& rit) throw();

//------------------------------------------------------------------------------
inline void CRIT::setCoefficient(const RIT_COEFF which, const double value ) 
throw()
{
  m_coefficients[which] = value;
} // void setCoefficient(const double) throw()

//------------------------------------------------------------------------------
inline double CRIT::getCoefficient(const RIT_COEFF which) const throw()
{
  return m_coefficients[which];
} // double getCoefficient(void) const throw()

//------------------------------------------------------------------------------
inline double CRIT::getMaxThrusterPower(void) const throw()
{
  return m_coefficients[MAX_THRUSTER_POWER];
} // double getMaxThrusterPower(void) const throw()

//------------------------------------------------------------------------------
inline double CRIT::getMinThrusterPower(void) const throw()
{
  return m_coefficients[MIN_THRUSTER_POWER];
} // double getMinThrusterPower(void) const throw()

//------------------------------------------------------------------------------
inline double CRIT::getOptThrusterPower(void) const throw()
{
  return m_coefficients[OPT_THRUSTER_POWER];
} // double getOptThrusterPower(void) const throw()

//------------------------------------------------------------------------------
inline void CRIT::update(void) throw()
{
  getThrust();
} // void update(void) throw()

#endif // CRIT_H_INCLUDED