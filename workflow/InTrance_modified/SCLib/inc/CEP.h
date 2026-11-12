/*!
\file CEP.h
\brief A class that provides the functionality for an electric propulsion (EP) system.

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

  29.10.2006 (AO)

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
#ifndef CEP_H_INCLUDED
#define CEP_H_INCLUDED

#include <CPropSys.h>

/// A class that provides the functionality for an electric propulsion (EP) system.
class CEP : public CPropSys
{
                  /// throttle in the range [0.0,1.0]
       private:   double m_throttle;
                  /// Get the throttle of the EP system
                  /** \return throttle in the range [0.0,1.0] */
       protected: double getThrottle(void) const { return m_throttle; };
                  /// Set the throttle of the EP system
                  /** \param throttle throttle in the range [0.0,1.0] */
       protected: void setThrottle(const double throttle);

                  /// Throttle type
       private:   ThrottleType m_throttletype;
                  /// Get the information, whether the throttle type of the EP system is "variable" or "bang-bang"
                  /** \return throttle type */
       protected: ThrottleType getThrottleType(void) const { return m_throttletype; };
                  /// Set the throttle type of the EP system to "variable" or to "bang-bang"
                  /** \param type throttle type (BANGBANG or VARIABLE) */
       protected: void setThrottleType(const ThrottleType type) { m_throttletype = type; };

protected:

  // ----- pure virtual functions -----

  /// Get actual specific impulse of the EP system
  /** \return actual specific impulse [s] */
  virtual double getIsp(void) = 0;

	/// Get the actual mass flow through the EP system
  /** \return actual mass flow through the EP system [kg/s] */
  virtual double getMassFlow(void) = 0;

  /// Get the actual maximum thrust of the EP system
  /** \return actual maximum thrust [N] */
  virtual double getMaxThrust(void) = 0;

  /// Get the actual thrust of the EP system
  /** \return actual thrust of the EP system [N] */
  virtual double getThrust(void) = 0;

  /// Check whether the remaining propellant is sufficient for a thrust phase of length dt; if not, set propellant tank empty
  /** \param dt length of the thrust phase [d] */
  virtual void CheckRemainingPropellant(const double dt) = 0;

public:

  /// Constructor with given host state, i.e. a CSC object 
  /// (which is also a state object)
  CEP(const PropulsionSystemType psystype, CState* state = 0) : CPropSys(psystype,state)
  {
     m_throttle     = 1.0;
     m_throttletype = UNDEFINED_TT;
  }

  /// Copy constructor
  CEP(const CEP& ep, CState* state = 0) : CPropSys(ep, state)
  {
    m_throttle     = ep.m_throttle;
    m_throttletype = ep.m_throttletype;
  }

  /// Assign operator
  CEP & operator=(const CEP & ep)
  {
    CPropSys::operator =(ep);
    m_throttle     = ep.m_throttle;
    m_throttletype = ep.m_throttletype;
    return *this;
  };

  /// Destructor
  ~CEP(void) {};

  // ----- virtual functions -----

  virtual void printVrmlAppearance(std::ostream& os,const std::string& sc_name) 
    const throw ();

  /// Virtual function that gives error if not implemented in derived class
  virtual PropulsionSystemType getPropulsionSystemType(void) const
  {
    FatalError("class CEP: Function getPropulsionSystemType() not defined for the selected propulsion system");
    return NEP; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getNumberOfThrusters(void) const
  {
    FatalError("class CEP: Function getNumberOfThrusters() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getNumberOfUsedThrusters(void) const
  {
    FatalError("class CEP: Function getNumberOfUsedThrusters() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getPowerOfUsedThruster(void) const
  {
    FatalError("class CEP: Function getPowerOfUsedThruster() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getPower(void) const
  {
    FatalError("class CEP: Function getPower() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getMaxPower(void) const
  {
    FatalError("class CEP: Function getMaxPower() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  double getSolarRadiationDose(void) const { return 0.0; };
  void setSolarRadiationDose(const double srd) {};

}; // end class CEP

#endif // CEP_H_INCLUDED
