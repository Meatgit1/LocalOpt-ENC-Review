/*!
\file CSEP.h
\brief A class that provides the functionality for a solar electric propulsion (SEP) system.

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
#ifndef CSEP_H_INCLUDED
#define CSEP_H_INCLUDED

class CSC;

/// A class that provides the functionality for a solar electric propulsion (SEP) system.
class CSEP : public CEP
{
  friend class CSC;

             /// Characteristic power [kW] of the solar arrays (at Earth distance)
  private:   double m_characteristic_power;
             /// Get characteristic power [kW] of the solar arrays
  public:    double getCharacteristicPower(void) const { return m_characteristic_power; };
             /// Set characteristic power [kW] of the solar arrays
  protected: void setCharacteristicPower(const double cpwr) { m_characteristic_power = cpwr; };

             /// Number of thrusters
  private:   double m_thrusters;
             /// Get number of thrusters
  public:    double getNumberOfThrusters(void) const { return m_thrusters; };
             /// Set Number of thrusters
  protected: void setNumberOfThrusters(const double n) { m_thrusters = n; };

             /// Variation exponent p of the solar array power with solar distance => P(r)/P0 ~ (1/r)^p
  private:   double m_power_exponent;
             /// Get variation exponent of the solar array power with solar distance
  public:    double getPowerExponent(void) const { return m_power_exponent; };
             /// Get variation exponent of the solar array power with solar distance
  protected: void setPowerExponent(const double p) { m_power_exponent = p; };

             /// Number of used thrusters (set by ThrusterAllocation())
  private:   double m_nt;
             /// Get number of used thruster
  public:    double getNumberOfUsedThrusters(void) const { return m_nt; };

             /// Available power for each used thruster (set by ThrusterAllocation())
  private:   double m_pt;
             /// Get available power for each used thruster
  public:    double getPowerOfUsedThruster(void) const { return m_pt; };

public:

  /// Equations of motion for solar electric propulsion system in RK-form
  void EOM_RK_POLAR(double t, const Vector& y, Vector& yp, void* pAux);

  virtual void printVrmlAppearance(std::ostream& os, const std::string& name)
    const throw ();

protected:

  /// Constructor
  CSEP(const PropulsionSystemType psystype, CState* state = 0) : CEP(psystype, state)
  {
    m_characteristic_power = 0.0;
    m_thrusters            = 0.0;
    m_power_exponent       = 2.0;
    m_nt                   = 0.0;
    m_pt                   = 0.0;
  };

  /// Copy constructor
  CSEP(const CSEP & sep, CState* state = 0) : CEP(sep,state)
  {
    m_characteristic_power = sep.m_characteristic_power;
    m_thrusters            = sep.m_thrusters;
    m_power_exponent       = sep.m_power_exponent;
    m_nt                   = sep.m_nt;
    m_pt                   = sep.m_pt;
  };

  // Assign operator
  CSEP& operator=(const CSEP & sep)
  {
    CEP::operator =(sep);
    m_characteristic_power = sep.m_characteristic_power;
    m_thrusters            = sep.m_thrusters;
    m_power_exponent       = sep.m_power_exponent;
    m_nt                   = sep.m_nt;
    m_pt                   = sep.m_pt;
    return *this;
  };

  /// Destructor
  ~CSEP(void) {};

  /// Allocate the optimal thruster combination for the specified power input
  /** \param p power input */
  void ThrusterAllocation(const double p);

  /// Get actual Isp of the SEP system (dependent on thrust and mass flow)
  double getIsp(void);

  /// Get actual power of the SEP system (dependent on solar distance and throttle)
  double getPower(void) const;

  /// Get actual maximum power of the SEP system (dependent on solar distance)
  double getMaxPower(void) const;

  /// Check whether propellant tank is empty
  void CheckRemainingPropellant(const double dt);

  // ----- virtual functions -----

  /// Virtual function that gives error if not implemented in derived class
  virtual PropulsionSystemType getPropulsionSystemType(void) const
  {
    FatalError("class CSEP: Function getPropulsionSystemType() not defined for the selected propulsion system");
    return RIT; // dummy
  }

  // ----- pure virtual functions -----

  /// Get the maximum power [kW] per thruster
  virtual double getMaxThrusterPower(void) const = 0;

  /// Get the minimum power [kW] per thruster
  virtual double getMinThrusterPower(void) const = 0;

  /// Get the optimum power [kW] per thruster (max Isp!)
  virtual double getOptThrusterPower(void) const = 0;

}; // end class CSEP


#endif // CSEP_H_INCLUDED