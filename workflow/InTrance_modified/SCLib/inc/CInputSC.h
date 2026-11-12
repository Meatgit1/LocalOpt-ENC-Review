/*!
\file CInputSC.h
\brief An auxiliarary class that contains the values of the spacecraft input file.

\version

  1.7

\author

  Bernd Dachwald (BD)\n
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
  andreas.ohndorf@dlr.de


\date

  30.06.2005 (AO)

\b Compiler:

  GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++

\b Copyright:

  (c) 2001-2004 by Bernd Dachwald, Andreas Ohndorf

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/
#ifndef CINPUTSPACECRAFT_H_INCLUDED
#define CINPUTSPACECRAFT_H_INCLUDED

// STL
#include <string>

// SpaceLib
#include <Global.h>
#include <MLVector.h>

// SCLib
#include <SCConst.h>

const std::string DFLT_SC_NAME = "SC";

/// spacecraft configuration describing object
class CInputSC
{
//=============================================================================
// CInputSC CLASS INTERFACE
//=============================================================================
public:

  //-----------------------------------------------------------------

  PropulsionSystemType  psystype;
                        // spacecraft propulsion system

  std::string name; ///< name of spacecraft

  /// s/c's moments of inertia about principle body axis
  MathLib::Vector    moments_of_inertia;

  /// maximum moments about principle body axes for steering purposes
  MathLib::Vector    max_steering_moments;

  // ONLY FOR SOLAR SAIL

  double             sailcraft_performance;
                     // sailcraft lightness number

  bool               ideal_or_real_sail;
                     // true,  if ideal sail
                     // false, if real  sail

  double             max_temperature;
                     // maximum tolerable temerature of the sail film (only used for real sails)

  double             reflection_coefficient_new;
                     // sail's reflection coefficient, when sail is new

  double             specular_reflection_coefficient_new;
                     // sail's coefficient for specular reflection, when sail is new

  double             emission_coefficient_front_new;
                     // sail's front side emission coefficient, when sail is new

  double             emission_coefficient_back_new;
                     // sail's back side emission coefficient, when sail is new

  double             reflection_coefficient_aged;
                     // sail's reflection coefficient, when sail is aged

  double             specular_reflection_coefficient_aged;
                     // sail's coefficient for specular reflection, when sail is aged

  double             emission_coefficient_front_aged;
                     // sail's front side emission coefficient, when sail is aged

  double             emission_coefficient_back_aged;
                     // sail's back side emission coefficient, when sail is aged

  double             reflection_coefficient_halflife;
                     // half life of sail's reflection coefficient

  double             specular_reflection_coefficient_halflife;
                     // half life of sail's coefficient for specular reflection

  double             emission_coefficient_front_halflife;
                     // half life of sail's front side emission coefficient

  double             emission_coefficient_back_halflife;
                     // half life of sail's back side emission coefficient

  double             non_lambertian_coefficient_front;
                     // sail's front side non-Lambertian coefficient

  double             non_lambertian_coefficient_back;
                     // sail's back side non-Lambertian coefficient

  // ONLY FOR NEP AND SEP

  double             sc_mass;
                     // spacecraft mass [kg]

  double             min_mass_initial_propellant;
                     // minimum intitial propellan mass [kg]

  double             max_mass_initial_propellant;
                     // maximum intitial propellan mass [kg]

  ThrottleType       throttletype;
                     // throttle type

  // ONLY FOR NEP

  double             thrust_max;
                     // maximum thrust [N]

  double             isp;
                     // specific impulse [s]

  // ONLY FOR SEP

  double             char_power;
                     // characteristic power P0 at Earth distance [kW]

  double             p_exponent;
                     // variation (exponent) p of solar array power with
                     // solar distance => P(r)/P0 ~ (1/r)^p

  double             thrusters;
                     // number of thrusters

private:
  void init(void); ///< initialize object
public:
  //-----------------------------------------------------------------

  /// standard constructor
  CInputSC(void) {init();};

  /// parameter list constructor
  CInputSC(ParamList& parameters, const std::string& trailing_chars = " \t");

  /// destructor
  ~CInputSC(void) {}

  /// resets input object
  void reset(void) {init();};

  /// prints current settings
  void printSettings(void) const;

  /// prints solar sail settings
  void printSolarSailSettings(void) const;

}; // end class CInputSC


#endif // CINPUTSPACECRAFT_H_INCLUDED
