/*!
\file Kepler.h
\version 1.7
\brief Miscellaneous functions for Keplerian elements.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Tanja D. Schmidt <tanja.schmidt@dlr.de>c (TS)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 19.08.2006 (AO)

\b Purpose:
   For the conversion of Keplerian elements certain mathematic
   formulaes, algorithms and operations are needed. The functions
   provided by this part of the SpaceLib - library are able to
   to these calculations and also do some error checking.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.
   This code is in parts based on the code from:
   Montenbruck, Gill: Satellite Orbits. Springer, 2000 and
   Montenbruck, Pfleger: Astronomie mit dem Personal Computer (3. Aufl.). Springer, 1999

\b Copyright:
   (c) 1999-2006 Oliver Montenbruck, Eberhard Gill, Thomas Pfleger, Andreas Ohndorf 
*/

#ifndef KEPLERSUPP_H
#define KEPLERSUPP_H

#include <string>
#include <algorithm>

#include <MathLib.h>

using namespace std;
namespace SpaceLib 
{

  /// Evaluates the Stumpff functions c1, c2 and c3 (for near-parabolic orbits)
  void Stumpff(double E2, double &c1, double &c2, double &c3);

  /*! \brief periapsis distance Rp = a(1-e)
   *
   * \param a semi-major axis of orbit
   * \param e eccentricity
   * \return periapsis distance in unit system of parameter a
   */
  template<class T, class U, class V>
  const T CalcPeriapsisDistance(const U a, const V e)
  {
    assert( ((a<0.0) && (e>1.0)) || ((a>=0.0) && (e<=1.0)) );
    return T(a) * ( T(1) - T(e) );
  }

  /*! \brief apoapsis distance Rp = a(1+e)
   *
   * \param a semi-major axis of orbit (>=0.0)
   * \param e eccentricity (<=1.0)
   * \return apoapsis distance in unit system of parameter a
   */
  template<class T, class U, class V>
  const T CalcApoapsisDistance(const U a, const V e)
  {
    assert( (a>=0.0) && (e<=1.0) );
    return T(a) * ( T(1) + T(e) );
  }

  /*! \brief periapsis velocity Vp = sqrt( GM*( (2/Rp) - (1/|a|) ) )
   *
   * \param a semi-major axis of orbit
   * \param e eccentricity
   * \param GM gravitational parameter of central body
   * \return Periapsis velocity in units depending on input parameter units
   */
  template<class T, class U, class V, class W> 
    const T CalcPeriapsisVelocity(const U a, const V e, const W GM)
  {
    return sqrt(T(GM)*(T(2.0)/CalcPeriapsisDistance<T>(a,e)-T(1.0)/T(abs(a))));
  }

  /*! \brief periapsis velocity Vp = sqrt( GM*( (2/Rp) - (1/|a|) ) )
   *
   * \param Rp periapsis distance
   * \param Ra apoapsis distance
   * \param GM gravitational parameter of central body
   * \return Periapsis velocity in units depending on input parameter units
   */
  template<class T, class U, class W>
    const T CalcPeriapsisVelocity2(const U Rp, const U Ra, const W GM)
  {
    return sqrt(T(2.0)*T(GM)*T(Ra)/(T(Rp)*(T(Ra)+T(Rp))));
  }


  /*! \brief mean motion = sqrt(GM/(a^3))
   * \param a semi-major axis
   * \param GM central body's gravitational parameter
   * \return mean motion (unit system is AUD or SI, depending on unit of GM)
   */
  template<class T,  class U>
  const T CalcMeanMotion(const T a, const U GM)
  {
    return sqrt(T(GM)/pow(a,int(3)));
  }

  /**
  * \brief Periapsis distance of orbit.
  *
  * \param orbitParameter Orbit-perifocus distance at true anomaly = PI/4
  * \param e The orbits numeric eccentricity
  */
  double PeriapsisDistance(double orbitParameter, double e);


  /// The orbit parameter from periapsis distance and eccentricity
  /**
    \param q Periapsis distance of Keplerian orbit.
    \param e The numeric eccentricity of Keplerian orbit.
  */
  double OrbitParameter(double q, double e);

  /// returns eccentric anomaly for given mean anomaly and eccentricity
  /**
     \param M mean anomaly
     \param e eccentricity
     \return the eccentric anomaly E
  */
  double EccentricAnomaly(double M, const double e);


  /// returns hyperbolic anomaly for given mean hyperbolic anomaly and eccentricity
  /**
     \param Mh mean anomaly
     \param e eccentricity
     \return the hyperbolic anomaly Eh
  */
  double HyperbolicAnomaly(double Mh, double e);


  /// Returns the hyperbolic anomaly for a given time \em t
  /**
   \param t0 The time of periapsis passage
   \param t  The desired point of time for which the anomaly is needed
   \param gm The gravitational parameter of the center body
   \param a  The semimajor axis of orbit (the absolute value will be used)
  */
  double HyperbolicMeanAnomaly(double t0, double t, double gm, double a);


  /// calculates eccentric anomaly and true anomaly (or hyperbolic values)
  /**
     Note: If the calculated state is parabolic, a fatal error will occur.

     \param a  semi major axis
     \param M  mean anomaly
     \param e  eccentricity
     \param EH calculated eccentric anomaly (reference)
     \param f  calculated true anomaly (reference)
     \return void
  */
  void CalcAnomaliesM(const double a, const double M, const double e, double &EH, double &f);


  /// calculates eccentric anomaly and mean anomaly (or hyperbolic values)
  /// from true anomaly
  /**
     Note: If the calculated state is parabolic, a fatal error will occur.

     \param a  semi major axis
     \param f  true anomaly
     \param e  eccentricity
     \param EH calculated eccentric anomaly (reference)
     \param M  calculated mean anomaly (reference)
     \return void
  */
  void CalcAnomaliesf(const double a, const double f, const double e, double &EH, double &M);

  /** Gives fatal error, if "theta" is out of range [-PI/2,+PI/2] */
  bool AssertElevationRange(const double theta);

  /** Gives fatal error, if "e" is < 0.0 */
  void AssertEccentricityRange(const double e);

  /** Gives fatal error, if "i" is out of range [0.0,pi] */
  void AssertInclinationRange(const double i);


  // ----------------------------------------------------------------------------
  /// identifier of Keplerian orbit element
  enum KEPLERIAN_ELEMENT
  {
    SEMIMAJORAXIS,       ///< semimajor axis
    ECCENTRICITY,        ///< eccentricity
    INCLINATION,         ///< inclination
    NODE,                ///< longitude of ascending node of orbit
    PERIAPSIS,           ///< argument of periapsis
    ANOMALY,             ///< mean anomaly angle
    N_KEPLERIAN_ELEMENTS
  };

  /// literal description of Keplerian elements
  const std::string KEPLERIAN_ELEMENT_NAME[N_KEPLERIAN_ELEMENTS] =
  {
    "semimajor axis a"
   ,"eccentricity e"
   ,"inclination i"
   ,"longitude of asc. node Omega"
   ,"argument of periapsis omega"
   ,"mean anomaly M"
  };

  // ----------------------------------------------------------------------------
  /*! \brief helper function for calculation of distance function values during proximity quotient calculation
  *
  * references.: 
  * [1] Anastassios E. Petropoulos, "Low-Thrust Orbit Transfers using Candidate
  * Lyapunov Functions with a Mechanism for Coasting", 2004, AIAA 2004-5089

  * \param element Keplerian element of which the difference w/r to the target element
  *should be calculated
  * \param osc_e current value of specified Keplerian element
  * \param target_e desired target value of specified Keplerian element
  * \param max_error maximum allowed absolute difference of osc_e to target_e
  * \return distance function value according to ref. [1]
  * \sa CalcScalingFactor
  */
  inline double CalcDistanceFunction(const KEPLERIAN_ELEMENT element, 
    const double osc_e, 
    const double target_e,
    const double max_error)
  {
    double diff = 0.0;
    switch (element)
    {
    case SEMIMAJORAXIS:
    case ECCENTRICITY:
    case INCLINATION:
      diff = (osc_e-target_e);      
      break;
    case NODE:
    case PERIAPSIS:
      diff = acos(cos(osc_e - target_e));
      break;
    default: return 0; break;
    }
    if (fabs(diff)<=max_error) return 0.0;
    else 
    {
      if (diff<0.0) return (diff+max_error);
      else return (diff-max_error);
    }
  }

  // ----------------------------------------------------------------------------
  /*! \brief helper function for calculation of temporary scaling factors during proximity quotient calculation
  *
  * references.: 
  * [1] Anastassios E. Petropoulos, "Low-Thrust Orbit Transfers using Candidate
  * Lyapunov Functions with a Mechanism for Coasting", 2004, AIAA 2004-5089

  * \param element Keplerian element of which the scaling factor should be calculated
  * \param osc_e current value of specified Keplerian element
  * \param target_e desired target value of specified Keplerian element
  * \return scaling factor according to ref. [1]
  * \sa CalcDistanceFunction
  */
  inline double CalcScalingFactor(const KEPLERIAN_ELEMENT element,
    const double osc_e, 
    const double target_e)
  {
    switch (element)
    {
    case SEMIMAJORAXIS:
      return(sqrt(1.0+pow( (osc_e-target_e)/(3.0*target_e), int(4) ) ) );
      break;
    case ECCENTRICITY:
    case INCLINATION:
    case NODE:
    case PERIAPSIS:
      return 1.0;
      break;
    default: return 1.0; break;
    }
  }

  // ----------------------------------------------------------------------------
  /*! \brief helper function for calculation of maximum rate of change of 
  semimajor axis depending on current state
  *
  * references.: 
  * [1] Anastassios E. Petropoulos, "Low-Thrust Orbit Transfers using Candidate
  * Lyapunov Functions with a Mechanism for Coasting", 2004, AIAA 2004-5089
  *\n
  * Calculates the maximum achievable rate of change of semimajor axis
  * according to (12) from [1]. The returned value must be 
  * multiplied with the s/c current absolute acceleration in order to obtain 
  * the final value.
  *
  * \param a semimajor axis [AU]
  * \param e eccentricity
  * \param mu gravitational parameter of central body
  * \return normed maximum rate of change of semimajor axis 
  *(multiply with s/c acceleration)
  * \sa CalcMaxElementRateOfChange, CalcMaxNodeRateOfChange
  * \sa CalcMaxPeriapsisRateOfChange, CalcMaxEccentricityRateOfChange*/
  inline double CalcMaxSemimajorAxisRateOfChange(const double e, const double a,
    const double mu)
  {
    return 2.0*sqrt(a*a*a*(1.0+e)/(mu*(1.0-e)));
  }


  // ----------------------------------------------------------------------------
  /*! \brief helper function for calculation of maximum eccentricity rate of 
  change depending on current state
  *
  * references.: 
  * [1] Anastassios E. Petropoulos, "Low-Thrust Orbit Transfers using Candidate
  * Lyapunov Functions with a Mechanism for Coasting", 2004, AIAA 2004-5089
  *\n
  * Calculates the maximum achievable eccentricity rate of change
  * according to (13) from [1]. The returned value must be multiplied with the 
  * s/c current absolute acceleration in order to obtain the final value.
  *
  * \param p semilatus rectum of current orbit
  * \param h specific orbit momentum
  * \return normed maximum eccentricity rate of change (multiply with s/c 
  acceleration)
  * \sa CalcMaxElementRateOfChange, CalcMaxNodeRateOfChange
  * \sa CalcMaxPeriapsisRateOfChange, CalcMaxInclinationRateOfChange
  * \sa CalcMaxSemimajorAxisRateOfChange*/
  inline double CalcMaxEccentricityRateOfChange(const double p,
    const double h)
  {
    return 2.0*p/h;
  }

  // ----------------------------------------------------------------------------
  /*! \brief helper function for calculation of maximum inclination rate of 
  change depending on current state
  *
  * references.: 
  * [1] Anastassios E. Petropoulos, "Low-Thrust Orbit Transfers using Candidate
  * Lyapunov Functions with a Mechanism for Coasting", 2004, AIAA 2004-5089
  *\n
  * Calculates the maximum achievable inclinationrate of change
  * according to (14) from [1]. The returned value must be multiplied with the 
  * s/c current absolute acceleration in order to obtain the final value.
  *
  * \param e orbit eccentricity
  * \param h specific orbit momentum
  * \param om orbit argument angle of perapsis
  * \param p orbit semilatus rectum
  * \return normed maximum inclination rate of change (multiply with s/c 
  acceleration)
  * \sa CalcMaxElementRateOfChange, CalcMaxNodeRateOfChange
  * \sa CalcMaxPeriapsisRateOfChange, CalcMaxEccentricityRateOfChange
  * \sa CalcMaxSemimajorAxisRateOfChange*/
  inline double CalcMaxInclinationRateOfChange(const double e,
    const double h, const double om, const double p)
  {
    const double eps = 1.0e-20;
    double t = h*( sqrt( 1.0-pow(e*sin(om),int(2)) ) - e*fabs(cos(om)) );
    if (fabs(t)<eps) t = signum(t)*eps;
    return p/t;
  }

  // ----------------------------------------------------------------------------
  /*! \brief helper function for calculation of maximum rate of change of 
  node angle depending on current state
  *
  * references.: 
  * [1] Anastassios E. Petropoulos, "Low-Thrust Orbit Transfers using Candidate
  * Lyapunov Functions with a Mechanism for Coasting", 2004, AIAA 2004-5089
  *\n
  * Calculates the maximum achievable rate of change of argument angle of 
  * ascending node according to (18) from [1]. The returned value must be multiplied
  * with the s/c current absolute acceleration in order to obtain the final value.
  *
  * \param e orbit eccentricity
  * \param h specific orbit momentum
  * \param i orbit inclination
  * \param om argument of periapsis
  * \param p semilatus rectum of orbit
  * \return normed maximum ascending node rate of change (multiply with s/c acceleration)
  * \sa CalcMaxElementRateOfChange, CalcMaxPeriapsisRateOfChange
  * \sa CalcMaxSemimajorAxisRateOfChange, CalcMaxEccentricityRateOfChange
  * \sa CalcMaxInclinationRateOfChange*/
  inline double CalcMaxNodeRateOfChange(const double e, const double h,
    const double i, const double om, const double p)
  {
    const double eps = 1.0e-20;
    double t = h*sin(i)*
      (sqrt(1.0-pow(e*cos(om),int(2)))-e*fabs(sin(om)));
    if ( fabs(t) < eps ) t = signum(t)*eps;
    return p/t;
  }

  // ----------------------------------------------------------------------------
  /*! \brief helper function for calculation of maximum rate of change of 
  periapsis angle depending on current state
  *
  * references.: 
  * [1] Anastassios E. Petropoulos, "Low-Thrust Orbit Transfers using Candidate
  * Lyapunov Functions with a Mechanism for Coasting", 2004, AIAA 2004-5089
  *\n
  * Calculates the maximum achievable rate of change of argument angle of periapsis
  * according to (15), (16) and (17) from [1]. The returned value must be multiplied
  * with the s/c current absolute acceleration in order to obtain the final value.
  *
  * \param e orbit eccentricity
  * \param h specific orbit momentum
  * \param p semilatus rectum of orbit
  * \return normed maximum periapsis rate of change (multiply with s/c acceleration)
  * \sa CalcMaxElementRateOfChange, CalcMaxNodeRateOfChange
  * \sa CalcMaxSemimajorAxisRateOfChange, CalcMaxEccentricityRateOfChange
  * \sa CalcMaxInclinationRateOfChange*/
  inline double CalcMaxPeriapsisRateOfChange(const double e, const double h,
    const double p)
  {
    const double t1 = (1-e*e)/(e*e*e);
    const double t2 = sqrt(t1*t1/4.0+1.0/27.0);
    const double cost = pow(t1/2.0+t2, 1.0/3.0)-pow(-t1/2.0+t2, 1.0/3.0)-1/e;
    const double sint2 = max(0.0, 1.0-cost*cost);
    const double r = p/(1.0+e*cost);
    return 1.0/(e*h)*sqrt(pow(p*cost,int(2))+pow((p+r),int(2))*sint2);
  }

} // namespace spacelib

#endif // of KEPLERSUPP_H
// end of file

