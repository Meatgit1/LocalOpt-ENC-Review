/*!
\file SpaceConst.h
\version 1.7
\brief Constants for astrodynamical calculations.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 17.10.2006 (AO)

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

#ifndef SPACECONST_H
#define SPACECONST_H

#include <string>

namespace SpaceLib
{
  /** SpaceLib unit system.
   */
  enum UnitType {
    SI = 0,    ///< distance in m, time in s
    AUD = 1,   ///< distance in AU, time in days [d]
    STORED = 2 ///< default, indicates that the stored unit type should be used
  };

  /** SpaceLib distance units.
   */
  const std::string DISTANCE_UNIT[] =
  {
    "m", ///< meter
    "AU" ///< astronomical unit
  };

  /** SpaceLib velocity units.
   */
  const std::string VELOCITY_UNIT[] =
  {
    "m/s", ///< meter per second
    "AU/d" ///< astronomical unit per day
  };

  // ----------------------------------------------------------------------------
  // solar system constants (part 1)
  // ----------------------------------------------------------------------------

  /** mean distance sun-earth (1 astronomical unit, 1 AU)\n
   *  Source: JPL Solar System Dynamics Group Website, www.jpl.nasa.gov\n
   *  Date  : 04.10.2003 (AO)
   */
  const double AU = 149597870691.0; ///< mean distance earth - sun [m]

  // ----------------------------------------------------------------------------
  // conversion factors (part 1)
  // ----------------------------------------------------------------------------

  // distance
  const double AU2KM = AU / 1000.0; ///< conversion factor: [AU]=>[km]
  const double KM2AU = 1.0 / AU2KM; ///< conversion factor: [km]=>[AU]
  const double AU2M = AU; ///< conversion factor: [AU]=>[m]
  const double M2AU = 1.0 / AU2M; ///< conversion factor: [m]=>[AU]

  // velocity
  const double AUD2MS = AU2M / 86400.0; ///< conversion factor: [AU/d]=>[m/s]
  const double MS2AUD = 1.0 / AUD2MS; ///< conversion factor: [m/s]=>[AU/d]
  const double AUD2KMS = AU / 1000.0 / 86400.0; ///< conversion factor: [AU/d]=>[km/s]
  const double KMS2AUD = 1.0 / AUD2KMS; ///< conversion factor: [km/s]=>[AU/d]

  // acceleration
  const double AUD22MS2 = AU / 86400.0 / 86400.0; ///< conversion factor: [AU/d²]=> [m/s²]
  const double MS22AUD2 = 1.0 / AUD22MS2; ///< conversion factor: [m/s²] => [AU/d²]

  // ----------------------------------------------------------------------------
  // physical constants
  // ----------------------------------------------------------------------------

  /** speed of light in vacuum\n
   *  Source: IAU 76
   */
  const double C_SI = 299792458.0;  ///< value of c in SI units [m/s]
  const double C_AUD = C_SI * 86400 / AU; ///< value of c in AUD units [AU/d]


  /** universal gravitational constant
   *  Source: JPL Solar System Dynamics Group Website, ssd.jpl.nasa.gov
   */
  const double G_SI = 6.67259e-11;  ///< gravitational constant G in SI units [m³/(kg*s²)]
  const double G_AUD = G_SI * 86400.0 * 86400.0 / (AU * AU * AU); ///< gravitational constant G in AUD units [AU³/(kg*d²)]

  /** Stefan-Boltzmann constant [m³/(kg*s²)]\n
   *  Source: HP48GX scientific calculator\n
   *  Date  : 22.12.2004 (BD)
   */
  const double SIGMA_AUD = 36573232.6048; // Stefan-Boltzmann constant in AUD [kg/(d^3*K^4)]

  // ----------------------------------------------------------------------------
  // conversion factors (part 2)
  // ----------------------------------------------------------------------------

  const double GM_AUD2GM_SI = (AU * AU * AU) / (86400.0 * 86400.0); ///< conversion factor: [AU³/d²]=> [m³/s²]
  const double GM_SI2GM_AUD = 1.0 / GM_AUD2GM_SI; ///< conversion factor: [m³/s²]=> [AU³/d²]

  // ----------------------------------------------------------------------------
  // solar system constants (part 2)
  // ----------------------------------------------------------------------------

  /** heliocentric gravitational parameter G*M
   *  Source: JPL Solar System Dynamics Group Website, ssd.jpl.nasa.gov
   */
  // our value:
  const double GM_SUN_SI = 1.32712440018e20;  ///< G*M (sun) in SI units [m³/s²]
//   const double GM_SUN_SI = 1.32713440018e20;  ///< G*M (sun) in SI units [m³/s²]
  // ACT's value:
  // const double GM_SUN_SI = 1.32712428e20;  ///< G*M (sun) in SI units [m³/s²]
  const double GM_SUN_AUD = GM_SUN_SI * GM_SI2GM_AUD; ///< G*M (sun) in AUD units [AU³/d²]

  /** geocentric gravitational parameter G*M
   *  Source: JGM-3
   */
  const double GM_EARTH_SI = 3.986004415e14;  ///< G*M (earth) in SI units [m³/s²]
  const double GM_EARTH_AUD = GM_EARTH_SI * GM_SI2GM_AUD; ///< G*M (earth) in AUD units [AU³/d²]

  /** selenocentric gravitational parameter
   *  Source: DE405
   */
  const double GM_MOON_SI = GM_EARTH_SI / 81.30056; ///< G*M (moon) in SI units [m³/s²]

  /** selenocentric gravitational parameter
   *  Source: DE405
   */
  const double GM_MOON_AUD = GM_MOON_SI * GM_SI2GM_AUD; ///< G*M (earth) in AUD units [AU³/d²]

  /** solar radiation pressure at 1 AU
   *  Source: IERS 96
   */
  const double P_SUN = 4.560E-6;  ///< solar radiation pressure at 1 AU in SI units [N/m²] (radiation power ~ 1367 W/m²)

  /** solar radiation flux at 1 AU
   */
  const double S0_SI = 1367.0; // solar constant in AUD units [W/m²=kg/s³] (radiation power ~ 1367 W/m²)
  const double S0_AUD = 8.82322440192E+17; // solar constant in AUD units [kg/d³] (radiation power ~ 1367 W/m²)

  /** earth (standard) gravitational acceleration [m/s²]\n
   *  Source: HP48GX scientific calculator
   */
  const double G0 = 9.80665;    ///< earth gravitational acceleration [m/s²]

  /** obiquity of earth orbit to ecliptic in [arcsec]\n
   *  Source : http://ssd.jpl.nasa.gov/astro_constants.html
   *  Date   : 04.10.2003 (AO)
   */
  const double EPS0_J2000 = 84381.412;  ///< J2000-obliquity of earth orbit [arcsec]

  typedef int BodyID; ///< tbd

  //  enum PlanetIndex {
  const BodyID SUN = 1;                       ///< Sun
  const std::string SUN_NAME = "SUN";         ///< Sun

  const BodyID MERCURY = 2;                   ///< Mercury
  const std::string MERCURY_NAME = "MERCURY"; ///< Mercury

  const BodyID VENUS = 3;                     ///< Venus
  const std::string VENUS_NAME = "VENUS";     ///< Venus

  const BodyID EARTH = 4;                     ///< Earth
  const std::string EARTH_NAME = "EARTH";     ///< Earth

  const BodyID MARS = 5;                      ///< Mars
  const std::string MARS_NAME = "MARS";       ///< Mars

  const BodyID JUPITER = 6;                   ///< Jupiter
  const std::string JUPITER_NAME = "JUPITER"; ///< Jupiter

  const BodyID SATURN = 7;                    ///< Saturn
  const std::string SATURN_NAME = "SATURN";   ///< Saturn

  const BodyID URANUS = 8;                    ///< Uranus
  const std::string URANUS_NAME = "URANUS";   ///< Uranus

  const BodyID NEPTUNE = 9;                   ///< Neptune
  const std::string NEPTUNE_NAME = "NEPUNE"; ///< Neptune

  const BodyID PLUTO = 10;                    ///< Pluto
  const std::string PLUTO_NAME = "PLUTO";     ///< Pluto

  const BodyID LUNA = 11;                     ///< our Moon
  const std::string LUNA_NAME = "moon";       ///< our Moon

  const BodyID SS_BARY = 0;                   
            ///< solar system's barycenter (the origin of the coordinate system)
  const std::string SS_BARY_NAME = "solar barycenter";
            ///< solar system's barycenter (the origin of the coordinate system)

  /*! \brief usable coordinate parameter types
   * It shows, which kind of coordinate elements are used. The main types are
   * cartesian coordinates and different forms of keplerian elements. Theangles
   * and angular velocities are expressed in rad and rad/s whereas the unit of
   * distance and velocity values is usually indicated by an UnitType variable.
   * \sa UnitType
   */
  enum CoordType 
  {
    /*! \brief Cartesian coordinates and velocities.
     * position vector (x,y,z)\n
     *   velocity vector (dx,dy,dz)\n\n */
    CARTESIAN,

    /*! \brief speherical coordinates and velocities.
     *
     * position vector (phi,theta,r)\n
     * velocity vector (dphi,dtheta,dr)\n
     * phi   (0,2PI]\n
     * theta (-PI/2,+PI/2)\n\n */
    POLAR,

    /*! \brief Kepler elements (mean anomaly M).
     *
     * semi-major axis a[described by a UnitType value]\n
     *  eccentricity e[]\n
     *  inclination i[rad]\n
     *  longitude of ascending node (capital Omega) [rad])\n
     *  argument of pericenter (small omega) [rad])\n
     *  mean anomaly M[rad]\n\n */
    KEPLERM,

    /*! \brief Kepler elements (true anomaly f).
     *
     * semi-major axis a[described by a UnitType value]\n
     * eccentricity e[]\n
     * inclination i[rad]\n
     * longitude of ascending node Omega[rad])\n
     * argument of pericenter[rad])\n
     * true anomaly f[rad]\n\n */
    KEPLERF,

    /*! \brief Kepler elements (mean longitude L).
     *
     * semi-major axis a[described by a UnitType value]\n
     * eccentricity e[]\n
     * inclination i[rad]\n
     * longitude of ascending node Omega[rad])\n
     * argument of pericenter[rad])\n
     * mean longitude L[rad]\n */
    ELEMENTS
  }; // end of CoordType

  /*! \brief Type of coordinate system for that given coordinates are valid for.
   *
   * It indicates, what reference coordinate system is valid for the given
   * or returned coordinates/elements.*/
  enum CoordSysType 
  {

    /*! \brief heliocentric, ecliptical, nonrotating (J2000 reference frame; 
     *  \brief default)
     *
     * The base frame is the earth orbit plane of the julian date 2451545.0
     *   The x-axis is determined by the vernal equinox and the z-axis by 
     * the direction of the Earth's angular momentum vector. This is the 
     * default coodinate system and used for internal representation.*/
    HEN,

    /*! \brief heliocentric, orbital(earth equatorial), nonrotating
     *
     * This reference frame is a special one, because it is dependend to the
     * Earth's equator. The x-axis points to the vernal equinox as it does 
     * in the heliocentric ecliptical reference frame but is rotated by 
     * the obliquity of the Earth. It function-dependent whether the 
     * obliquity will be of date of J2000.*/
    HON,

    /*! \brief bodycentric, ecliptical, nonrotating
     *
     * The basis plane of this reference frame still the ecliptic of J2000
     * but the origin of the coordinate system is not the Sun or the solar
     * system's barycenter. The new coordinate center is the position of date
     * of the mentioned body.*/
    BEN,

    /*! bodycentric, ecliptical, rotating
     *
     * The reference base plane is the J2000 eclliptic. The x-axis is the 
     * projection of the position vector onto the base plane.\n
     * Transformation: Rot_z(atan(position(t).y,position(t).x))*/
    BER,

    /*! bodycentric, orbital, nonrotating
     *
     * The reference plane is the bodies orbit plane and the x-axis direction
     * is determined by the periapsis direction in the body's orbit. The
     * z-axis is the direction of the angular momentum vector and thus
     * perpendicular to the base plane. The coordinate frame should be not
     * but through the possibility of using osculating Keplerian elements, it
     * can not be assured that the frame is time invariant. (the rotating
     * angles are often build from the Keplerian elements of date!!!)\n
     * Transformation:\n
     * Rot_z(argument_of_periapsis(t)) *\n
     * Rot_x(inclination(t)) *\n
     * Rot_z(longitude_of_ascending_node(t)) */
    BON,

    /*! \brief bodycentric, orbital, rotating (with true anomaly)
     *
     * The reference base plane is the orbit plane and the x-axis direction
     * is the vector from the center body to the actual body. So, different
     * times also mean different x-axis. The coordinate system is rotating
     * with the true anomaly in the orbit plane.\n
     * Transformation:\n
     *  Rot_z(argument_of_periapsis(t)+true_anomaly(t)) *\n
     * Rot_x(inclination(t)) *\n
     *  Rot_z(longitude_of_ascending_node(t))*/
    BORT
  }; // end of CoordSysType


  /*! \brief Possible state types
   *
   * In context of celestial mechanics a state my be virtually everything
   * that can be described by means of distance/velocity vectors and/or 
   * Keplerian elements.
   */
  enum StateType
  {
    /*! \brief state is expressed by a celestial body
     *
     * The state is incorporated by an existing body. "Existing" in the context
     * of "SpaceLib" means the body exists within a CSolarSystem object. The
     * body's position and velocity can be obtained from this object and so
     * the state basically is the reference body's state of date.
     */
    BODY_STATE,

    /*! \brief state expressed by the orbit of a celestial orbit
     *
     * The reference body's orbit is used to describe the state. So the 
     * state parameters differ from those of the reference body only in 
     * the in-orbit angular position, generally described with help of the
     * anomaly (mean, true).*/
    BODY_ORBIT_STATE,

    /*! \brief state expressing a user-defined orbit about a given central body
     *
     * In this case a free orbit is described by position/velocity vector or
     * a set of Keplerian elements w.r.t. a given central body.
     * \sa ORBITAL_ELEMENTS_STATE
     */
    FREE_ORBIT_STATE,

    /*! \brief state expressing a user-defined orbit defined by Keplerian 
    elements
     * This state type differs from FREE_ORBIT_STATE in that it may only be used
     * in conjunction with Keplerian elements, which must'nt be complete. The
     * use of position/velocity vector combination is not allowed with it.
     *
     * \sa FREE_ORBIT_STATE
     */
    ORBITAL_ELEMENTS_STATE,

    /*! \brief state that is described by position/velocity vector
     * 
     * The state is expressed by a possibly not allowed combination
     * of a position and a velocity vector. Not allowed means in this
     * case that this state is not from an orbit conforming to Newton's
     * Law of gravitation.
     */
    VECTOR_STATE,

    /*! \brief reserved*/
    RESERVED_01,

    /*! \brief reserved*/
    RESERVED_02,

    /*! \brief reserved*/
    RESERVED_03
  }; // enum StateType

  /// target state type description
  const std::string StateTypeDesc[] =
  {
     "body"
    ,"body's orbit"
    ,"free orbit"
    ,"orbital elements"
    ,"state vector"
    ,"reserved 01"
    ,"reserved 02"
    ,"reserved 03"
  };

} // end of namespace SpaceLib

#endif // SPACECONST_H_INCLUDED
