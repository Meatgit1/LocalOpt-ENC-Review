/*!
\file CSolarSystem.h
\version 1.7
\brief A class that provides the functionality of the whole solar system (CSolarSystem).
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Tanja D. Schmidt <tanja.schmidt@dlr.de>c (TS)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 07.11.2006 (AO)

\b Purpose:
   A CSolarSystem object is intended to keep the main bodies of the solar
   system in its internal data structures and provide easy access to all the
   bodies data. Additionaly there are some useful functions implemented, e.g.
   the HILL-Radius of a solar system body.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.
   This code is based on the code from:
   Montenbruck, Gill: Satellite Orbits. Springer, 2000 and
   Montenbruck, Pfleger: Astronomie mit dem Personal Computer (3. Aufl.). Springer, 1999

\b Copyright:
   (c) 2001-2006 Bernd Dachwald, Tanja D. Schmidt, Andreas Ohndorf
*/

#ifndef CSOLARSYSTEM_H
#define CSOLARSYSTEM_H

// STL
#include <stdexcept>
#include <string>

namespace SpaceLib
{

  /**
    \brief Number of tries for inserting a body into the CSolarSystem - object

    This is the maximum number of attempts to insert a CBody object into
    the internal data structures of CSolarSystem.

    \sa insertbody
  */
  const int N_O_TRS = 1000;

  // internal data structure
  struct CSolarSystemData;

  /// The class that contains the bodies of the solar system.
  class CSolarSystem
  {

  private:

    CSolarSystemData * solarSystemData;

  public:

    /**
    * \brief The constructor.
    *
    * The constructor initializes the solar system object with the nine planets.
    * They can be accessed with an ID that can be obtained from 'GetPlanetIndex'
    * with the planets names in 'SpaceConst.h'. The file whose name is given by
    * the parameter DE405FileName will be used to get the position and velocity
    * informations of these planets. If no filename is provided, these
    * information will be derived from Keplerian elements.
    *
    * \param DE405FileName The absolute path and filename of the data file that
    *                      contains the DE405 data from JPL
    */
    CSolarSystem(const std::string & DE405FileName = "")
      throw (std::bad_alloc&);

    /// Destructor.
    ~CSolarSystem(void) 
    {
      SafeDelete(solarSystemData);
    }

    /** Tbd. */
    DE405* getEphemObj(void);

   /**
   * \brief Sets the JPL's comet data file.
    The mentioned file, usually with the name 'ELEMENTS.COMET', is obtainable
    at the JPL's internet site:
    http://ssd.jpl.nasa.gov.
    The file should not be modified by the user because the function of this
    method strongly depends on the the untouched data file. It is a plain text
    file so that you can also use the values of this file for a specific comet,
    initialize a CState object with dynamic memory allocation and pass the
    resulting pointer to the method insertBody. The easier way is to use
    insertComet with the identifying text string you can get from the data file.
    \param fileName The absolute path to a file containing the Keplerian orbit
                    parameters of many comets.
    \sa insertComet
   */
    void setCometDataFile(const std::string & fileName) throw (invalid_argument);


   /**
   * \brief Sets the JPL's asteroid data file.
   *
    The mentioned file, usually with the name 'ELEMENTS.NUMBR' or
    'ELEMENTS.UNNUM', can be downloaded from JPL's internet site:
    http://ssd.jpl.nasa.gov.
    The file should not be modified by the user because the function of this
    method strongly depends on the the untouched data file. It is a plain text
    file so that you can also use the values of this file for a specific comet,
    initialize a CState object with dynamic memory allocation and pass the
    resulting pointer to the method insertBody. The easier way is to use
    insertAsteroid with the identifying text string you can get from the data
    file.
    \param fileName The absolute path to a file containing the Keplerian orbit
                    parameters of many comets.
    \sa insertAsteroid
   */
    void setAsteroidDataFile(const std::string & fileName) throw (invalid_argument);


   /**
   * \brief Inserts a certain body into the solar system class.
     This member method allows to insert a user provided space object into
     the data structures of CSolarSystem. All data to be used for the body,
     e.g. orbit parameters, vectors, must be in the CBody - object the argument
     pointer points to.

    \param body  Pointer to a CBody object to insert in solar system.
    \param force Usually the method fails if the given CBody - object contains
                 an ID that is already in the internal data structures of the
                 CSolarSystem - object. That means another body owns the the
                 desired ID. The result of an insertBody function call will
                 normally by a return value of 0 and therefore an error. You
                 may override this by setting the parameter \em force to true.
                 The CBody - object's ID will be ignored and another one will
                 be tried. This is done "attempts" times before a value
                 of 0 will be returned as indication of an error.
    \param attempts The number of tries, the function should trie to insert
                    the given object, before a possible error is returned. This
                    parameter is only useful in connection with force.
    \return A unique ID that identifies the given body (0 if insertation failed).
   */
    BodyID insertBody(CBody * body, bool force=false, int attempts = N_O_TRS);


    /*! \brief remove a body from solar system
     *
     *  \param body ID of the body to be removed from the CSolarSystem object.*/
    void removeBody(const BodyID body);

    /// number of bodies in current solar system
    const int size(void) const;


   /** Inserts the desired comet into the solar system and therefor makes it accessible.
       \param stringID The identifiying string in data file. */
    BodyID insertComet(std::string stringID);


   /** Inserts the desired asteroid into the solar system and thus makes it accessible.
       \param stringID The identifiying string in data file. */
    BodyID insertAsteroid(std::string stringID);

   /**
   * \brief Returns a pointer to the desired CBody object that is
   *        indicated by a given ID.
   *
   *  If the body indicated by the parameter bodyID is existent within the CSolarSystem
   *  class object, a pointer to the corresponding CBody class object will be returned.
   *  Otherwise the method returns 0.
   *
   * \param bodyID The ID of the planet you want to receive a pointer of.
   * \return A pointer to the CBody class object of the body, zero if the
   *         body is not existent within the CSolarSystem class object.
   */
    CBody* getBody(const BodyID bodyID) const;


    /// Returns the number of bodies stored in the CSolarSystem-object.
    int NumberOfPlanets(void) const;

    /// returns the unique ID of the body in the solar system
    int Index(const BodyID planet) const {
      return getBody(planet)->getId();
    };

    /// Returns the name of a object stored in CSolarSystem-object.
    const char *Name(const BodyID planet) const;

    /// name of celestial body
    const std::string& getBodyName(const BodyID body_id) const;

    /// Returns the gravitational parameter of the desired body in SI - units.
    double GM_SI(const BodyID planet) const {
      return getBody(planet)->getGM(SI);//GM_SI();
    };

    /// Returns the gravitational parameter of the desired body in AUD - units.
    double GM_AUd(const BodyID planet) const {
      return getBody(planet)->getGM(AUD);//GM_AUd();
    };

    /// Returns a reference to a planets data subset incorporated
    /// by a CState-class.
    /**
      \param planet The ID of the desired body in solar system.
      \param mjd The date for which the CState data should be valid for.
    */
    CState & State(const BodyID planet, const double mjd) {
      return getBody(planet)->State(mjd);
    };


    /// Returns the x-coordinate of the desired planet.
   /**
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
   */
    double x(const BodyID planet, const double mjd) {
      return getBody(planet)->x(mjd);
    };


    /// Returns the y-coordinate of the desired planet.
   /**
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
   */
    double y(const BodyID planet, const double mjd) {
      return getBody(planet)->y(mjd);
    };


    /// Returns the z-coordinate of the desired planet.
   /**
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
   */
    double z(const BodyID planet, const double mjd) {
      return getBody(planet)->z(mjd);
    };


    /// Returns the planets velocity x-component at the given date.
   /**
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
   */
    double dx(const BodyID planet, const double mjd) {
      return getBody(planet)->dx(mjd);
    };


    /// Returns the planets velocity y-component at the given date.
   /**
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
   */
    double dy(const BodyID planet, const double mjd) {
      return getBody(planet)->dy(mjd);
    };


    /// Returns the planets velocity z-component at the given date.
   /**
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
   */
    double dz(const BodyID planet, const double mjd) {
      return getBody(planet)->dz(mjd);
    };


    /// Returns the distance of the desired body from its reference system origin at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
   */
    double r(const BodyID planet, const double mjd) {
      return getBody(planet)->r(mjd);
    };


    /// Returns the azimuth angle of the desired body in its reference system at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    /return Azimuth angle phi (0, 2PI]
   */
    double phi(const BodyID planet, const double mjd) {
      return getBody(planet)->phi(mjd);
    };


    /// Returns the elevation angle of the desired body in its reference system at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return Elevation angle theta (-PI/2,PI/2).
   */
    double theta(const BodyID planet, const double mjd) {
      return getBody(planet)->theta(mjd);
    };


    /// Returns the radial velocity of the desired body in its reference system at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return Radial velocity in the units of the CState-object.
   */
    double dr(const BodyID planet, const double mjd) {
      return getBody(planet)->dr(mjd);
    };


    /// Returns the angular velocity phi dot of the desired body in its reference system at the given date.

   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return Angular velocity phi dot in the units of the CState-object (rad/s, rad/d).
   */
    double dphi(const BodyID planet, const double mjd) {
      return getBody(planet)->dphi(mjd);
    };


    /// Returns the angular velocity theta dot of the desired body in its reference system at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return Angular velocity theta dot in the units of the CState-object (rad/s, rad/d).
   */
    double dtheta(const BodyID planet, const double mjd) {
      return getBody(planet)->dtheta(mjd);
    };


    /// Returns the semimajor axis a of the desired body's orbit in its reference system at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return semimajor axis a in the units of the CState-object (m, AU).
   */
    double a(const BodyID planet, const double mjd) {
      return getBody(planet)->a(mjd);
    };


    /// Returns the numerical eccentricity e of the desired body's orbit in its reference system at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return Numerical eccentricity e (dimensionless value).
   */
    double e(const BodyID planet, const double mjd) {
      return getBody(planet)->e(mjd);
    };


    /// Returns the orbit inclination i of the desired body in its reference system at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return Orbit's inclination i in rad.
   */
    double i(const BodyID planet, const double mjd) {
      return getBody(planet)->i(mjd);
    };


    /// Returns the longitude of the body's orbits ascending node in its reference system at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return Longitude of ascending node in rad.
   */
    double Omega(const BodyID planet, const double mjd) {
      return getBody(planet)->Omega(mjd);
    };


    /// Returns the argument of the body's orbits pericenter in its reference system at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return Argument of pericenter in rad.
   */
    double omega(const BodyID planet, const double mjd) {
      return getBody(planet)->omega(mjd);
    };

    /// Returns Omega + omega
    double curlypi(const BodyID planet, const double mjd) {
      return getBody(planet)->curlypi(mjd);
    };


    /// Returns the mean anomaly of the body's orbit position in its reference system at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return Mean anomaly in rad.
   */
    double M(const BodyID planet, const double mjd) {
      return getBody(planet)->M(mjd);
    };


    /// Returns the eccentric anomaly of the body's orbit position in its reference system at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return Eccentric anomaly in rad.
   */
    double E(const BodyID planet, const double mjd) {
      return getBody(planet)->E(mjd);
    };


    /// Returns the true anomaly of the body's orbit position in its reference system at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return True anomaly in rad.
   */
    double f(const BodyID planet, const double mjd) {
      return getBody(planet)->f(mjd);
    };


    /// Returns mean longitude = Omega+omega+M
   /**
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return Mean longitude in rad.
    \sa Omega, omega, M
   */
    double L(const BodyID planet, const double mjd) {
      return getBody(planet)->L(mjd);
    };


   /**
   * \brief The specific angular momentum.
   *
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return Specific angular momentum in in the units of the CState-object.
   */
    double h(const BodyID planet, const double mjd) {
      return getBody(planet)->h(mjd);
    };

    /** Tbd. */
    double k(const BodyID planet, const double mjd) {
      return getBody(planet)->k(mjd);
    };

    /** Tbd. */
    double p(const BodyID planet, const double mjd){
      return getBody(planet)->p(mjd);
    };

    /** Tbd. */
    double q(const BodyID planet, const double mjd){
      return getBody(planet)->q(mjd);
    };

    /** Tbd. */
    double l(const BodyID planet, const double mjd) {
      return getBody(planet)->l(mjd);
    };

    /** Tbd. */
    double Velocity(const BodyID planet, const double mjd) {
      return getBody(planet)->Velocity(mjd);
    };

    /** Tbd. */
    double Momentum(const BodyID planet, const double mjd) {
      return getBody(planet)->Momentum(mjd);
    };

    /** Tbd. */
    double KineticEnergy(const BodyID planet, const double mjd) {
      return getBody(planet)->KineticEnergy(mjd);
    };

    /** Tbd. */
    double PotentialEnergy(const BodyID planet, const double mjd) {
      return getBody(planet)->PotentialEnergy(mjd);
    };


    /**
    * \brief The body's mass specific total energy.
    */
    double Energy(const BodyID planet, const double mjd) {
      return getBody(planet)->Energy(mjd);
    };

    /** Tbd. */
    double OrbitPeriod(const BodyID planet, const double mjd) {
      return getBody(planet)->OrbitPeriod(mjd);
    };

    /** Tbd. */
    double MeanMotion(const BodyID planet, const double mjd) {
      return getBody(planet)->MeanMotion(mjd);
    };


    /// Returns the semiminor axis a of the desired body's orbit in its reference system at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with a origin depending
    on the reference body. Usually this reference body is the sun and therefor the reference
    system is often a heliocentric, eclliptical one. But it mustn't, the earth's natural
    satellite, the moon, for example uses the earth as reference system.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return semiminor axis b in the units of the CState-object (m, AU).
   */
    double SemiMinorAxis(const BodyID planet, const double mjd) {
      return getBody(planet)->SemiMinorAxis(mjd);
    };


    /// Returns a vector containing the cartesian position and velocity at the given date.
   /**
    The mentioned reference system is a right handed cartesian system with an origin depending
    on the reference body. Usually this reference body is the sun and therefore the reference
    system is often a heliocentric, eclliptical one with the x-axis pointing towards the vernal
    equinox. But it mustn't, the earth's natural satellite, the moon, for example uses a
    geocentric reference frame.
    \param planet The ID of the desired body in solar system.
    \param mjd The date for which the data should be valid for.
    \return A vector with the cartesian position and velocity (m, AU).
   */
    Vector Cart(const BodyID planet, const double mjd) {
      return getBody(planet)->Cart(mjd);
    };



   /**
   * \brief State variables with reference to a certain body, date and unit type.
   *
   * The function returns a vector containing the three position and velocity
   * components of a body in solar system indicated by the parameter "target".
   * The reference frame of the returned cartesian coordinates is centered at
   * a body indicated by the parameter "center". Usually this will be the J2000
   * reference frame. The returned values are guilty for a point of time given
   * by the parameter "date".
   *
   * \param target The "target" body whose coordinates and velocities are
   * desired.
   * \param center The "center" body and therefor the origin of the coodinate
   * system.
   * \param date Point of time for which the returned data should be valid for.
   * \param unit_type The desired type of units the vector should be returned in.
   * \return A Vector containing the desired position and velocity information.
   * \sa example cart.cpp
    */
    Vector& Cart(const BodyID target, const BodyID center,
                const double date, const UnitType unit_type);


/**
* \brief State variables with reference to a certain frame, date and unit type.
*
* The Vector returned by this member function consists of cartesian position
* and velocity information. The are with respect to a certain body, date (epoch)
* and reference frame (coordinate system) which might even have a different
* epoch. Basically the function tries to calculate the geometric difference
* vector from the body "center" to "target". If this has been performed
* successfully, it tries to convert it into the desired unit system and
* reference frame.
*
* \param target ID of target body
* \param center ID of center body
* \param date Epoch of target and center bodies
* \param unit_type Unit system of returned state content
* \param ref_frame Desired reference frame to express the result in
* \param ref_frame_date Epoch of reference frame
*/
    Vector Cart(const BodyID target,
                const BodyID center,
                const double date,
                const UnitType unit_type,
                const CoordSysType ref_frame,
                const double ref_frame_date);

   /**
   * \brief Returns a vector containing the position and velocity of given
   * date in a polar coordinates.
   *
   * The mentioned reference system is a right handed cartesian system with a
   * origin depending on the reference body. Usually this reference body is the
   * sun and therefor the reference system is often a heliocentric, eclliptical
   * one. But it mustn't, the earth's natural satellite, the moon, for example
   *  uses the earth as reference system.
   *
   * \param planet The ID of the desired body in solar system.
   * \param mjd The date for which the data should be valid for.
   * \return A vector containing the polar position and velocity (m, AU, rad).
   */
    Vector Polar(const BodyID planet, const double mjd) {
      return getBody(planet)->Polar(mjd);
    };



   /**
   * \brief Returns a vector containing the Keplerian elements of the body's
   * orbit (mean anomaly).
   *
   * The mentioned reference system is a right handed cartesian system with a
   * origin depending on the reference body. Usually this reference body is the
   * sun and therefor the reference system is often a heliocentric, eclliptical
   * one. But it mustn't, the earth's natural satellite, the moon, for example
   * uses the earth as reference system.
   *
   * \param planet The ID of the desired body in solar system.
   * \param mjd The date for which the data should be valid for.
   * \return A vector containing the Keplerian elements (with mean anomaly).
   */
    Vector KeplerM(const BodyID planet, const double mjd) {
      return getBody(planet)->KeplerM(mjd);
    };



   /**
   * \brief Returns a vector containing the Keplerian elements of the body's
   * orbit (true anomaly).
   *
   * The mentioned reference system is a right handed cartesian system with a
   * origin depending on the reference body. Usually this reference body is the
   * sun and therefor the reference system is often a heliocentric, eclliptical
   * one. But it mustn't, the earth's natural satellite, the moon, for example
   * uses the earth as reference system.
   *
   * \param planet The ID of the desired body in solar system.
   * \param mjd The date for which the data should be valid for.
   * \return A vector containing the Keplerian elements (with true anomaly).
   */
    Vector Keplerf(const BodyID planet, const double mjd){
      return getBody(planet)->Keplerf(mjd);
    };


    /// Unit vector in radial direction (orbit reference frame).
    Vector rVector(const BodyID planet, const double mjd){
      return getBody(planet)->rVector(mjd);
    };


    /// Unit vector in transversal direction (orbit reference frame).
    Vector tVector(const BodyID planet, const double mjd){
      return getBody(planet)->tVector(mjd);
    };


    /// Unit vector in direction of angular momentum (perpendicular to orbit).
    Vector hVector(const BodyID planet, const double mjd){
      return getBody(planet)->hVector(mjd);
    };


    /// The GAUSSIAN P-Unit Vector of the body's orbit of
    /// date (points to periapsis).
    Vector PVector(const BodyID planet, const double mjd){
      return getBody(planet)->PVector(mjd);
    };


    /// The GAUSSIAN Q-Unit Vector of the body's orbit
    /// of date (RxP).
    Vector QVector(const BodyID planet, const double mjd){
      return getBody(planet)->QVector(mjd);
    };


    /// The GAUSSIAN R-Unit Vector of the body's orbit of
    /// date (normal to orbit plane).
    Vector RVector(const BodyID planet, const double mjd){
      return getBody(planet)->RVector(mjd);
    };

    // hill radius and sphere of influence

    /** Tbd. */
    double RHill(const BodyID planet, const double mjd){
      return getBody(planet)->RHill(mjd);
    };

    /** Tbd. */
    double RSOI(const BodyID planet, const double mjd){
      return getBody(planet)->RSOI(mjd);
    };

    /**
    * \brief The enumeration that holds the possible return values of
    * 'GravityAssistPossible'
    *
    * The abbreviation 'GAC' means Gravity Assist Check. The different return
    * values of the method 'GravityAssistPossible' indicate what the calling
    * function could do in order to improve the solution or to get closer to
    * the sphere of influence.
    */
    enum GACResult {
      SUCCESS,                  ///< manuever successfully performed, step size can be increased
      FAILURE_NO_EPS_DIST,      ///< no epsilon for the difference between SOI-radius and distance from planet
      OUTSIDE_SOA,              ///< to far away from body, therefore time step size decrease is needed
      INSIDE_SOA_MOVING_AWAY,   ///< repeat the last integration step with a smaller time step size
      INSIDE_SOA_MOVING_TOWARDS,  ///< continue integration with the last time step size
      INSIDE_SOI                ///< repeat last integration step with a smaller time step size
    };


    /**! \brief Check method for gravity assist maneuvers.
       This function is intended to be called after every integration step to
       check if a given spacecraft (s/c) is close enough to a certain body
       to become dominated by the body's gravitational force. This is true
       if the s/c is within the so called sphere of influence.
       The return value describes in more detail, what further modifications
       should be done to the parameter 'step size' in the integration process.
      */
    GACResult GravityAssistPossible(const BodyID planet,  ///< index of center body for gravity assist
      CState & sc,              ///< state of the s/c, that should do the gravity assist
      double &dT,               ///< time step size used for last integration step
      const double epsDist      ///< maximum distance from sphere of influence as a fraction of the SOI radius
      );


    /// Calculates (numerical integration) a gravity assist maneuver
    /// with respect to a certain body of Solar System
    bool GravityAssistNumerical(const BodyID planet,  ///< index of center body for gravity assist
      const double end_time,    ///< maximum time, the gravity assist could last
      CState & dest             ///< state of the s/c, that should do the gravity assist
      );

    /// Calculates (analytical solution) a gravity assist maneuver
    /// with respect to certain body in Solar System
    bool GravityAssistAnalytical(const BodyID planet, ///< index of center body for gravity assist
      CState& dest,             ///< state of the s/c, that should do the gravity assist
      double i,                 ///< 1st gravity assist parameter
      double j                  ///< 2nd gravity assist parameter
      );


    /// Calculates (analytical solution) a gravity assist maneuver
    /// with respect to certain body in Solar System
    bool GravityAssistAnalytical(const BodyID planet, ///< index of center body for gravity assist
      CState & dest             ///< state of the s/c, that should do the gravity assist
      );


    /**
    * \brief Calculates the distance to another celestial bodies or orbital position
    *
    * The method returns the distance of the desired body to another body
    * incorporated by a state class object. The date on which that distance
    * estimation should be done will be taken from the given parameter object
    * named body.
    *
    * \param planet The index of the body/planet in solar system that should used
    * \param target_state Estimates whether a free orbit or the position of
    * the given boy will be used for comparison
    * \param a Semimajor axis of the free orbit (depending on target_state)
    * \param e eccentricity of the free orbit (depending on target_state)
    * \param i inclination of the free orbit (depending on target_state)
    * \param Omega longitude of ascending node of the free orbit (depending on target_state)
    * \param omega argument of periapsis of the free orbit (depending on target_state)
    * \param body reference to a CState class object that describes the body whose distance
    * to the solar system object or free orbit should be estimated.
    * \return The distance between the solar system object/free orbit and the given body in AU.
    */
    double Distance(const BodyID planet,  // mjd is body's
                    const StateType target_state,
                    const double a,
                    const double e,
                    const double i,
                    const double Omega,
                    const double omega,
                    const CState & body);

    /** Tbd. */
    double Distance_km(const BodyID planet, // mjd is body's
                       const StateType target_state,
                       const double a,
                       const double e,
                       const double i,
                       const double Omega,
                       const double omega,
                       const CState & body);

    /** Tbd. */
    double Velocity(const BodyID planet,  // mjd is body's
                    const StateType target_state,
                    const double a,
                    const double e,
                    const double i,
                    const double Omega,
                    const double omega,
                    const CState & body);

    /** Tbd. */
    double Velocity_kms(const BodyID planet,  // mjd is body's
                        const StateType target_state,
                        const double a,
                        const double e,
                        const double i,
                        const double Omega,
                        const double omega,
                        const CState & body);

    /** Tbd. */
    double SensitivityFactor(const BodyID planet,
                             const int n,
                             const double dr_max,
                             const double dv_max)
    {

#ifdef DEBUG
      cout << "CSolarSystem::SensitivityFactor(const BodyID, const int, const double, const double)"
           << " => entered!\n";
      cout.flush();
#endif

      return getBody(planet)->SensitivityFactor(n, dr_max, dv_max);
    }

    // input / output

    /** Tbd. */
    void PrintHeader(const BodyID planet,
                     std::ofstream & out,
                     const char *comment = NULL) const {

      getBody(planet)->PrintHeader(out, comment);

    }

    /** Tbd. */
    void Print(const BodyID planet,
               const double mjd,
              std::ofstream & out) {

      getBody(planet)->Print(mjd, out);
    }


    /// Retuns ID of body if it is in solar system. Otherwise 0.
    BodyID getPlanetIndex(const string & planet_name) const;

  };                            // end class CSolarSystem

}                               // end namespace SpaceLib


#endif // CSOLARSYSTEM_H

