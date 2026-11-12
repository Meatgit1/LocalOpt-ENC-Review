/*!
\file CBody.h
\version 1.7
\brief A class for a celestial body (CBody).
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Tanja D. Schmidt <tanja.schmidt@dlr.de>c (TS)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 06.11.2006 (AO)

\b Purpose:
   This module represents the functionality of a class of bodies
   in space that actually have a mass whose attracting force is
   driving them arounf a center body that also has a mass. The
   planets motion describing elements can be of different kinds.
   Either you use a polynomial series of Keplerian elements in
   order to use osculating elements, or you use the simple
   elements of the undisturbed motion of the planet. But you can
   also use the undisturbed simple Keplerian elements by omitting
   the series elements of a power higher than one. You can also
   provide initial conditions of a distance an a velocity vector
   out of which also non osculating Keplerian elements will be
   computed and used further. In connection with CSolarSystem class
   it is also possible to use a ephemeris data file from JPL/NASA
   which offers exact position and velocity informations of the
   solar systems planets and the earth's satellite moon.

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
   (c) 2001-2005 Bernd Dachwald, Tanja D. Schmidt, Andreas Ohndorf
*/

#ifndef CPLANET_H
#define CPLANET_H

#include <MLConst.h>
#include <MLGlobal.h>
#include <CState.h>   // base class

using namespace std;
using namespace MathLib;

namespace SpaceLib 
{
  class CBodyData;                      // forward declaration of data structure
  class  DE405;

  /// A class for planets and other celestial bodies.
  class CBody : virtual public CState
  {
    friend class CSolarSystem; // forward declaration of friend class

  protected:

    CBodyData* planetData;         ///< data structure of additionally body data

    /// pointer to internal data structure
    CBodyData* getData(void);

    /// pointer to internal data structure (const version)
    const CBodyData* getData(void) const;

    /** Tbd. */
    void Initialize(
      const BodyID index,
      const string a_name,
      const double mass,
      const double epoch,
      const double a,
      const double e,
      const double i,
      const double Omega,
      const double curlypi,
      const double L,
      const double da,
      const double de,
      const double di,
      const double dOmega,
      const double dcurlypi,
      const double dL,
      const CoordType _CSt = ELEMENTS,
      CState * _CS = 0);


    /// get the state information from DE405 data file
    Vector DE405_State(
      DE405StateIndex _index,                // body index in DE405 file
      double _mjd                            // modified julian date of epoch
    );


    /// initialization with the DE405 ephemeris database
    void InitializeDE405(
      const BodyID index,
      const BodyID center,
      CState * aCS,
      const DE405StateIndex indexDE405,
      const string a_name,
      CSolarSystem * _ref,
      CStateType _type
      ,DE405* _DE = 0
    );

    /** Sets the state with position and velocity informations gained from DE405. */
    void setStateDE405(const double mjd);

    /** Sets the state with actual kepler elements gained from polynom development. */
    void setStateElements(const double mjd);

    /** sets the state according to the given epoch. */
    virtual void SetState(const double mjd);

    //--------------------------------------------------------------------------
    // CBody CLASS INTERFACE
    //--------------------------------------------------------------------------
  public:

    static void initStatics(void);

    CBody(void);

    /** Copy constructor */
    CBody(const CBody &);

    /** Assign operator */
    CBody & operator=(const CBody &);

    /** Compare operator.*/
    bool operator==(const CBody& rhs);

    /** destructor */
    virtual ~CBody(void) 
    {
      SafeDelete(planetData);
    }

     /*!
       \brief Constructor to initialize a body in solar system.

       This is the standard way to construct and initialize a space
       object that has a mass.

       \param a_name The name of the space object.
       \param mass Body's mass.
       \param epoch Reference data for osculating Keplerian elements, if used.
       \param a Semimajor axis
       \param e Eccentricity
       \param i Inclination
       \param Omega Longitude of ascending node
       \param curlypi Longitude of pericenter (Omega+omega)
       \param L Mean Longitude = curlypi + M (mean anomaly)
       \param da increment of a per day
       \param de increment of e per day
       \param di increment of i per day
       \param dOmega increment of Omega per day
       \param dcurlypi increment of curlypi per day
       \param dL increment of mean longitude per day
       \param _CSt type of coordinates entered in other parameters
       \param _CS Pointer to center state.

     */
     CBody(const string a_name,
             const double mass,
             const double epoch,
             const double a,
             const double e,
             const double i,
             const double Omega,
             const double curlypi,
             const double L,
             const double da,
             const double de,
             const double di,
             const double dOmega,
             const double dcurlypi,
             const double dL,
             const CoordType _CSt,
             CState * _CS
         );

    /// Changes the object's unique ID. !!Make sure that you know what you do!!
    /// Usually, you don't need this function.
    void setId(const BodyID id);

    /// unique body ID within a CSolarSystem object
    BodyID getId(void) const;

    /// Returns the body's name (must be unique).
    const std::string& getName(void) const;

    /*! \brief set Keplerian elements (a,e,i, Omega, Omega+omega, Omega+omega+M)
     *
     * The user provided input structure must have an index operator [] and
     * at least six data values. There is no scope checking implemented so
     * be sure to give valid data structures.\n
     * Currently there's only an implementation for Mathlib::Vector structures.
     *  \brief param input data structure holding the needed six parameters
     *  \brief epoch modified julian date (MJD) of parameters epoch
     *  \sa setKeplerianElementsRate*/
    template<class T>
    void setKeplerianElements(const T& input, const double epoch);

    /*! \brief set Keplerian elements rate (da/dt, de/dt,di/dt, dOmega/dt, d(Omega+omega)/dt,
    d(Omega+omega+M/dt) )
     *
     * Given rates have to be expressed as change rate per day and
     * the user provided input structure must have an index operator [] and
     * at least six data values. There is no scope checking implemented so
     * be sure to give valid data structures.\n
     * Currently there's only an implementation for Mathlib::Vector structures.
     *  \brief param input data structure holding the needed six parameters
     *  \brief epoch modified julian date (MJD) of parameters epoch
     *  \sa setKeplerianElements*/
    template<class T>
    void setKeplerianElementsRate(const T& input, const double epoch);

    /// Sets the body's name.
    void setName(const string &);

    /// Returns a reference of the CState object of the body.
    CState & State(const double mjd);

    // Cartesian x-coordinate of the planet.
    /** \param mjd The desired date of x - coordinate.
    */
    double x(const double mjd); // x-coordinate


    // Cartesian y-coordinate of the planet.
    /** \param mjd The desired date of y - coordinate.
    */
    double y(const double mjd); // y-coordinate


    // Cartesian z-coordinate of the planet.
    /** \param mjd The desired date of z - coordinate.
    */
    double z(const double mjd); // z-coordinate

    // Cartesian x-component of the planet's velocity.
    /** \param mjd The desired date of the velocities x-component.
    */
    double dx(const double mjd);  // x-velocity


    // Cartesian y-component of the planet's velocity.
    /** \param mjd The desired date of the velocities y-component.
    */
    double dy(const double mjd);  // y-velocity


    // Cartesian z-component of the planet's velocity.
    /** \param mjd The desired date of the velocities z-component.
    */
    double dz(const double mjd);  // z-velocity


    // Distance from the origin of coordinate system
    /** This is the first component of the position vector described in
        polar coordinates. It is actually the norm or absolute value
        of the distance vector.
        \param mjd The desired date the returned distance should be valid for.
    */
    double r(const double mjd); // radius


    // Polar angle phi.
    /** This is the second component of the position vector described in
       polar coordinate system. It is the azimuthal angle between
       the x-axis and the projection of the distance vector on the xy - plane.
       Its values ranges from 0 to 2*Pi.
       \param mjd The desired date the returned angle should be valid for.
    */
    double phi(const double mjd); // azimuth


    // Polar angle theta.
    /** This is the third component of the position vector described in
       polar coordinate sytem. It is the elevation angle between the base plane
       (xy-plane) and the position vector. It ranges from -Pi/2 to +Pi/2, where
       +Pi/2 matches the z-direction of the cooordinate system and an angle of
       0 would be in the base plane.
       \param mjd The desired date the returned angle should be valid for.
    */
    double theta(const double mjd); // elevation


    // The planets radial velocity.
    /** This is the first derivation of the distance r with respect to time.
        \param mjd The desired date the returned value should be valid for.
    */
    double dr(const double mjd);


    /// Angular velocity dphi (=dPhi/dt) for a given date.
    double dphi(const double mjd);


    /// Angular velocity dtheta (=dTheta/dt) for a given date.
    double dtheta(const double mjd);


    /// The semimajor axis of the planets orbit for a given date.
    double a(const double mjd);


    /// The eccentricity of the planets orbit for a given date.
    double e(const double mjd);


    /// The inclination of the planets orbit for a given date.
    double i(const double mjd);


    /// The longitude of the ascending node of the planets orbit for a given date.
    double Omega(const double mjd);


    /// The argument of periapsis of the planets orbit for a given date.
    double omega(const double mjd);


    /// The longitude of periapsis of the planets orbit for a given date.
    double curlypi(const double mjd);


    /// The mean anomaly of the planets orbit for a given date.
    double M(const double mjd);


    /// The eccentric anomaly of the planets orbit for a given date.
    double E(const double mjd);


    /// The true anomaly of the planets orbit for a given date.
    double f(const double mjd);


    /// The mean longitude of the planets orbit for a given date.
    double L(const double mjd);

    /// Tbd.
    double h(const double mjd);

    /// Tbd.
    double k(const double mjd);

    /// Tbd.
    double p(const double mjd);

    /// Tbd.
    double q(const double mjd);

    /// Tbd.
    double l(const double mjd);

    /// The absolute velocity of the planet for a given date.
    double Velocity(const double mjd);


    /// The planets mass specific angular momentum at a given date.
    double Momentum(const double mjd);


    /// The planets mass specific kinetic energy at a given date.
    double KineticEnergy(const double mjd);


    /// The planets mass specific potential energy at a given date.
    double PotentialEnergy(const double mjd);


    /// The planets mass specific total energy at a given date.
    double Energy(const double mjd);


    /// The orbit period of the planet according to the state data at the desired date.
    double OrbitPeriod(const double mjd);


    /// The mean motion of the planet at the desired date.
    double MeanMotion(const double mjd);


    /// The semiminor axis calculated by the  values of position and velocity at the given date.
    double SemiMinorAxis(const double mjd);


    /// The cartesian elements of the position and velocity vectors at the given date.
    Vector Cart(const double mjd);


    /// Returns the position and the velocity at the given time in polar coordinates.
    Vector Polar(const double mjd);

    /**
    * \brief Returns state vector of desired Lagrange point with respect to
    * central body.
    *
    * Guilty values for parameter LPoint are in range of [1,5]. The
    * Lagrangian equilibrium points are enumerated as follows:\n 1 is the point
    * between the center body 1 and the smaller body 2. Point 2 is the one
    * outside of body 2 on the line connecting the central body and body 2.
    * Point 3 is the corresponding one outside of the central body. These three
    * points are known as the colinear equilibrium points. The remaining two
    * triangular points are point 4 and 5 whereas point 4 is the leading one
    * and 5 is the trailing one. Leading in this context means, it seems that
    * point 4 is "before" body 2.\n
    * The resulting vector contains position and velocity information of desired
    * Lagrangian point in the coordinate system of body 1. In general, this
    * is a ecliptical one with a unit system of AUD.
    *
    * \param LPoint The desired Lagrange equilibrium point.
    * \param mjd Date of body 2 and thus the Lagrangian point.
    * \return A state vector with position/velocity data of Langrangian point.
    * In case of no existing central body or error during calculation, a zero
    *  vector will be returned.
    */
    Vector LagrangePoint(const int LPoint, const double mjd);


    /// Returns a vector with the Keplerian elements of the planets orbit at the given time where
    /// the position in the orbit is described by the mean anomaly.
    Vector KeplerM(const double mjd);


    /// Returns a vector with the Keplerian elements of the planets orbit at the given time where
    /// the position in the orbit is described by the true anomaly.
    Vector Keplerf(const double mjd);


    /// Returns a unit vector pointing in the direction of the distance vector at the given date.
    Vector rVector(const double mjd);


    /// Returns a unit vector laying in the orbit plane. It is the cross product of the
    /// normalized distance and angular momemtum vectors.
    Vector tVector(const double mjd);


    /// Returns a unit vector pointing in the direction of the angular momentum vector.
    /// It is perpendicular to the orbit plane and the distance vector.
    Vector hVector(const double mjd);


    /// Returns the Gaussian P unit vector that points to periapsis.
    Vector PVector(const double mjd);


    /// Returns the Gaussian Q unit vector (R x P).
    Vector QVector(const double mjd);


    /// Returns the Gaussian R unit vector that is normal to the planets orbit.
    Vector RVector(const double mjd);


    /// Returns the Hill-Radius of the planet at the given date.
    double RHill(const double mjd);


    /// Returns the radius of the planet's sphere of influence at a given date.
    double RSOI(const double mjd);


    /**
    * \brief Returns the distance to another space object (in AU)
    *
        Calculates the distance to another object whose orbit information
        is stored in the given CState reference. According to the parameter
        target_state this distance will a real distance to a real body or
        a free orbit in space. In this case it is not known for what true
        anomaly and therefore what position in orbit of the artificial "body"
        this distance should be calculated. The true anomaly will be
        calculated so that the projection of the resulting orbit position
        vector is in line with the base plane projection of the "body"
        position vector. The date of calculation will be the date of body.
        \param body The body the distance to should be calculated
        \param target_state Determines whether to use a real body or "artificial" one
        \return The distance between the target state and the method's CState body/orbit.
        \sa Distance_km
    */
    double Distance(const CState & body, const StateType target_state);



    /**
    * \brief Returns the distance to another space object (in KM)
    *
        Calculates the distance to another object whose orbit information
        is stored in the given CState reference. According to the parameter
        target_state this distance will a real distance to a real body or
        a free orbit in space. In this case it is not known for what true
        anomaly and therefore what position in orbit of the artificial "body"
        this distance should be calculated. The true anomaly will be
        calculated so that the projection of the resulting orbit position
        vector is in line with the base plane projection of the "body"
        position vector. The date of calculation will be the date of body.
        \param body The body the distance to should be calculated
        \param target_state Determines whether to use a real body or "artficial" one
        \return The distance between the target state and the method's CState body/orbit.
        \sa Distance
    */
    double Distance_km(const CState & body, const StateType target_state);



    /**
    * \brief Returns the absolute velocity to another space object (in AU/D)
        Calculates the velocity relative to another object whose orbit information
        is stored in the given CState reference. According to the parameter
        target_state this velocity will be a real velocity to a real body or
        a free orbit in space. In this case it is not known for what true
        anomaly and therefore what position in orbit of the artificial "body"
        this vecocity should be calculated. The true anomaly will be
        calculated so that the projection of the resulting orbit position
        vector is in line with the base plane projection of the "body"
        position vector. The date of calculation will be the date of body.
        \param body The body the relative velocity to should be calculated to
        \param target_state Determines whether to use a real body or "artificial" one
        \return The absolute velocity between the target state and the method's CState body/orbit.
        \sa Velocity_kms
    */
    double Velocity(const CState & body, const StateType target_state);



    /**
    * \brief Returns the absolute velocity to another space object (in KM/S)
        Calculates the velocity relative to another object whose orbit information
        is stored in the given CState reference. According to the parameter
        target_state this velocity will be a real velocity to a real body or
        a free orbit in space. In this case it is not known for what true
        anomaly and therefore what position in orbit of the artificial "body"
        this vecocity should be calculated. The true anomaly will be
        calculated so that the projection of the resulting orbit position
        vector is in line with the base plane projection of the "body"
        position vector. The date of calculation will be the date of body.
        \param body The body the relative velocity to should be calculated to
        \param target_state Determines whether to use a real body or "artificial" one
        \return The absolute velocity between the target state and the method's CState body/orbit.
        \sa Velocity
    */
    double Velocity_kms(const CState & body, const StateType target_state);

    /** Tbd. */
    double SensitivityFactor( const int n, const double& dr_max, const double& dv_max);

    /// Prints the position, velocity and some other data for the given date
    /// to an output stream.
    void Print(const double mjd, std::ofstream & out);

    /// print VRML script code to output stream describing the body's orbit
    /// with color given in components and given # of steps
    void printVrmlTrajectory(std::ostream& os, const int width
      , const int prec, const double time, const int n_steps
      , const double red, const double green, const double blue);

  }; // class CBody

  //----------------------------------------------------------------------------
  inline CBodyData* CBody::getData(void) 
  {
    return planetData;
  }

  //----------------------------------------------------------------------------
  inline const CBodyData* CBody::getData(void) const 
  { 
    return planetData;
  }

}                               // namespace SpaceLib

#endif                          // CPLANET_H

