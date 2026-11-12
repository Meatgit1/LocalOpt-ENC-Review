/*!
\file CState.h
\version 1.7
\brief A class for an astrodynamical state (CState).
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Tanja D. Schmidt <tanja.schmidt@dlr.de>c (TS)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 25.08.2006 (AO)

\b Purpose:
   The declaration of the elemental base class CState is done
   in this header file. CState contains the basic information
   of a body/object in the used space. The more specialised
   class CBody is derived from CState.

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

#ifndef CSTATE_H
#define CSTATE_H

// STL
#include <iostream>
#include <cassert>

// SpaceLib
#include <Kepler.h> // supplementary functions for calculation of orbit elements
#include <DE405.h>  // DE 405 access functions and constants


/*! \brief type of CState class object
 * It indicates, which kind of state is incorperated by the class CState. The
 * main difference ist the fact, that there are special points of interest
 * which have no mass and regular space bodies with a mass.*/
enum CStateType
{
  SPECIAL, ///< e.g. barycenter, Lagrange equilibrium point (no mass)
  BODY     ///< e.g. planets, SUN, spacecraft, asteroid, comet (mass)
};


using namespace MathLib;

namespace SpaceLib 
{

  const double E_M_SMALL = 1.0E-9;  ///< lowest eccentric mean anomaly
  const double H_M_SMALL = 1.0E-9;  ///< lowest hyperbolic mean anomaly

  // forward declarations
  struct CStateData;            // internal data structure
  class CSolarSystem;           // solar system object class
  class CAttitude;              // attitude class

  // ---------------------------------------------------------------------------
  /// base class of SpaceLib celestial bodies
  class CState
  {
  protected:

    CStateData* stateData; ///< internal data

    /** Tbd. */
    void Set(const double v0 = 0.0,  // parameter 1 (x,r,a)
      const double v1 = 0.0,         // parameter 2 (y,phi,e)
      const double v2 = 0.0,         // parameter 3 (z,theta,i)
      const double v3 = 0.0,         // parameter 4 (dx,dr,Omega)
      const double v4 = 0.0,         // parameter 5 (dy,dphi,omega)
      const double v5 = 0.0,         // parameter 6 (dz,dtheta,M,f,L)
      const CoordType st = CARTESIAN,// type of parameters
      const double mjd = MJD_J2000,  // date of data
      const BodyID cb = SUN,         // index of center body
      const UnitType ut = AUD)throw (std::runtime_error);      // unit type of given data

    /** \brief Cartesian to polar state. */
    void Cart2Polar(void);

    /** \brief Polar to cartesian state. */
    void Polar2Cart(void);

    /** \brief Cartesian to keplerian state. */
    void Cart2Kepler(void);

    /** \brief Keplerian to cartesian state. */
    void Kepler2Cart(void);

// to be removed
    
    /** GM of central body in AUd units (old function). */
    double GM_AUd(void) const;

    /** GM of central body in SI units (old function). */
    double GM_SI(void) const;
    
    /** GM from stored values (old function). */
    double GM_old(void) const;

    /** Tbd. */
    double E_PARABOLIC1(void) const;

    /** Tbd. */
    double E_PARABOLIC2(void) const;

    /** Tbd. */
    double H_PARABOLIC1(void) const;

    /** Tbd. */
    double H_PARABOLIC2(void) const;

    /** Tbd. */
    void useDE405(bool);

//-----------------------------------------------------------------------------
// CState CLASS INTERFACE
//-----------------------------------------------------------------------------
  public:

    /// default constructor

   /** \hideinitializer
       Initializes a state in the given unit system at the given time,
       whose center body is the given one (default one is sun). The coordinates
       will be cartesian and set to zero.
    */
    CState(const double mjd = MJD_J2000,  ///< the modified julian date of state
           const BodyID cb = SUN,         ///< index of center state
           const UnitType ut = AUD        ///< the unit system to use
    );


    CState(const CState & toCopy);  ///< Copy constructor


    /// Vector constructor
   /** Initializes a state with the given parameters stored in a Vector, with the
       given parameter type system (default is CARTESIAN), at the given time (default is
       01.01.2000 12:00 TBD), with the given center body (default is the sun) and the
       unit type system  (unit type is AUD).
    */
    CState(
      /// A vector of the state parameters to initialize the CState with
      /** These can be either be cartesian coordinates (x,y,z,dx,dy,dz),
          sphere coordinates (r,phi,theta,dr,dphi,dtheta) or Keplerian
          elements (a,e,i,Omega,omega,M[E,L]).
      */
      const Vector & V,
      const CoordType st = CARTESIAN, ///< the desired type of elements
      const double mjd = MJD_J2000, ///< the epoch of state
      const BodyID cb = SUN,    ///< center body
      const UnitType ut = AUD   ///< the type of units to use
      );


    /// Array constructor
    /**
      Initializes a state with the given data in a memory area described by the pointer.
      As in the Vector constructor, the coordinate system type, the epoch, the center body
      and the used unit type system can be specified or ommited. Then the given default
      values will be used.
    */
    CState(
      /// State parameters are stored in the memory area pointed to by this pointer.
      /** These can be either be cartesian data (x,y,z,dx,dy,dz), sphere data
          (r,phi,theta,dr,dphi,dtheta) or keplerian data (a,e,i,Omega,omega,M[E,L]).
      */
      const double *p, const CoordType st = CARTESIAN,  ///< the desired type of elements
      const double mjd = MJD_J2000, ///< the epoch of state
      const BodyID cb = SUN,    ///< the index of center body
      const UnitType ut = AUD   ///< the type of units to use
      );


    /// Value constructor
    /**
      Initializes a state with the given data in double values.
      As in the Vector constructor, the coordinate system type, the epoch, the index
      of the center body and the used unit type system can be specified or ommited.
      Then the given default values will be used.
    */
    CState(
      const double v0,          ///< parameter 1 (x,r,a)
      const double v1,          ///< parameter 2 (y,phi,e)
      const double v2,          ///< parameter 3 (z,theta,i)
      const double v3,          ///< parameter 4 (dx,dr,Omega)
      const double v4,          ///< parameter 5 (dy,dphi,omega)
      const double v5,          ///< parameter 6 (dz,dtheta,M,f,L)
      const CoordType st = CARTESIAN, ///< the desired type of elements
      const double mjd = MJD_J2000, ///< the epoch of state
      const BodyID cb = SUN,    ///< the index of center body
      const UnitType ut = AUD   ///< the type of units to use
      );


    const CState& operator=(const CState& toCopy);          ///< Assign operator
    bool operator==(const CState& rhs);                      ///< equal operator
    ~CState(void);                                               ///< destructor

    /// sets the CState object with the data from the given CState object.
    /** The object given by a reference will be copied. Therefor our object
        represents an exact copy of the argument.

        \param s the reference to a CState object.
    */
    virtual void SetState(const CState & s) {(*this)=s;};


    /// attitude expressed as CAttitude class object
    const CAttitude& getAttitude(void) const;

    /// angular velocity
    const MathLib::Vector& CState::getAngularVelocity(void) const;


    /**
    * \brief Sets the CState's data with a vector.
    * This is the one of the different setState functions that covers the
      need for an initializer when you've got a vector (see MathLib)
      containing the data. A call of this method updates all the class'
      data fields.

      \param V          a reference to a vector class (six elements), containing the data
      \param st         The kind of elements used in the vector.
      \param a_mjd      The epoch, for which the data is valid.
      \param cb         the index of center body
      \param ut         the type of units to use
      \param aST        the type of state (a real body in space or not)
      \param useDE405   this integer estimates, whether the data comes from
                        DE405(JPL) or not
      \param indexDE405 the body's in DE405 database
      \param aCS        a pointer to body's center CState object
    */
    void SetState(
      const Vector & V,
      const CoordType st = CARTESIAN,
      const double a_mjd = MJD_J2000,
      const BodyID cb = SUN,
      const UnitType ut = AUD,
      const CStateType aST = BODY,
      const bool useDE405 = false,
      const DE405StateIndex indexDE405 = DE405_SUN,
      CState * aCS = 0
    );


    /**
    * \brief Sets the CState's data with the data from a pointer
    * This is one of the different setState functions that covers the
    * need for an initializer when you've got an array of double values
    * containing the data for describing an orbit. The call of this
    * method updates all the class' data fields.\n
    * Attention: It is important to set the right central body which means
    * the index of the central body must match to a real body (a major
    * planet for example). The reason for this is obvious if you try to set
    * the CState with Keplerian elements. To calculate the Cartesian velocities
    * SpaceLib needs the gravitational parameter of the central body. For
    * the major planets the values are "hardwired" in SpaceLib.

      \param p          a pointer to an array of double values (six elements)
                        containing the data
      \param st         The kind of elements used in the vector.
      \param a_mjd      The epoch, for which the data is valid.
      \param cb         the index of center body
      \param ut         the type of units to use
    */
    void SetState(
      const double *p,
      const CoordType st = CARTESIAN,
      const double a_mjd = MJD_J2000,
      const BodyID cb = SUN,
      const UnitType ut = AUD
    );



    /// Set the CState's data with the given data of six double values

  /** This is the one of different setState functions that covers the
      need for an initializer when you've got the data for describing an orbit
      as double values. Every call of this method updates all the data fields of
      the class.
   */
    void SetState(
      const double v0,          ///< parameter 1 (x,r,a)
      const double v1,          ///< parameter 2 (y,phi,e)
      const double v2,          ///< parameter 3 (z,theta,i)
      const double v3,          ///< parameter 4 (dx,dr,Omega)
      const double v4,          ///< parameter 5 (dy,dphi,omega)
      const double v5,          ///< parameter 6 (dz,dtheta,M,f,L)
      const CoordType st = CARTESIAN, ///< the type of the given elements
      const double mjd = MJD_J2000, ///< the epoch, the data is valid for
      const BodyID cb = SUN,    ///< the index of the center body
      const UnitType ut = AUD   ///< the unit type of the given data
      );



    // accessor methods

    /// reference to CAttitude object representing current attitude
    CAttitude& attitude(void);

    /*!
    \brief Stores the given pointer of the CSolarSystem object for which the
    CState object is valid in a private variable.

    The given pointer is a CSolarSystem pointer and thus a reference to
    the environment the CState object is located in.
    \param solsys A pointer to a CSolarSystem object.
    \sa getRefSolSys
    */
    void setRefSolSys(CSolarSystem* solsys);


    /*!
    \brief Returns the stored pointer to a CSolarSystem object.

    The returned pointer is the CSolarSystem object for which the CState object
    is valid for. That means, that the the information held within this CState
    object is \em not suitable for using in another CSolarSystem object.
    \return Pointer to a CSolarSystem object
    \sa setRefSolSys
    */
    CSolarSystem* getRefSolSys(void) const;


    /// Returns whether the DE405 data file is used (true) or not (false).
    bool useDE405(void);

    /// Returns the index of the state's body in DE405 data file.

  /** If the JPL's DE405 data file is used and the state incorporates one of these
      bodies, the index in the data file can be obtained with this method.
      \return The desired index or DE405_NO_BODY.
      \sa useDE405
  */
    DE405StateIndex getIndexDE405(void);


    /// Sets the type of the CState object.

  /** The type of a CState object can be one of the values of CStateType. So you should
      look at the possible values of this enumeration.
      \sa CStateType
      \return nothing
   */
    void setType(CStateType);    // sets the type (see def. of CStateType)



    /// Gets the type of the CState object.
  /** The type of a CState object can be one of the values of CStateType. So you should
      look at the possible values of this enumeration.
      \sa CStateType
      \return The type of CState as CStateType value.
   */
    CStateType getType(void) const;  // returns it



    /// Stores the given poiner to an CSolarSystem object in a local variable.

  /** The stored CSolarSystem-pointer can be used to obtain the data of other
      bodies in a star system.
      \returns nothing
   */
//    void setRefSys(CSolarSystem *);


    /// Returns the pointer to a reference CSolarSystem object.

  /** The stored CSolarSystem-pointer can be used to obtain the data of other
      bodies in a star system.
      \returns the stored pointer to the CSolarSystem object used as reference.
   */
    //  CSolarSystem *getRefSys(void);



   /// Stores a pointer to a CState object incorporating the center object.
   /** Every state normally describes an orbit around a central body. In "SpaceLib"
       this central body is also represented by a CState object. The pointer to such
       an central state is stored in every "child" state. The purpose of this method
       is to pass the mentioned pointer to the object.
       \return nothing
    */
    void setCenterState(CState *);  // stores the CState of center



    /// Returns a pointer to a CState object incorporating the center object.

  /** Every state normally describes an orbit around a central body. In "SpaceLib"
      this central body is also represented by a CState object. The pointer to such
      an central state is stored in every "child" state. The purpose of this method
      is to return the mentioned pointer to the object.
      \return Pointer to the CState obejct, describing the central body
   */
    CState *getCenterState(void); // pointer to center CState-object



    /// sets the value of the state's gravitational parameter ("its mass")

  /** The gravitational parameter mentioned here is the product of the gravitational
      constant gamma and the mass of the body and thus the most important parameter
      influencing orbits of other bodies who are influenced by this state.

      \param aGM The gravitational parameter to be set for the state
      \param aType the unit type system, the value of GM is given in

      \return nothing
   */
    void setGM(double aGM, UnitType aType); // sets the gravitational parameter



    /// returns the value of the state's gravitational parameter ("its mass")

  /** The gravitational parameter mentioned here is the product of the gravitational
      constant gamma and the body's mass and thus the most important parameter,
      which influences orbits of other bodies around this state.

      \param aType the unit type system, the value of GM should be returned in

      \return the stored gravitational parameter of THIS body, not the orbit determining
              gravitational parameter of the central body
      \sa G_SI,G_AUd
   */
    double getGM(UnitType aType) const; // returns it



    /// Sets the index of the central body.

  /** With this index and a CSolarSystem index, this value is valid for, you
      can access this body and all the data he stores. This method stores the
      given index in the class object.

      \param aPI The index of the central body in the star system.

      \return nothing
   */
    void setCenterBodyId(const BodyID aPI);



    /// Gets the index of the central body.

  /** With this index and a CSolarSystem index, this value is valid for, you
      can access this body and all the data he stores. This method returns the
      stored index in the class object.

      \return The index(as BodyID) of the central body of the state.
   */
    BodyID getCenterBodyId(void) const;



    /// Returns the unit system type the data in a CState object is stored in.

  /**
     \return The used unit type system as an UnitType value.
     \sa UnitType
   */
    UnitType Unit(void) const;


    /// Returns the epoch for that the data in a CState is valid.
    /** Naturally some values in a CState object are only guilty for a
        specific moment of time called the epoch. Such data values are
        for example the cartesian and spheric coordinates and velocities.
        This method returns the mentioned epoch as modified julian date.

       \return The epoch of CState data as modified julian date.
    */
    double MJD(void) const;


//    Vector VecFromReferenceKOS(const BodyID reference_body,
//                               CSolarSystem & ss,
//                               CoordSysType _KOS = HEN) const;


    /// returns the cartesian X-coordinate of the state

  /** This coordinate is referenced to the center state, which must
      not be the central star of a star system. For example the coordinates
      of the state incorporating the moon would probably be geocentric.
      In order to get the heliocentric coordinates, you've got to add the
      heliocentric coordinates of the earth. Therefor the coordinate accessor
      methods have a variation, that is able to provide the desired data relative
      to an arbitrary body in a CSolarSystem object.

      \return the X-coordinate in an cartesian coordinate system
   */
    double x(void) const;

    /// returns the cartesian X-coordinate of the state, relative to another body

  /** This coordinate is referenced to the given body in the also given CSolarSystem
      object.

      \param reference_body The index of the body in a CSolarSystem object
                            that should be used as reference
      \param star_system A reference to the CSolarSystem object which should
                         include both bodies (this state an the reference_body
      \param _KOS The type of coordinate system that should be used
                  (heliocentric, bodycentric,...).
      \return the X-coordinate in an cartesian coordinate system relative to
              the given body
   */
//  double x(const BodyID    reference_body,
//                 CSolarSystem & star_system,
//                 CoordSysType   _KOS = HEN) const;


    /// returns the cartesian Y-coordinate of the state

  /** This coordinate is referenced to the center state, which must
      not be the central star of a star system. For example the coordinates
      of the state incorporating the moon would probably be geocentric.
      In order to get the heliocentric coordinates, you've got to add the
      heliocentric coordinates of the earth. Therefor the coordinate accessor
      methods have a variation, that is able to provide the desired data relative
      to an arbitrary body in a CSolarSystem object.

      \return the Y-coordinate in an cartesian coordinate system
   */
    double y(void) const;

    /// returns the cartesian Y-coordinate of the state, relative to another body

  /** This coordinate is referenced to the given body in the also given CSolarSystem
      object.

      \param reference_body The index of the body in a CSolarSystem object that
                            should be used as reference
      \param star_system A reference to the CSolarSystem object which should
                         include both bodies (this state an the reference_body
      \param _KOS The type of coordinate system that should be used
                  (heliocentric, bodycentric,...).
      \return the Y-coordinate in an cartesian coordinate system relative to the given body
   */
//  double y(const BodyID    reference_body,
//                 CSolarSystem & star_system,
//                 CoordSysType   _KOS = HEN) const;


    /// returns the cartesian Z-coordinate of the state

  /** This coordinate is referenced to the center state, which must
      not be the central star of a star system. For example the coordinates
      of the state incorporating the moon would probably be geocentric.
      In order to get the heliocentric coordinates, you've got to add the
      heliocentric coordinates of the earth. Therefor the coordinate accessor
      methods have a variation, that is able to provide the desired data relative
      to an arbitrary body in a CSolarSystem object.

      \return the Z-coordinate in an cartesian coordinate system
   */
    double z(void) const;


    /// returns the cartesian Z-coordinate of the state, relative to another body

  /** This coordinate is referenced to the given body in the also given CSolarSystem
      object.

      \param reference_body The index of the body in a CSolarSystem object, that should be used as reference
      \param star_system A reference to the CSolarSystem object, which should include both bodies (this state an the reference_body
      \param _KOS The type of coordinate system, that should be used (heliocentric, bodycentric,...).
      \return the Z-coordinate in an cartesian coordinate system relative to the given body
   */
//  double z(const BodyID    reference_body,
//           CSolarSystem & star_system,
//           CoordSysType   _KOS = HEN) const;



    /// returns the cartesian X-component of body's velocity

  /** This velocity is referenced to the center state, which must
      not be the central star of a star system. For example the velocity component
      of the state incorporating the moon would probably be geocentric.
      In order to get the moon's heliocentric velocity, you've got to add the
      heliocentric velocity of the earth. Therefor the velocity accessor
      methods have a variation, that is able to provide the desired data relative
      to an arbitrary body in a CSolarSystem object.

      \return X-component of body's velocity in an cartesian coordinate system
   */
    double dx(void) const;

    /// returns the cartesian X-component of body's velocity relative to another body

  /** This velocity is referenced to the given body in the also given CSolarSystem
      object.

      \param reference_body The index of the body in a CSolarSystem object, that should be used as reference
      \param star_system A reference to the CSolarSystem object, which should include both bodies (this state an the reference_body
      \param _KOS The type of coordinate system, that should be used (heliocentric, bodycentric,...).
      \return X-component of body's velocity in an cartesian coordinate system
              relative to the given body
   */
//  double dx(const BodyID    reference_body,
//                  CSolarSystem & star_system,
//            CoordSysType   _KOS = HEN) const;


    /// returns the cartesian Y-component of body's velocity

  /** This velocity is referenced to the center state, which must
      not be the central star of a star system. For example the velocity component
      of the state incorporating the moon would probably be geocentric.
      In order to get the moon's heliocentric velocity, you've got to add the
      heliocentric velocity of the earth. Therefor the velocity accessor
      methods have a variation, that is able to provide the desired data relative
      to an arbitrary body in a CSolarSystem object.

      \return Y-component of body's velocity in an cartesian coordinate system
   */
    double dy(void) const;


    /// returns the cartesian Y-component of body's velocity relative to another body

  /** This velocity is referenced to the given body in the also given CSolarSystem
      object.

      \param reference_body The index of the body in a CSolarSystem object, that should be used as reference
      \param star_system A reference to the CSolarSystem object, which should include both bodies (this state an the reference_body
      \param _KOS The type of coordinate system, that should be used (heliocentric, bodycentric,...).
      \return Y-component of body's velocity in an cartesian coordinate system
              relative to the given body
   */
//  double dy(const BodyID    reference_body,
//                  CSolarSystem & star_system,
//                  CoordSysType   _KOS = HEN) const;


    /// returns the cartesian Z-component of body's velocity

  /** This velocity is referenced to the center state, which must
      not be the central star of a star system. For example the velocity component
      of the state incorporating the moon would probably be geocentric.
      In order to get the moon's heliocentric velocity, you've got to add the
      heliocentric velocity of the earth. Therefor the velocity accessor
      methods have a variation, that is able to provide the desired data relative
      to an arbitrary body in a CSolarSystem object.

      \return Z-component of body's velocity in an cartesian coordinate system
   */
    double dz(void) const;


    /// returns the cartesian Z-component of body's velocity relative to another body

  /** This velocity is referenced to the given body in the also given CSolarSystem
      object.

      \param reference_body The index of the body in a CSolarSystem object that
                            should be used as reference
      \param star_system A reference to the CSolarSystem object which should
                         include both bodies (this state an the reference_body
      \param _KOS The type of coordinate system, that should be used
                  (heliocentric, bodycentric,...).
      \return Z-component of body's velocity in an cartesian coordinate system
              relative to the given body
   */
//  double dz(const BodyID    reference_body,
//                  CSolarSystem & star_system,
//                  CoordSysType   _KOS = HEN) const;


    /// The distance from the origin of reference system (polar coordinates).

  /** This function is rather simple and does nothing more than to return
      the stored value of the distance to the origin of reference coordinate
      system.

      \return the distance value with respect to the reference system
   */
    double r(void) const;

    /// Value of azimuth angle (polar coordinates).

  /** In a polar coordinate system this is the azimuth angle. It is
      measured in the xy-plane, starts from the x-axis counterclockwise
      with respect to the positive z-axis to the projection of the desired vector
      on the xy-plane and can take values from 0 to 2*PI.

      \return Azimuth angle phi
   */
    double phi(void) const;

    /// Returns the Value of elevation angle theta (polar coordinates).

  /** In a polar coordinate system this is the elevation angle. It
      is measured from the xy-plane to the desired vector (positive values
      in direction of positive z-axis) and takes values from (-Pi/2,+Pi/2].
    */
    double theta(void) const;

    /// Returns the value of radial velocity (polar coordinates).

  /** This coordinate is defined in the cartesian coordinate system and
      describes the radial velocity of the state, that means the
      change of the 'length' of the position vector per time unit.
    */
    double dr(void) const;

    /// Returns the angular velocity Phi dot.

  /** This velocity element is bound to a cartesian coordinate system
      and describes the change in azimuth angle per time unit.
    */
    double dphi(void) const;

    /// Returns the angular velocity Theta dot.

  /** This velocity element is bound to a cartesian coordinate system
      and describes the change in elevation angle per time unit.
   */
    double dtheta(void) const;

    /// Returns the semi major axis a of date [AU,m].
    double a(void) const;

    /// Returns the value of eccentricity e of date.
    double e(void) const;

    /// Returns the value of the inclination i of date [rad].
    double i(void) const;

    /// Returns the longitude of the ascending node of date [rad].
    double Omega(void) const;

    /// Returns the argument of pericenter of date [rad].
    double omega(void) const;

    /// Returns the longitude of pericenter = Omega + omega [rad].
    double curlypi(void) const;

    /// Returns mean anomaly [rad].
    double M(void) const;

    /// Returns the eccentric anomaly [rad].
    double E(void) const;

    /// Returns the true anomaly [rad].
    double f(void) const;

    /// longitude L = Omega + omega + M
    double L(void) const;

    /// velocity vectors radial component
    double vr(void) const;


    /// velocity vectors azimuthal component
    double vphi(void) const;


    /// velocity vectors elevational component
    double vtheta(void) const;

//!
    double h(void) const; ///< equinoctial element
    double k(void) const; ///< equinoctial element
    double p(void) const; ///< equinoctial element
    double q(void) const; ///< equinoctial element
    double l(void) const; ///< equinoctial element
//!

    /// Tbd.
    double Momentum(void) const;

    /// Returns the value of the semi major axis [AU, m].
    double SemiMinorAxis(void) const;

    /// Returns the absolute value of the state's velocity [AU/D, m/s].
    double Velocity(void) const;

    /// Returns the mass specific kinetic energy [AU^2/D^2, m^2/s^2].
    double KineticEnergy(void) const;

    /// Returns the mass specific potential energy [AU^2/D^2, m^2/s^2].
    double PotentialEnergy(void) const;

    /// Returns the mass specific total energy [AU^2/D^2, m^2/s^2].
    double Energy(void) const
    {
      return KineticEnergy() + PotentialEnergy();
    };

    /// Returns the amount of time, the state needs for one orbit [D, s].
    double OrbitPeriod(void) const;

    /// Returns the value of the fictional uniform velocity the state has [rad/D, rad/s].

  /** For the case of a circular orbit, the mean motion is equal to the true motion.
      The mean motion is the velocity the state would have, if it moves uniformly
      along the trajectory of the orbit. Because of Kepler's second law it can be shown
      that for an elliptical orbit the angular velocity is not uniform for the
      complete orbit. It is very easy, to calculate the mean anomaly with the mean
      motion and the time since passage of pericenter and therefor the true anomaly.
    */
    double MeanMotion(void) const;


    /// Returns the value of the semi latus rectum [AU, m].
    double SemilatusRectum(void) const;


    /// Returns the rotation angle Zeta [rad].

  /** Zeta is the rotation angle between the orbit reference frame
      (E_r,E_t,E_h) and the reference frame of sphere
      coordinates (E_r,E_phi,E_teta) [rad].
    */
    double Zeta(void) const;


    /// Returns a vector containing the cartesian components of position and velocity.

  /** This velocity element is bound to a cartesian coordinate system
      and contains the elements of the position and velocity vector.
      They describe the state of motion of the body or spacecraft in space.
   */
    Vector Cart(UnitType unit_type = STORED) const;


    /// Returns a vector containing the polar coordinates components of position and velocity.

  /** This vector is bound to a polar coordinates system
      and contains the elements of the position and angular velocities.\n
      The elements of the vector are therefor:\n
      r -----> distance from origin\n
      phi ---> azimuth angle\n
      theta -> elevation angle\n
      dr ----> first derivative of distance -> radial velocity\n
      dphi --> first derivative of azimuth angle ->azimuthal angular velocity\n
      dtheta > first derivative of elevation angle ->elevational angular velocity\n

      \sa Cart, Polarv, KeplerM, Keplerf

   */
    Vector Polar(void) const;
    // vector containing the state's position and velocity information
    // (r, phi, theta, dr, dphi, dtheta) in sphere coordinate elements



    /// polar coordinates components of position and angular velocities.

  /** This returned vector contains the polar coordinates of the
      position vector and the components of the velocity vector expressed
      in the (r,phi,theta) reference frame.
      r -----> distance from origin\n
      phi ---> azimuth angle\n
      theta -> elevation angle\n
      vr ----> first derivative of vector r = radial velocity\n
      vphi --> first derivative of vector r in direction of phi- unit vector\n
      vtheta > first derivative of vector r in direction of theta- unit vector\n

      \sa Cart, Polar, KeplerM, Keplerf
      \return a vector containing the state's position and velocity elements in sphere coordinates

   */
    Vector Polarv(void) const;

    /**
    * \brief vector with the ceplerian elements (mean anomaly)
    * The returned vector contains the classic Keplerian elements
      whereas the bodiy's position in the orbit plane is described by the
      parameter "mean anomaly". The resulting vector is:\n
        a --->semi major axis\n
        e --->eccentricity\n
        i --->inclination of orbit\n
        Omega longitude of ascending node\n
        omega argument of pericenter\n
        M --->mean anomaly
      \sa Cart, Polar, Polarv, Keplerf
    */
    Vector KeplerM(void) const;


    /// vector with ceplerian elements (true anomaly)
    Vector Keplerf(void) const;


    /// Returns the unit vector along the radius vector.

  /** This is the first unit vector of the orbit reference frame. This frame
      is the result of a positive x-rotation of the sphere-coordinates
      reference frame with the angle zeta.

      \sa Zeta, tVector, hVector
    */
    Vector rVector(void) const;


    /// Returns the unit vector along the velocity vector.

  /** This is the second unit vector of the orbit reference frame. This frame
      is the result of a positive x-rotation of the sphere-coordinates
      reference frame with the angle zeta.

      \sa Zeta, rVector, hVector
    */
    Vector tVector(void) const;



    /// Returns the unit vector along the vector of angular momentum.

  /** This is the third unit vector of the orbit reference frame. This frame
      is the result of a positive x-rotation of the sphere-coordinates
      reference frame with the angle zeta.

      \sa Zeta, rVector, tVector
    */
    Vector hVector(void) const;



    /// Gaussian unit vector P (points to periapsis).
    /** The returned vector's elements are in units of a cartesian reference frame.
    */
    Vector PVector(void) const;



    /// Gaussian unit vector Q (= R x P).

  /** The returned vector's elements are in units of a cartesian reference frame.
    */
    Vector QVector(void) const;



    /// Gaussian unit vector R (normal to orbit).

  /** The returned vector's elements are in units of a cartesian reference frame.
    */
    Vector RVector(void) const;


    //------------------------ functionality-----------------------------------

    /// the gravitational parameter of the central body

  /** The method returns the gravitational parameter of the central body
      the state is orbiting around.

      \return G*M = gravitational constant multiplied with the central body mass
   */
    double GM(void) const;



    /// Returns true, if the state's orbit is elliptic.
    bool IsElliptic(void) const;



    /// Returns true, if the state's orbit is parabolic.
    bool IsParabolic(void) const;



    /// Returns true, if the state's orbit is hyperbolic.
    bool IsHyperbolic(void) const;

    /// yes if state has a central body different from the sun, false otherwise
    bool isBodycentric(void) const;

    /// yes, if the sun is the central body
    bool isHeliocentric(void) const;

	   
	//----------------------- input / output------------------------------------

    /// Prints out a header with the column description of the corresponding Print() method.

  /** This is used in conjunction with the Print()-method to produce comma seperated value-
      files (CSV) which can be read by every usual table calculation program. If you
      want to know, what the meaning of the colums is, you've to call this method. It
      prints (comma sperated) the description of the columns. Usually it is followed
      by many calls to Print().

      \sa Print
    */
    void PrintHeader(std::ofstream & out, const char *comment = NULL) const;

    /// Prints the states's data fields to standard output (CSV).

  /** The output is done in the order, that is descroibed by the PrintHeader - method,
      which should be called before, but must not. The values are seperated by a comma,
      one by one. Therefor it is easy to use them with a table calculation program.
    */
    void Print(std::ofstream & out) const;


    /* TestCState debug functions
       friend void PrintStateError(CState &);
       friend void CalcTestCStateError(CState &, CState &, const double);
       friend void     xyzError(CState &, CState &);
       friend void    dxyzError(CState &, CState &);
       friend void       rError(CState &, CState &);
       friend void     phiError(CState &, CState &);
       friend void   thetaError(CState &, CState &);
       friend void      drError(CState &, CState &);
       friend void    dphiError(CState &, CState &);
       friend void  dthetaError(CState &, CState &);
       friend void       aError(CState &, CState &);
       friend void       eError(CState &, CState &);
       friend void       iError(CState &, CState &);
       friend void   OmegaError(CState &, CState &);
       friend void   omegaError(CState &, CState &);
       friend void curlypiError(CState &, CState &);
       friend void       MError(CState &, CState &);
       friend void       EError(CState &, CState &);
       friend void       fError(CState &, CState &);
       friend void       LError(CState &, CState &); */



  };                            // end class CState

  // ----------------------------------------------------------------------------
  /*! \brief helper function for calculation of maximum rates of change of 
  Keplerian orbital elements depending on current state
  *
  * references.: 
  * [1] Anastassios E. Petropoulos, "Low-Thrust Orbit Transfers using Candidate
  * Lyapunov Functions with a Mechanism for Coasting", 2004, AIAA 2004-5089
  *\n
  * Calculates the maximum achievable rates of change desired Keplerian element
  * with formulaes given in [1]. The returned value must be multiplied
  * with the s/c current absolute acceleration in order to obtain the final value.
  *
  * \param element the Keplerian element for which you want to calculate the
  maximum rate of change
  * \param s current orbit described with CState class object
  * \return desired normed maximum orbital element rate of change (multiply with s/c acceleration)
  * \sa CalcMaxElementRateOfChange, CalcMaxPeriapsisRateOfChange
  * \sa CalcMaxSemimajorAxisRateOfChange, CalcMaxEccentricityRateOfChange
  * \sa CalcMaxInclinationRateOfChange*/
  inline double CalcMaxElementRateOfChange(
    const KEPLERIAN_ELEMENT element,
    CState& s)
  {
    switch (element)
    {
    case SEMIMAJORAXIS:
      return CalcMaxSemimajorAxisRateOfChange(s.e(), s.a(), (s.getCenterState())->GM());
      break;
    case ECCENTRICITY:
      return CalcMaxEccentricityRateOfChange(s.SemilatusRectum(), s.Momentum());
      break;
    case INCLINATION:
      return CalcMaxInclinationRateOfChange(s.e(), s.Momentum(), s.omega(), 
        s.SemilatusRectum());
      break;
    case NODE:
      return CalcMaxNodeRateOfChange(s.e(), s.Momentum(), s.i(), s.omega()
        , s.SemilatusRectum());
      break;
    case PERIAPSIS:
      return CalcMaxPeriapsisRateOfChange(s.e(), s.Momentum(), s.SemilatusRectum());
      break;
    default: assert(false); return 0.0; break; // raise assertation failed exception for debugging
    }
  }

}                               // namespace SpaceLib

#endif                          // CSTATE_H

