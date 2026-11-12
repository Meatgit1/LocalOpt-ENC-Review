// file     : CStateData.h
// author   : Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
// compiler : GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++ 
// date     : 15.08.2006 (AO)
// purpose  : Definition class CState's internal data structure class


#ifndef CSTATEDATA_H
#define CSTATEDATA_H

#include <CState.h>
#include <CAttitude.h>
#include <MLVector.h>

namespace SpaceLib {

  class CSolarSystem;

  //----------------------------------------------------------------------------
  struct CoordinateState 
  {
    MathLib::Vector position;
    MathLib::Vector velocity;

    CoordinateState(void):position(3), velocity(3) {};

    const struct CoordinateState &operator=(const struct CoordinateState &_cs) 
    {
      position = _cs.position;
      velocity = _cs.velocity;
      return (*this);
    }

    bool operator==(const struct CoordinateState& rhs)
    {
      return ((position == rhs.position) &&
        (velocity == rhs.velocity));
    }
  }; // CoordinateState


  //----------------------------------------------------------------------------
  struct CStateData 
  {
  public:
    // administration data
    bool useDE405;              // a switch, whether JPL's DE405 should be used, or not
    DE405StateIndex indexDE405; // if DE405 is used, this is the index in the database
    CStateType type;             // determines what kind our state is made of(real or artificial)
    CSolarSystem *refSys;       // reference system of state (optional)                          // should be removed
    double GM_AUD;              // the gravitational parameter of the state (its mass) in AUd units [AU^3/d^2]
    double GM_SI;               // and in SI units [m^3/s^2]                                     // should be removed
    CState *m_central_state;    // pointer to a CState - object of the center body
    BodyID m_central_body;      // central body for which the state is valid                     // should be removed
    UnitType unitType;          // used unit type system (AUd units or SI units )

    // physical data
    double m_mjd;               // modified julian date of epoch the state's components are valid for

    struct CoordinateState cartesian; ///< position/velocity expressed (cartesian coordinates)
    struct CoordinateState sphere;    ///< position/velocity expressed (sphere coordinates)

    CAttitude attitude;               ///< attitude
 
    MathLib::Vector angular_velocity; ///< rotational state (angular velocity)

    double m_a;                 ///< semi-major axis
    double m_e;                 ///< eccentricity
    double m_i;                 ///< inclination
    double m_Omega;             ///< longitude of the ascending node
    double m_omega;             ///< argument of perihelion
    double m_curlypi;           ///< longitude of perihelion
    double m_M;                 ///< mean / mean hyperbolic anomaly
    double m_E;                 ///< eccentric / hyperbolic anomaly
    double m_f;                 ///< true anomaly
    double m_L;                 ///< mean longitude

    //--------------------------------------------------------------------------
    CStateData(void)
      : useDE405(false)
       ,indexDE405(DE405_NO_BODY)
       ,type(BODY)
       ,refSys(0)
       ,GM_AUD(0.0)
       ,GM_SI(0.0)
       ,m_central_state(0)
       ,m_central_body(SS_BARY)
       ,unitType(AUD)
       ,m_mjd(MJD_J2000)
       ,cartesian()
       ,sphere()
       ,attitude()
       ,angular_velocity(3)
       ,m_a(0.0)
       ,m_e(0.0)
       ,m_i(0.0)
       ,m_Omega(0.0)
       ,m_omega(0.0)
       ,m_curlypi(0.0)
       ,m_M(0.0)
       ,m_E(0.0)
       ,m_f(0.0)
       ,m_L(0.0)
    {}

    //--------------------------------------------------------------------------
    ~CStateData() 
    {
      refSys=0;
      m_central_state=0;
    }

    //--------------------------------------------------------------------------
    CStateData(const CStateData& to_copy)
      : useDE405(to_copy.useDE405)
      ,indexDE405(to_copy.indexDE405)
      ,type(to_copy.type)
      ,refSys(to_copy.refSys)
      ,GM_AUD(to_copy.GM_AUD)
      ,GM_SI(to_copy.GM_SI)
      ,m_central_state(to_copy.m_central_state)
      ,m_central_body(to_copy.m_central_body)
      ,unitType(to_copy.unitType)
      ,m_mjd(to_copy.m_mjd)
      ,cartesian(to_copy.cartesian)
      ,sphere(to_copy.sphere)
      ,attitude(to_copy.attitude)
      ,angular_velocity(to_copy.angular_velocity)
      ,m_a(to_copy.m_a)
      ,m_e(to_copy.m_e)
      ,m_i(to_copy.m_i)
      ,m_Omega(to_copy.m_Omega)
      ,m_omega(to_copy.m_omega)
      ,m_curlypi(to_copy.m_curlypi)
      ,m_M(to_copy.m_M)
      ,m_E(to_copy.m_E)
      ,m_f(to_copy.m_f)
      ,m_L(to_copy.m_L)
    {}

    //--------------------------------------------------------------------------
    const CStateData &operator=(const CStateData &_cs) 
    {
      useDE405 = _cs.useDE405;
      indexDE405 = _cs.indexDE405;
      type = _cs.type;
      refSys = _cs.refSys;
      GM_AUD = _cs.GM_AUD;
      GM_SI = _cs.GM_SI;
      m_central_state = _cs.m_central_state;
      m_central_body = _cs.m_central_body;
      unitType = _cs.unitType;
      m_mjd = _cs.m_mjd;

      cartesian = _cs.cartesian;
      sphere = _cs.sphere;
      attitude = _cs.attitude;
      angular_velocity = _cs.angular_velocity;

      m_a = _cs.m_a;
      m_e = _cs.m_e;
      m_i = _cs.m_i;
      m_Omega = _cs.m_Omega;
      m_omega = _cs.m_omega;
      m_curlypi = _cs.m_curlypi;
      m_M = _cs.m_M;
      m_E = _cs.m_E;
      m_f = _cs.m_f;
      m_L = _cs.m_L;

      return (*this);
    }   // operator=

    //--------------------------------------------------------------------------
    bool operator==(const CStateData& rhs)
    {
      return(
        (useDE405 == rhs.useDE405) &&
        (indexDE405 == rhs.indexDE405) &&
        (type == rhs.type) &&
        (refSys == rhs.refSys) &&
        (GM_AUD == rhs.GM_AUD) &&
        (GM_SI == rhs.GM_SI) &&
        (m_central_state == rhs.m_central_state) &&
        (m_central_body == rhs.m_central_body) &&
        (unitType == rhs.unitType) &&
        (m_mjd == rhs.m_mjd) &&
  
        (cartesian == rhs.cartesian) &&
        (sphere == rhs.sphere) &&
        (attitude == rhs.attitude) &&
        (angular_velocity == rhs.angular_velocity)&&

        (m_a == rhs.m_a) &&
        (m_e == rhs.m_e) &&
        (m_i == rhs.m_i) &&
        (m_Omega == rhs.m_Omega) &&
        (m_omega == rhs.m_omega) &&
        (m_curlypi == rhs.m_curlypi) &&
        (m_M == rhs.m_M) &&
        (m_E == rhs.m_E) &&
        (m_f == rhs.m_f) &&
        (m_L == rhs.m_L)
      );
    };

  };                            // struct CStateData

};                               // namespace SpaceLib

#endif // #ifndef

// end of file

