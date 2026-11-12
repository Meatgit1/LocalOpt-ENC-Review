/*!
\file CBody.cpp
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
   (c) 2001-2006 Bernd Dachwald, Tanja D. Schmidt, Andreas Ohndorf
*/

//SpaceLib
#include <configSpaceLib.h>

// C/C++ standard lib
#include <math.h>
#include <cassert>

// STL
#include <string>
#include <sstream>
#include <algorithm>

// SpaceLib
#include <MLMatrix.h>
#include <Kepler.h>
#include <CStateData.h>
#include <CBodyData.h>
#include <CBody.h>

using namespace std;
using namespace Global;
using namespace MathLib;


namespace SpaceLib 
{

  static bool factors_are_calculated = false;
  static double factor_max[3] = { 0.0, 0.0, 0.0 };
  static double factor_min[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
  static double factor_avg[3] = { 0.0, 0.0, 0.0 };


  void CBody::initStatics(void) {

    factors_are_calculated = false;

    for (int i=0; i < sizeof(factor_max)/sizeof(double);++i) 
    {
      factor_max[i] = factor_avg[i] = 0.0;
      factor_min[i] = HUGE_VAL;
    }
    
  }

  // destructor
  // --------------------------------------------------------------------------
  CBody::CBody(void):CState(), planetData(0) {

    planetData = new CBodyData;
    if (planetData == 0)
      FatalError("CBody:CBody -> error with memory allocation");

  }


  // copy constructor
  // --------------------------------------------------------------------------
  CBody::CBody(const CBody & to_copy)
    :CState(to_copy)
    , planetData(new CBodyData(*(to_copy.planetData)))
  {
      if (planetData == 0)
        FatalError("CBody:CBody(const CBody&) -> error during memory \
                    allocation (CBodyData)");
  }


  // assign operator
  // --------------------------------------------------------------------------
  CBody & CBody::operator=(const CBody & rhs) 
  {
    if (&rhs != this)
    {
      CState::operator=(rhs);    
      *planetData = *(rhs.planetData);
    }
    return (*this);
  }
    
  // --------------------------------------------------------------------------
  bool CBody::operator==(const CBody& rhs)
  {
    return *planetData == *(rhs.planetData);
  }


  //---------------------------------------------------------------------------
  // note:
  //
  //   this function initializes all members according to their epoch of
  //   osculation
  //
  //---------------------------------------------------------------------------
  void CBody::Initialize(const BodyID index,
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
    const CoordType _CSt, CState * _CS) {

    planetData->m_index = index;
    planetData->m_name = a_name;
    planetData->m_epoch = epoch;

    // set state gravitational parameter
    CState::setGM(G_SI * mass, SI);
    // store the pointer to the centers CState object
    CState::setCenterState(_CS);

    planetData->m_elements = Vector(a, e, i, Omega, curlypi, L);
    planetData->m_delements = Vector(da, de, di, dOmega, dcurlypi, dL);

    CState::SetState(planetData->m_elements, _CSt,
                     planetData->m_epoch,
                     SUN,
                     AUD,
                     BODY,
                     0,
                     DE405_SUN,
                     _CS);
  };



  //---------------------------------------------------------------------------
  CBody::CBody(const string a_name,
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
    CState * _CS) 
    : CState() , planetData(0)
  {

    planetData = new CBodyData;
    if (planetData == 0)
      FatalError("CBody:CBody(const CBody&) -> error with memory \
                  allocation (CBodyData)");

    planetData->m_name = a_name;
    upper(planetData->m_name);
    planetData->m_epoch = epoch;

    // set state gravitational parameter
    CState::setGM(G_SI * mass, SI);
    // store the pointer to the centers CState object
    CState::setCenterState(_CS);
    
    
    if (_CSt==CARTESIAN)
    {
      const MathLib::Vector s(a,e,i,Omega,curlypi,L);
      CState::SetState(s,
                     _CSt,
                     planetData->m_epoch,
                     SUN,
                     AUD,
                     BODY,
                     0,
                     DE405_SUN,
                     _CS);
      planetData->m_elements(0) = CState::a();
      planetData->m_elements(1) = CState::e();
      planetData->m_elements(2) = CState::i();
      planetData->m_elements(3) = CState::Omega();
      planetData->m_elements(4) = CState::omega()+CState::omega();
      planetData->m_elements(5) = CState::omega()+CState::omega()+CState::M();
      planetData->m_delements.Set(0.0);
      planetData->m_delements(5) = CalcMeanMotion(CState::a(),
        _CS->GM());
    }
    else if ( _CSt==KEPLERM)
    {
      planetData->m_elements = Vector(a, e, i, Omega, curlypi, L);
      planetData->m_delements = Vector(da, de, di, dOmega, dcurlypi, dL);

      CState::SetState(planetData->m_elements,
        _CSt,
        planetData->m_epoch,
        SUN,
        AUD,
        BODY,
        0,
        DE405_SUN,
        _CS);
    }
  };




  //---------------------------------------------------------------------------
  Vector CBody::DE405_State(
    DE405StateIndex _index,
    double          _mjd) {

    double EM_RatioGM = planetData->m_DE->getConstant("EMRAT ");
    double _JD = _mjd + JD_MJD_0;  // convert to julian date

    Vector tmp;

    if (_index == DE405_SS_BARY) return Vector(0.0,0.0,0.0,0.0,0.0,0.0);
    else {
     try {tmp = planetData->m_DE->getState(_index, _JD);}
     catch (DE405::TimeNotCoveredByFileEXC) {
       return Vector(0.0,0.0,0.0,0.0,0.0,0.0);
     }
    }

    if (_index== DE405_EARTH) {
      Vector _moon = planetData->m_DE->getState(DE405_MOON, _JD);
      tmp = (tmp    - _moon * 1.0 / (1.0 + EM_RatioGM));
    } else if (_index == DE405_MOON) {
      Vector _earth = planetData->m_DE->getState(DE405_EARTH, _JD);
      tmp = (_earth + tmp   * 1.0 / (1.0 + 1.0 / EM_RatioGM));
    }

    return tmp;
  }



  //---------------------------------------------------------------------------
  void CBody::InitializeDE405(
    const BodyID _index,                // index of actual body
    const BodyID _center,               // index of body we are turning around
    CState * _aCS,                      // state of center body
    const DE405StateIndex _indexDE405,  // index of actual body in DE405
    const string _name,                 // the name of the body
    CSolarSystem * _ref,                // reference to "parent" star system
    CStateType _type
    ,DE405 * _DE
  ) {

    useDE405(true);
    planetData->m_index = _index; // set our index
    planetData->m_name = _name;   // the name of the planet
    stateData->indexDE405 = _indexDE405;  // store DE405 index

    planetData->m_DE = _DE;       // store pointer to access object
    // set the gravitational parameter
    CState::setGM(planetData->m_DE->getGM(_indexDE405), AUD);
    CState::setRefSolSys(_ref);    // store the pointer to a solar system object
    CState::setType(_type);       // store the desired type of state
    CState::setCenterBodyId(_center); // store index of center body
    CState::setCenterState(_aCS); // store the pointer of center object

    CBody::SetState(MJD_J2000); // initialize from DE405
  }



  //---------------------------------------------------------------------------
  void CBody::setStateDE405(const double mjd)
  {

    Vector target = DE405_State(stateData->indexDE405, mjd);
    Vector center;

    // get center planet; valid for moon whose center is planet earth
    if (stateData->m_central_state) {
      center = DE405_State(stateData->m_central_state->getIndexDE405(), mjd);
      target -= center;
    }

    // but the coordinates still have to be transformed into
    // eccliptic reference frame
    // therefor we build up a six-dimensional transformation matrix
    // to change position and velocity information with one operation
    Matrix equatorial_to_ecliptical = R_x(EPS0_J2000 * ARCSEC2RAD);
    // DE405 data are equatorial (reference frame is IRCF), so
    // they must be rotated about the x-axis (angle is J2000 obliquity)

    // compose the six-dimensional transformation matrix
    Matrix transformation(6,6);
    transformation.SetSubMatrix(0, 0, equatorial_to_ecliptical * KM2AU);
    transformation.SetSubMatrix(3, 3, equatorial_to_ecliptical * KMS2AUD);

    // now, we do the transformation and store the new data in CState

    CState::SetState(
      (transformation*target),  // apply conversion/transformation
      CARTESIAN,
      mjd,
      stateData->m_central_body,
      AUD,
      stateData->type,
      true,
      stateData->indexDE405,
      stateData->m_central_state);
  } // void CBody::setStateDE405(const double mjd)



  //---------------------------------------------------------------------------
  void CBody::setStateElements(const double mjd) {
    const double dt = mjd - planetData->m_epoch;
    const double a = planetData->m_elements(0) + dt * planetData->m_delements(0);

    const double e = planetData->m_elements(1) + dt * planetData->m_delements(1);

    AssertEccentricityRange(e);

    bool inclination_was_negative = false;
    double i = planetData->m_elements(2) + dt * planetData->m_delements(2);

    if (i < 0.0) {
      i = fabs(i);
      inclination_was_negative = true;
    }
    AssertInclinationRange(i);

    double Omega = planetData->m_elements(3) + dt * planetData->m_delements(3);

    if (inclination_was_negative) Omega += PI;
    AssertAngularRange2pi(Omega);

    double curlypi = planetData->m_elements(4) + dt * planetData->m_delements(4);
    AssertAngularRange2pi(curlypi);

    double L = planetData->m_elements(5) + dt * planetData->m_delements(5);
    AssertAngularRange2pi(L);

    CState::SetState(a,
                     e,
                     i,
                     Omega,
                     curlypi,
                     L,
                     ELEMENTS,
                     mjd,
                     CState::getCenterBodyId());
  }



  //---------------------------------------------------------------------------
  void CBody::SetState(const double mjd) {
    if (useDE405()) {
      setStateDE405(mjd);
    }
    else {
      setStateElements(mjd);
    }
  };


  //---------------------------------------------------------------------------
  // public functions
  //---------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  BodyID CBody::getId(void) const 
  {
    return planetData->m_index;
  }

  //----------------------------------------------------------------------------
  void CBody::setId(const BodyID id) 
  {
    planetData->m_index = id;
  }

  //----------------------------------------------------------------------------
  /// returns body name
  const std::string& CBody::getName(void) const
  {
    return getData()->m_name;
  } // const std::string& getName(void) const

  //---------------------------------------------------------------------------
  template<>
  void CBody::setKeplerianElements(const MathLib::Vector& input
    , const double epoch)
  {
    const unsigned int n_elements = 6;
    assert(input.size() >= n_elements);
    for (unsigned int i=0;i<n_elements;++i) planetData->m_elements[i]=input[i];
    planetData->m_epoch = epoch;
    stateData->m_mjd    = planetData->m_epoch+1.0;  // enforce setState call
    // on next member function call by setting a current date different from
    // epoch date
  }

  //---------------------------------------------------------------------------
  template<>
  void CBody::setKeplerianElementsRate(const MathLib::Vector& input
    , const double epoch)
  {
    const unsigned int n_elements = 6;
    assert(input.size() >= n_elements);
    for (unsigned int i=0;i<n_elements;++i) planetData->m_delements[i]=input[i];
    planetData->m_epoch = epoch;
    stateData->m_mjd    = planetData->m_epoch+1.0;  // enforce setState call
    // on next member function call by setting a current date different from
    // epoch date
  }

  //---------------------------------------------------------------------------
  /// sets the name of the CBody class
  void CBody::setName(const string & a_name) {
    planetData->m_name = a_name;
  }



  //---------------------------------------------------------------------------
  CState & CBody::State(const double mjd) {
    if (stateData->m_mjd != mjd) SetState(mjd);
    return (*this);
  };



  //---------------------------------------------------------------------------
  double CBody::x(const double mjd) {
    if (stateData->m_mjd != mjd) SetState(mjd);
    return stateData->cartesian.position(0);
  };



  //---------------------------------------------------------------------------
  double CBody::y(const double mjd) {
    if (stateData->m_mjd != mjd) SetState(mjd);
    return stateData->cartesian.position(1);
  };



  //---------------------------------------------------------------------------
  double CBody::z(const double mjd) {
    if (stateData->m_mjd != mjd) SetState(mjd);
    return stateData->cartesian.position(2);
  };



  //---------------------------------------------------------------------------
  double CBody::dx(const double mjd) {
    if (stateData->m_mjd != mjd) SetState(mjd);
    return stateData->cartesian.velocity(0);
  };



  //---------------------------------------------------------------------------
  double CBody::dy(const double mjd) {
    if (stateData->m_mjd != mjd) SetState(mjd);
    return stateData->cartesian.velocity(1);
  };



  //---------------------------------------------------------------------------
  double CBody::dz(const double mjd) {
    if (stateData->m_mjd != mjd) SetState(mjd);
    return stateData->cartesian.velocity(2);
  };



  //---------------------------------------------------------------------------
  double CBody::r(const double mjd) {
    if (stateData->m_mjd != mjd) SetState(mjd);
    return stateData->sphere.position(0);
  };



  //---------------------------------------------------------------------------
  double CBody::phi(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->sphere.position(1);
  };



  //---------------------------------------------------------------------------
  double CBody::theta(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->sphere.position(2);
  };



  //---------------------------------------------------------------------------
  double CBody::dr(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->sphere.velocity(0);
  };

  //---------------------------------------------------------------------------
  double CBody::dphi(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->sphere.velocity(1);
  };



  //---------------------------------------------------------------------------
  double CBody::dtheta(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->sphere.velocity(2);
  };



  //---------------------------------------------------------------------------
  double CBody::a(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->m_a;
  };



  //---------------------------------------------------------------------------
  double CBody::e(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->m_e;
  };



  //---------------------------------------------------------------------------
  double CBody::i(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->m_i;
  };



  //---------------------------------------------------------------------------
  double CBody::Omega(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->m_Omega;
  };




  //---------------------------------------------------------------------------
  double CBody::omega(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->m_omega;
  };



  //---------------------------------------------------------------------------
  double CBody::curlypi(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->m_curlypi;
  };



  //---------------------------------------------------------------------------
  double CBody::M(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->m_M;
  };



  //---------------------------------------------------------------------------
  double CBody::E(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->m_E;
  };



  //---------------------------------------------------------------------------
  double CBody::f(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->m_f;
  };



  //---------------------------------------------------------------------------
  double CBody::L(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return stateData->m_L;
  };



  //---------------------------------------------------------------------------
  double CBody::h(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::h();
  };



  //---------------------------------------------------------------------------
  double CBody::k(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::k();
  };



  //---------------------------------------------------------------------------
  double CBody::p(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::p();
  };



  //---------------------------------------------------------------------------
  double CBody::q(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::q();
  };



  //---------------------------------------------------------------------------
  double CBody::l(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::l();
  };



  //---------------------------------------------------------------------------
  double CBody::Velocity(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::Velocity();
  };



  //---------------------------------------------------------------------------
  double CBody::Momentum(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::Momentum();
  };



  //---------------------------------------------------------------------------
  double CBody::KineticEnergy(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::KineticEnergy();
  };



  //---------------------------------------------------------------------------
  double CBody::PotentialEnergy(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::PotentialEnergy();
  };



  //---------------------------------------------------------------------------
  double CBody::Energy(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::Energy();
  };



  //---------------------------------------------------------------------------
  double CBody::OrbitPeriod(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::OrbitPeriod();
  };



  //---------------------------------------------------------------------------
  double CBody::MeanMotion(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::MeanMotion();
  };



  //---------------------------------------------------------------------------
  double CBody::SemiMinorAxis(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::SemiMinorAxis();
  };



  //---------------------------------------------------------------------------
  Vector CBody::Cart(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::Cart();
  };



  //---------------------------------------------------------------------------
  Vector CBody::Polar(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::Polar();
  };

  //---------------------------------------------------------------------------
  Vector CBody::LagrangePoint(const int LPoint, const double mjd) {

    Vector result(6);
    double GM1,GM2;
    double mu1,mu2;

    if ((LPoint>0) && (LPoint<6)) {

      if (stateData->m_mjd != mjd) SetState(mjd);

      CBody* b1 = dynamic_cast<CBody*>(CState::getCenterState());

      if (b1 != 0) {

        GM1 = b1->CState::getGM(AUD);
        GM2 = CState::getGM(AUD);

        // mu1 and mu2 are the ratios between mass of sun or jupiter and
        // the sum of sun's and jupiters mass
        mu1 = GM1/(GM1+GM2);
        mu2 = GM2/(GM1+GM2);

        const Vector R_b2 = CState::Cart().Slice(0,2);
        const Vector V_b2 = CState::Cart().Slice(3,5);
//        const double nR_b2 = Norm(R_b2);

        const Vector E_r = R_b2/Norm(R_b2);  // unit Vector towards Jupiter (from sun)
        const Vector E_v = V_b2/Norm(V_b2);  // unit vector towards Jupiters velocity
        const Vector E_h = Cross(E_r,E_v);   // orbit normal unit vector

        // the following two vector are unit vectors in Jupiter's
        // orbit plane and will be used to estimate the locations of
        // Langrangian equilibrium points.
        // Ref.:
        //   C.D. Murray, Stanley F. Dermott; Solar System Dynamics,
        //      Cambridge University Press, 1999
        const Vector E_x = E_r;
        const Vector E_y = Cross(E_h,E_r);
        const Vector E_z = E_h;

        double STARTX,STARTY;
        double STARTXp=0.0,STARTYp=0.0;

        double alpha = pow( mu2/(3.0*mu1),(1.0/3.0));
        double beta,ratio;
        double r2;

        switch (LPoint)  {
          case 1:
            // approximation formulae taken from Murray, Dermott [1999]
               r2 =  alpha
                    -(1.0/3.0)  * pow(alpha,2.0)
                    -(1.0/9.0)  * pow(alpha,3.0)
                    -(23.0/81.0)* pow(alpha,4.0);
               STARTX = mu1-r2;
               STARTY = 0.0;
          break;

          case 2:
            // approximation formulae taken from Murray, Dermott [1999]
            r2 =  alpha
                 +(1.0/3.0)  * pow(alpha,2.0)
                 -(1.0/9.0)  * pow(alpha,3.0)
                 -(31.0/81.0)* pow(alpha,4.0);
            STARTX = mu1+r2;
            STARTY = 0.0;
          break;

          case 3:
            // approximation formulae taken from Murray, Dermott [1999]
            ratio = mu2/mu1;
            beta =
              -(7.0/12.0)*ratio
              +(7.0/12.0)*pow(ratio,2.0)
              -(13223.0/20736.0)*pow(ratio,3.0);

            r2 = 2.0+beta;

            STARTX = mu1-r2;
            STARTY = 0.0;
          break;

          case 4:
               STARTX = 0.5-mu2;
               STARTY = sqrt(3.0)/2.0;
          break;

          case 5:
               STARTX = 0.5-mu2;
               STARTY = -sqrt(3.0)/2.0;
          break;
        }

        const Vector R0 = mu2*Norm(R_b2)*E_x;
        const Vector omega = Norm(V_b2)/Norm(R_b2)*E_z;  // angular velocity vector
        const Vector V0 = Cross(omega, R0);

        const Vector r = (STARTX*E_x + STARTY*E_y)*Norm(R_b2);
        const Vector R = (R0 + r);
        const Vector v = (STARTXp*E_x + STARTYp*E_y)*Norm(R_b2)*TWOPI/CState::OrbitPeriod();
        const Vector V = V0 + v + Cross(omega,r);

        // initial position and velocity of test body
        result = Stack(R, V);
      }
    }
    return result;
  }


  //---------------------------------------------------------------------------
  Vector CBody::KeplerM(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::KeplerM();
  };



  //---------------------------------------------------------------------------
  Vector CBody::Keplerf(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::Keplerf();
  };



  //---------------------------------------------------------------------------
  Vector CBody::rVector(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::rVector();
  };



  //---------------------------------------------------------------------------
  Vector CBody::tVector(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::tVector();
  };



  //---------------------------------------------------------------------------
  Vector CBody::hVector(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::hVector();
  };



  //---------------------------------------------------------------------------
  Vector CBody::PVector(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::PVector();
  };



  //---------------------------------------------------------------------------
  Vector CBody::QVector(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::QVector();
  };



  //---------------------------------------------------------------------------
  Vector CBody::RVector(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return CState::RVector();
  };



  //---------------------------------------------------------------------------
  double CBody::RHill(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return r(mjd) * pow(stateData->GM_AUD / (3.0 * GM_SUN_AUD), 1.0 / 3.0);
  };



  //---------------------------------------------------------------------------
  double CBody::RSOI(const double mjd) {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    return r(mjd) * pow(stateData->GM_AUD / GM_SUN_AUD, 0.4);
  };



  //---------------------------------------------------------------------------
  double CBody::Distance(const CState & body,
                          const StateType target_state) {

    const double t = body.MJD();
    double distance;

    if (target_state == BODY_STATE) {
      const Vector planet_r(State(t).Cart(AUD).Slice(0, 2));
      const Vector body_r(body.Cart(AUD).Slice(0, 2));
      distance = Norm(planet_r - body_r);
    }
    else if (target_state == BODY_ORBIT_STATE) {

      CState virt_state(CState::a(),
                        CState::e(),
                        CState::i(),
                        CState::Omega(),
                        CState::omega(),
                        body.phi() - CState::curlypi(),
                        KEPLERF);

      const Vector planet_r(virt_state.x(), virt_state.y(), virt_state.z());
      const Vector body_r(body.Cart(AUD).Slice(0, 2));
      distance = Norm(planet_r - body_r);
    } else {
      FatalError("CBody::Distance(): target state type not valid");
    }
    return distance;
  };



  // --------------------------------------------------------------------------
  double CBody::Distance_km(const CState & body,
                              const StateType target_state) {
    return Distance(body, target_state) * AU2KM;
  };



  // --------------------------------------------------------------------------
  double CBody::Velocity(const CState & body,
                           const StateType target_state) {

    const double t = body.MJD();
    double velocity;

    if (target_state == BODY_STATE) {
      const Vector planet_v(State(t).Cart(AUD).Slice(3, 5));
      const Vector body_v(body.Cart(AUD).Slice(3, 5));

      velocity = Norm(planet_v - body_v);

    } else if (target_state == BODY_ORBIT_STATE) {

      CState virt_state(CState::a()*(stateData->unitType==AUD?1.0:M2AU),
                        CState::e(),
                        CState::i(),
                        CState::Omega(),
                        CState::omega(),
                        body.phi() - CState::curlypi(),
                        KEPLERF);

      const Vector planet_v(virt_state.dx(), virt_state.dy(), virt_state.dz());
      const Vector body_v(body.Cart(AUD).Slice(3, 5));

      velocity = Norm(planet_v - body_v);
    }
    else {
      FatalError("CBody::Velocity(): target state type not valid");
    }
    return velocity;
  };



  //---------------------------------------------------------------------------
  double CBody::Velocity_kms(const CState & body,
                              const StateType target_state) {
    return Velocity(body, target_state) * AUD2KMS;
  };



  //---------------------------------------------------------------------------
  double CBody::SensitivityFactor( const int n, const double& dr_max, const double& dv_max)
  {

#ifdef DEBUG
    cout << " CBody::SensitivityFactor() entered " << endl;cout.flush();
#endif

    if (!(n >= 0 && n <= 2)) FatalError("CBody::SensivityFactor  => index out of range\n");

    if (factors_are_calculated) {
      return factor_max[n];
    }
    else {

#ifdef DEBUG
  cout << " CBody::SensitivityFactor() calculation needed" << endl;cout.flush();
#endif

      cout << endl;
#ifdef DEBUG
  cout << " CBody::SensitivityFactor() adress of dr_max:" << hex << &dr_max
       << endl;cout.flush();
#endif

      ostringstream mess;
      mess << "maximum distance at target         : "
           << dr_max << " km" << endl;
      mess << "maximum relative velocity at target: "
           << dv_max * 1000.0 << " m/s" << endl;
      mess << endl;
      mess << "... calculating sensitivity factors ";
	    Output(mess.str());
	  
      const double T = OrbitPeriod(MJD_J2000)*(stateData->unitType==AUD?1.0:SEC2DAY);
      const double number_of_steps = 360.0;
      const double t_inc(T/number_of_steps);
      const double t_start(MJD_J2000);
      const double t_end(t_start + T);
      const double timespace(t_inc*number_of_steps/10);
      const double fac(stateData->unitType==AUD?1.0:M2AU);

#ifdef DEBUG
  cout << " CBody::SensitivityFactor() constants calculated" << endl;cout.flush();
#endif

      double diff_max;//, diff_min; //AO
      double dr_diff_max, dv_diff_max;
      double dr_diff_plus, dv_diff_plus;
      double dr_diff_minus, dv_diff_minus;
      double _a_, ref_a, ref_e, ref_Omega, ref_i, ref_curlypi, ref_L;

      for (double time = t_start; time < t_end; time += t_inc) {

        if (fmod(time-t_start,timespace)<1.0) Output(string("."));

        SetState(time);

        _a_ = 0.0;
        ref_a = a(time)*fac;
        ref_e = e(time);
        ref_i = i(time);
        ref_Omega = Omega(time);
        ref_curlypi = curlypi(time);
        ref_L = L(time);

        do {
          _a_ += 10.0;

          CState state_plus(ref_a + 1.0 / _a_,
                            ref_e,
                            ref_i,
                            ref_Omega,
                            ref_curlypi,
                            ref_L,
                            ELEMENTS,
                            time,
                            SUN,
                            AUD);
          dr_diff_plus = Distance_km(state_plus, BODY_STATE) / dr_max;
          dv_diff_plus = Velocity_kms(state_plus, BODY_STATE) / dv_max;
          CState state_minus(ref_a - 1.0 / _a_,
                             ref_e,
                             ref_i,
                             ref_Omega,
                             ref_curlypi,
                             ref_L,
                             ELEMENTS,
                             time,
                             SUN,
                             AUD);
          dr_diff_minus = Distance_km(state_minus, BODY_STATE) / dr_max;
          dv_diff_minus = Velocity_kms(state_minus, BODY_STATE) / dv_max;
          dr_diff_max = std::max(dr_diff_plus, dr_diff_minus);
          dv_diff_max = std::max(dv_diff_plus, dv_diff_minus);
          diff_max = std::max(dr_diff_max, dv_diff_max);
//AO          diff_min = min(dr_diff_max, dv_diff_max);

        } while (diff_max >= 1.0);

        factor_avg[0] += _a_;
        if (_a_ > factor_max[0]) factor_max[0] = _a_;
        if (_a_ < factor_min[0]) factor_min[0] = _a_;

        double _e_ = 0.0;

        do {
          _e_ += 10.0;
          double checked_e_plus = bounded(e(time) + 1.0 / _e_, 0.0, 1.0);
          double checked_e_minus = bounded(e(time) - 1.0 / _e_, 0.0, 1.0);

          CState state_plus(ref_a,
                            checked_e_plus,
                            ref_i,
                            ref_Omega,
                            ref_curlypi,
                            ref_L,
                            ELEMENTS,
                            time,
                            SUN,
                            AUD);
          dr_diff_plus = Distance_km(state_plus, BODY_STATE) / dr_max;
          dv_diff_plus = Velocity_kms(state_plus, BODY_STATE) / dv_max;

          CState state_minus(ref_a,
                             checked_e_minus,
                             ref_i,
                             ref_Omega,
                             ref_curlypi,
                             ref_L,
                             ELEMENTS,
                             time,
                             SUN,
                             AUD);
          dr_diff_minus = Distance_km(state_minus, BODY_STATE) / dr_max;
          dv_diff_minus = Velocity_kms(state_minus, BODY_STATE) / dv_max;
          dr_diff_max = std::max(dr_diff_plus, dr_diff_minus);
          dv_diff_max = std::max(dv_diff_plus, dv_diff_minus);
          diff_max = std::max(dr_diff_max, dv_diff_max);
//AO          diff_min = min(dr_diff_max, dv_diff_max);

        } while (diff_max >= 1.0);

        factor_avg[1] += _e_;
        if (_e_ > factor_max[1]) factor_max[1] = _e_;
        if (_e_ < factor_min[1]) factor_min[1] = _e_;

        double _i_ = 0.0;

        do {
          _i_ += 10.0;

          double checked_i_plus = bounded(i(time) + 1.0 / _i_, 1e-9, PI);
          double checked_i_minus = bounded(i(time) - 1.0 / _i_, 1e-9, PI);

          CState state_plus(ref_a,
                            ref_e,
                            checked_i_plus,
                            ref_Omega,
                            ref_curlypi,
                            ref_L,
                            ELEMENTS,
                            time,
                            SUN,
                            AUD);

          dr_diff_plus = Distance_km(state_plus, BODY_STATE) / dr_max;
          dv_diff_plus = Velocity_kms(state_plus, BODY_STATE) / dv_max;

          CState state_minus(ref_a,
                             ref_e,
                             checked_i_minus,
                             ref_Omega,
                             ref_curlypi,
                             ref_L,
                             ELEMENTS,
                             time,
                             SUN,
                             AUD);

          dr_diff_minus = Distance_km(state_minus, BODY_STATE) / dr_max;
          dv_diff_minus = Velocity_kms(state_minus, BODY_STATE) / dv_max;
          dr_diff_max = std::max(dr_diff_plus, dr_diff_minus);
          dv_diff_max = std::max(dv_diff_plus, dv_diff_minus);
          diff_max = std::max(dr_diff_max, dv_diff_max);
//AO          diff_min = min(dr_diff_max, dv_diff_max);
        } while (diff_max >= 1.0);

        factor_avg[2] += _i_;
        if (_i_ > factor_max[2]) factor_max[2] = _i_;
        if (_i_ < factor_min[2]) factor_min[2] = _i_;

      }                         // end for

      const double elem[3] = { a(MJD_J2000)*(stateData->unitType==AUD?1.0:M2AU),
                               e(MJD_J2000),
                               i(MJD_J2000) };

      for (int loop = 0; loop < 3; loop++) {
        factor_max[loop] = factor_max[loop] / std::max(0.01, elem[loop]);
        factor_min[loop] = factor_min[loop] / std::max(0.01, elem[loop]);
        factor_avg[loop] = factor_avg[loop] / std::max(0.01, elem[loop])
                           / number_of_steps;
      }

	    mess.str("");
      mess << "\n\nsensitivity factors for a: max: "
        << factor_max[0] << " min: " << factor_min[0]
        << " avg: " << factor_avg[0] << endl;

      mess << "sensitivity factors for e: max: "
        << factor_max[1] << " min: " << factor_min[1]
        << " avg: " << factor_avg[1] << endl;

      mess << "sensitivity factors for i: max: "
        << factor_max[2] << " min: " << factor_min[2]
        << " avg: " << factor_avg[2] << endl;

	    Output(mess.str());

      factors_are_calculated = true;

      return SensitivityFactor(n, dr_max, dv_max);
    }
  }                             // end CBody::SensitivityFactor


  // -------------------------------------------------------------------------
  void CBody::Print(const double mjd, ofstream & out)
  {
    if (stateData->m_mjd != mjd)
      SetState(mjd);
    CState::Print(out);
  } // void CBody::Print(const double mjd, ofstream & out)

  // ---------------------------------------------------------------------------
  void CBody::printVrmlTrajectory(std::ostream& os, const int width
   , const int prec, const double time, const int n_steps, const double red
   , const double green, const double blue)
  {
    CBody::SetState(time);
    const double _a = CBody::a(time);
    const double _e = CBody::e(time);
    const double _i = CBody::i(time);
    const double _O = CBody::Omega(time);
    const double _o = CBody::omega(time);
    int i;

    if (IsElliptic())
    {
      os << "DEF " << getName() <<"Trajectory Shape\n";
      os << "{\n";
      os << " appearance Appearance\n";
      os << " {\n";
      os << "  material Material\n";
      os << "  {\n";
      os << "   emissiveColor " << red << " " << green << " " << blue << endl;
      os << "  } # material\n";
      os << " } # appearance\n";
      os << " geometry IndexedLineSet\n";
      os << " {\n";
      os << "  coord Coordinate\n";
      os << "  {\n";
      os << "   point\n";
      os << "   [";
      double period = OrbitPeriod(time);
      const double t_inc = (period*1.01)/double(n_steps>0?n_steps:1);
      double t = time;
      if (getId()==PLUTO) t = 1000.0;
      os << scientific << setprecision(prec) << showpos;
      for (i=0; i<(n_steps+1); ++i,t+=t_inc)
      {
        if ( (i%2) == 0 ) os << "\n    ";
        else os << "  ";
        os << setw(width) << x(t) << " ";
        os << setw(width) << y(t) << " ";
        os << setw(width) << z(t) << " ";
      } // for (i=0; i<(n_steps+1); ++i,t+=t_inc)
      os << "\n   ] # point\n";
      os << "  } # coord\n";
      os << "  coordIndex\n";
      os << "  [";
      os << noshowpos;
      const int w = int(log10(float(n_steps)))+1;
      for (i=0;i<(n_steps);++i)
      {
        if (i%10==0) os << "\n   ";
        os << setw(w) << i << ",";          
      } // for (i=0;i<steps;++i)
      os << n_steps << endl;
      os << "  ] # coordIndex\n";
      os << " } # geometry\n";
      os << "}\n\n";
    } // if (IsElliptic())
    else if (IsParabolic())
    {
    } // if (isParabolic())
    else if (IsHyperbolic())
    {
    } // if (isHyperbolic())
    else 
      throw std::runtime_error("CBody::printVrmlTrajectory(...): invalid orbit type!\n");
  } // void printVrmlTrajectory(...)

} // of namespace SpaceLib

// end of file

