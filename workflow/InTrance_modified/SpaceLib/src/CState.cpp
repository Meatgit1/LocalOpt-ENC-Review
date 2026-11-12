/*!
\file CState.cpp
\version 1.7
\brief A class for an astrodynamical state (CState).
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Tanja D. Schmidt <tanja.schmidt@dlr.de> (TS)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 25.09.2006 (AO)

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

//SpaceLib
#include <configSpaceLib.h>

//STL
#include <iomanip>

// C/C++ standard math functions
#include <math.h>

//SpaceLib
#include <MLConst.h>
#include <MLGlobal.h>
#include <CStateData.h>     // structure to hold the state's data in
#include <MLMatrix.h>
#include <Kepler.h>         // support functions to calculate Keplerian elements
#include <CState.h>         // declaration of CState - class


using namespace std;
using namespace Global;
using namespace MathLib;


namespace SpaceLib 
{

  // global temp. variables
  Matrix m_Rx(3,3), m_Ry(3,3), m_Rz(3,3), _PQW(3,3), T(3,3);
  Vector _r(3), _v(3), _h(3);

  // --------------------------------------------------------------------------
  CAttitude& CState::attitude(void) {return stateData->attitude;};

  // --------------------------------------------------------------------------
  void CState::setRefSolSys(CSolarSystem*solsys) { stateData->refSys=solsys; }


  // --------------------------------------------------------------------------
  CSolarSystem*CState::getRefSolSys(void) const { return stateData->refSys; }


  // --------------------------------------------------------------------------
  void CState::Set(
   const double v0,
   const double v1,
   const double v2,
   const double v3,
   const double v4,
   const double v5,
   const CoordType st,   /* = cartesian */
   const double a_mjd,   /* = MJD_J2000 */
   const BodyID cb, /* = SUN */
   const UnitType ut) throw (std::runtime_error)   /* = AUD */
  {

         stateData->cartesian.position=0.0;
         stateData->cartesian.velocity=0.0;

         stateData->sphere.position=0.0;
         stateData->sphere.velocity=0.0;

         stateData->m_a=stateData->m_e=stateData->m_i=0.0;
         stateData->m_Omega=stateData->m_omega=stateData->m_curlypi=0.0;
         stateData->m_M=stateData->m_E=stateData->m_f=stateData->m_L=0.0;

         stateData->m_mjd=a_mjd;
         stateData->m_central_body=cb;
         stateData->unitType=ut;

         if (!(v0==0.0&&v1==0.0&&v2==0.0&&v3==0.0&&v4==0.0&&v5==0.0))
         {

           switch (st) {

             case CARTESIAN:

               stateData->cartesian.position(0)=v0;
               stateData->cartesian.position(1)=v1;
               stateData->cartesian.position(2)=v2;

               stateData->cartesian.velocity(0)=v3;
               stateData->cartesian.velocity(1)=v4;
               stateData->cartesian.velocity(2)=v5;

               Cart2Polar();

               // the central body must have a mass to
               // calculate valid Keplerian elements
               if (cb!=SS_BARY) Cart2Kepler();
               break;

             case POLAR:
               
               stateData->sphere.position(0)=v0; // radius
               if (!bounded(stateData->sphere.position(0),-10.0,10.0))
               {
                 flush(cout<<" v0(2):"<<v0);
                 FatalError("CState::Set() => negative radius");
               }

               stateData->sphere.position(1)=v1; // angle phi (longitude)
               AssertAngularRange(stateData->sphere.position(1));

               stateData->sphere.position(2)=v2; // angle theta (azimuth)
               if (AssertElevationRange(stateData->sphere.position(2))==false)
               {
                 cout << "v0:" << v0 << endl;
                 cout << "v1:" << v1 << endl;
                 cout << "v2:" << v2 << endl;
                 cout << "v3:" << v3 << endl;
                 cout << "v4:" << v4 << endl;
                 cout << "v5:" << v5 << endl;
                 throw std::runtime_error("exc:Elevation out of Range!\n");
                 FatalError("Elevation out of Range!\n");
               }

               stateData->sphere.velocity(0)=v3; // radial velocity
               stateData->sphere.velocity(1)=v4; // angular velocity dphi
               stateData->sphere.velocity(2)=v5; // angular velocity dtheta

               Polar2Cart();

               // the central body must have a mass to
               // calculate valid Keplerian elements
               if (cb!=SS_BARY) Cart2Kepler();
               break;

             case KEPLERM:

               stateData->m_a=v0;
               if (!(stateData->m_a!=0.0))
               FatalError("CStata::Set() => a = 0.0");

               stateData->m_e=v1;
               AssertEccentricityRange(stateData->m_e);

               if (stateData->m_e==1.0) stateData->m_e-=1e-12;
               // just a trick to prevent exact parabolicity

               stateData->m_i=v2;
               AssertInclinationRange(stateData->m_i);

               if (stateData->m_i!=0.0) stateData->m_Omega=v3; else
               stateData->m_Omega=0.0;
               AssertAngularRange2pi(stateData->m_Omega);

               if (stateData->m_e!=0.0) stateData->m_omega=v4; else
               stateData->m_omega=0.0;
               AssertAngularRange2pi(stateData->m_omega);

               stateData->m_M=v5;
               AssertAngularRange2pi(stateData->m_M);

               stateData->m_curlypi=stateData->m_Omega+stateData->m_omega;
               AssertAngularRange2pi(stateData->m_curlypi);

               stateData->m_L=stateData->m_M+stateData->m_curlypi;
               AssertAngularRange2pi(stateData->m_L);

               CalcAnomaliesM(stateData->m_a, stateData->m_M, stateData->m_e,
                    stateData->m_E, stateData->m_f);
               Kepler2Cart();
               Cart2Polar();
               break;

             case KEPLERF:
               stateData->m_a=v0;
               if (!(stateData->m_a!=0.0))
               FatalError("CState::Set() => a = 0.0");

               stateData->m_e=v1;
               AssertEccentricityRange(stateData->m_e);
               if (stateData->m_e==1.0) stateData->m_e-=1e-12;
               // just a trick to prevent exact parabolicity

               stateData->m_i=v2;
               AssertInclinationRange(stateData->m_i);

               if (stateData->m_i!=0.0) stateData->m_Omega=v3; else
               stateData->m_Omega=0.0;
               AssertAngularRange2pi(stateData->m_Omega);

               if (stateData->m_e!=0.0) stateData->m_omega=v4; else
               stateData->m_omega=0.0;
               AssertAngularRange2pi(stateData->m_omega);

               stateData->m_f=v5;
               AssertAngularRange2pi(stateData->m_f);
               stateData->m_curlypi=stateData->m_Omega+stateData->m_omega;
               AssertAngularRange2pi(stateData->m_curlypi);

               CalcAnomaliesf(stateData->m_a, stateData->m_f, stateData->m_e,
                    stateData->m_E, stateData->m_M);

               stateData->m_L=stateData->m_M+stateData->m_curlypi;
               AssertAngularRange2pi(stateData->m_L);

               Kepler2Cart();
               Cart2Polar();
               break;

             case ELEMENTS:
               stateData->m_a=v0;
               if (!(stateData->m_a!=0.0))
               FatalError("CState::Set() => a = 0");

               stateData->m_e=v1;
               AssertEccentricityRange(stateData->m_e);
               if (stateData->m_e==1.0) stateData->m_e-=1e-12;
               // just a trick to prevent exact parabolicity

               stateData->m_i=v2;
               AssertInclinationRange(stateData->m_i);

               if (stateData->m_i!=0.0) stateData->m_Omega=v3; else
               stateData->m_Omega=0.0;
               AssertAngularRange2pi(stateData->m_Omega);

               if (stateData->m_e!=0.0) stateData->m_curlypi=v4; else
               stateData->m_curlypi=stateData->m_Omega;
               AssertAngularRange2pi(stateData->m_curlypi);

               stateData->m_L=v5;
               AssertAngularRange2pi(stateData->m_L);

               stateData->m_omega=stateData->m_curlypi-stateData->m_Omega;
               AssertAngularRange2pi(stateData->m_omega);

               stateData->m_M=stateData->m_L-stateData->m_curlypi;
               AssertAngularRange2pi(stateData->m_M);

               CalcAnomaliesM(stateData->m_a, stateData->m_M, stateData->m_e,
                    stateData->m_E, stateData->m_f);

               Kepler2Cart();
               Cart2Polar();
           } // switch(st)
         } //if (!(V0==0.0 ...
  }; // end CState::Set


  // ---------------------------------------------------------------------------
  void CState::Cart2Polar(void) 
  {

    stateData->sphere.position=Norm(stateData->cartesian.position);

    const double rho=
         sqrt(stateData->cartesian.position(0)*
         stateData->cartesian.position(0)+stateData->cartesian.position(1)*
         stateData->cartesian.position(1));

    if (rho!=0.0)
      stateData->sphere.position(2)=atan2(stateData->cartesian.position(2),rho);
    else 
    {
      if (stateData->cartesian.position(2)>0.0)
        stateData->sphere.position(2)=HALFPI;
      else if (stateData->cartesian.position(2)<0.0)
        stateData->sphere.position(2)=-HALFPI; else
        stateData->sphere.position(2)=0.0;
    }

    const double p=
         2.0*atan2(stateData->cartesian.position(1),
         fabs(stateData->cartesian.position(0))+rho);

    if (stateData->cartesian.position(0)==0.0&&
         stateData->cartesian.position(1)==0.0)
           stateData->sphere.position(1)=0.0;
    else if (stateData->cartesian.position(0)>=0.0&&
         stateData->cartesian.position(1)>=0.0)
           stateData->sphere.position(1)=p;
    else if (stateData->cartesian.position(0)>=0.0&&
         stateData->cartesian.position(1)<0.0)
           stateData->sphere.position(1)=TWOPI+p; else
      stateData->sphere.position(1)=PI-p;

    R_y(-stateData->sphere.position(2), m_Ry);
    R_z(stateData->sphere.position(1), m_Rz);
    MultiplyMatrices(m_Ry, m_Rz, T);
    MultiplyMatrixVector(T, stateData->cartesian.velocity, _v); 

    // dr
    stateData->sphere.velocity(0) = _v(0);
    // dphi
    stateData->sphere.velocity(1) =
      _v(1)/(stateData->sphere.position(0)*cos(stateData->sphere.position(2)));
    //dtheta
    stateData->sphere.velocity(2) = _v(2)/stateData->sphere.position(0);

    AssertAngularRange(stateData->sphere.position(1)); // angle phi
  };


  // ---------------------------------------------------------------------------
  void CState::Cart2Kepler(void)
  {

    
    //const Vector r  =stateData->cartesian.position;
    //const Vector v  =stateData->cartesian.velocity;
    //const Vector h  =Cross(r, v); // areal velocity
    Cross(stateData->cartesian.position, stateData->cartesian.velocity, _h);
    //const double _r_=Norm(r); // normalized radius vector
    //const double _h_=Norm(h); // normalized areal velocity vector
    //const double v2 =Dot(v, v); // vector product v*v
    //const double h2 =Dot(h, h); // vector product h*h
    const double _r_=Norm(stateData->cartesian.position);
    const double _h_=Norm(_h); // normalized areal velocity vector
    const double v2 =Dot(stateData->cartesian.velocity, stateData->cartesian.velocity); // vector product v*v
    const double h2 =Dot(_h, _h); // vector product h*h

    if (h2==0.0)
      FatalError("CState::Cart2Kepler(): state corresponds to a straight line");

    // inclination
    stateData->m_i=atan2(sqrt(_h(0)*_h(0)+_h(1)*_h(1)), _h(2));

    // calculation of longitude of ascending node
    // and argument of latitude (longitude?)

    double u;

    if (stateData->m_i==0.0)
    { // non inclined orbit
      stateData->m_Omega=0.0;
//      u=atan2(r(1), r(0));
      u=atan2(stateData->cartesian.position(1),stateData->cartesian.position(0));
    } else { // inclined orbit
//      stateData->m_Omega=atan2(h(0), -h(1));
//      u=atan2(r(2)*_h_, -r(0)*h(1)+r(1)*h(0));
      stateData->m_Omega=atan2(_h(0), -_h(1));
      u=atan2(stateData->cartesian.position(2)*_h_, -stateData->cartesian.position(0)*_h(1)+stateData->cartesian.position(1)*_h(0));
    }

    // calculation of semilatus rectum

    double p=0.0;
    double lGM=0.0;

    if (stateData->m_central_state!=0) 
    {
      if (stateData->unitType==AUD) 
      {
        lGM=stateData->m_central_state->getGM(AUD);
      } else 
      { 
        lGM=stateData->m_central_state->getGM(SI); 
      }
    } else
    {
      lGM=GM_old();
    }

    if (lGM!=0.0) {

      p=h2/lGM;
      stateData->m_a=1.0/(2.0/_r_-v2/lGM); // semi-major axis
      stateData->m_e=0.0; // init eccentricity

      // eccentricity for elliptic and hyperbolic orbits
      if ((1.0-p/stateData->m_a)>1e-15)
        stateData->m_e=sqrt(1.0-p/stateData->m_a);

      if (stateData->m_e==0.0) { // if orbit is circular,
        stateData->m_M=stateData->m_E=stateData->m_f=u; // set different anomalies
        stateData->m_omega=0.0;
      } else if (stateData->m_e<=1.0) { // orbit is elliptic
        if (stateData->m_e==1.0) stateData->m_e-=1e-12;
        // just a trick to prevent exact parabolicity

        // calculation of mean anomaly
        const double eCosE=1.0-_r_/stateData->m_a; // e*cos(E)
//        const double eSinE=Dot(r, v)/sqrt(lGM*stateData->m_a);
        const double eSinE=Dot(stateData->cartesian.position, stateData->cartesian.velocity)/sqrt(lGM*stateData->m_a);

        // eccentric anomaly
        stateData->m_E=atan2(eSinE, eCosE);

        // mean anomaly
        stateData->m_M=stateData->m_E-eSinE;
        AssertAngularRange2pi(stateData->m_M);

        // calculation of argument of perihelion

        // true anomaly
        stateData->m_f=
             atan2(sqrt(1.0-stateData->m_e*stateData->m_e)*eSinE,
             eCosE-stateData->m_e*stateData->m_e);

        // argument of perihelion
        stateData->m_omega=u-stateData->m_f;
        AssertAngularRange2pi(stateData->m_omega);

      } else { // orbit must be hyperbolic

        // calculation of mean anomaly

        const double k=sqrt(lGM*fabs(stateData->m_a));

        // e*cosH(H)
        const double eCoshH=1.0+_r_/fabs(stateData->m_a);

        // e*sinH(H)
        //const double eSinhH=Dot(r, v)/k;
        const double eSinhH=Dot(stateData->cartesian.position, stateData->cartesian.velocity)/k;
        
        const double TanhH=eSinhH/eCoshH;

        // hyperbolic anomaly (= artanh x)
        if (fabs(TanhH)<1.0)
          stateData->m_E=0.5*log((1.0+TanhH)/(1.0-TanhH)); else
          FatalError("Hyperb(): hyperbolic anomaly out of range");

        stateData->m_M=eSinhH-stateData->m_E; // mean hyperbolic anomaly

        // calculation of argument of perihelion

        // true anomaly
        stateData->m_f=
             2.0*atan(sqrt((stateData->m_e+1.0)/(stateData->m_e-1.0))*
             tanh(0.5*stateData->m_E));

        // argument of perihelion
        stateData->m_omega=u-stateData->m_f;
        AssertAngularRange2pi(stateData->m_omega);
      }

      stateData->m_curlypi=stateData->m_Omega+stateData->m_omega;
      AssertAngularRange2pi(stateData->m_curlypi);

      stateData->m_L=stateData->m_M+stateData->m_curlypi;
      AssertAngularRange2pi(stateData->m_L);
    }
  } // end CState::Cart2Kepler


  // calculates from sphere polar coordinates/velocities to
  // cartesian position/velocity coordinates
  // --------------------------------------------------------------------------
  void CState::Polar2Cart(void)
  {

    // x = r*cos(phi)*cos(theta)
    stateData->cartesian.position(0)=
         stateData->sphere.position(0)*cos(stateData->sphere.position(1))*
         cos(stateData->sphere.position(2));

    // y = r*sin(phi)*cos(theta)
    stateData->cartesian.position(1)=
         stateData->sphere.position(0)*sin(stateData->sphere.position(1))*
         cos(stateData->sphere.position(2));

    // z = r*sin(theta)
    stateData->cartesian.position(2)=
         stateData->sphere.position(0)*sin(stateData->sphere.position(2));

    _v(0) = stateData->sphere.velocity(0);
    _v(2) = stateData->sphere.position(0)*stateData->sphere.velocity(2);
    _v(1) = stateData->sphere.position(0)*stateData->sphere.velocity(1)*
         cos(stateData->sphere.position(2));

    R_z(-stateData->sphere.position(1), m_Rz);
    R_y(stateData->sphere.position(2), m_Ry);
    MultiplyMatrices(m_Rz, m_Ry, T);
    MultiplyMatrixVector(T, _v, stateData->cartesian.velocity);
  } // void CState::Polar2Cart(void)


  // --------------------------------------------------------------------------
  void CState::Kepler2Cart(void) 
  {
    double lGM=0.0;

    if (stateData->m_central_state!=0) 
    {
      if (stateData->unitType==AUD) 
      {
        lGM=stateData->m_central_state->getGM(AUD);
      } else 
      { 
        lGM=stateData->m_central_state->getGM(SI); 
      }
    } else lGM=GM_old();

    if (IsParabolic()) 
    { // parabolic state
      const int MAXIT = 25;
      const double p  =stateData->m_a*(1.0-stateData->m_e*stateData->m_e);
      const double q  =p/(1.0+stateData->m_e);
      double fac      =0.5*stateData->m_e;
      const double k  =sqrt(lGM/(q*(1.0+stateData->m_e)));
      const double tau=
           stateData->m_M*sqrt(stateData->m_a*stateData->m_a*stateData->m_a);
      double E20, E2=0.0;
      double A, B, u, u2, c1, c2, c3;

      int loop=0;

      do 
      {
        E20=E2;

        A=1.5*sqrt(fac/(q*q*q))*tau;
        B=pow(sqrt(A*A+1.0)+A, 1.0/3.0);
        u=B-1.0/B;
        u2=u*u;
        E2=u2*(1.0-stateData->m_e)/fac;

        Stumpff(E2, c1, c2, c3);

        fac=3.0*stateData->m_e*c3;

        if (++loop == MAXIT) 
        {
          cerr<<" Konvergenzprobleme in Parab" <<endl;
          break;
        }
      } while (fabs(E2-E20)>=EPS);

      const double R=q*(1.0+u2*c2*stateData->m_e/fac);

      _r(0) = q*(1.0-u2*c2/fac);
      _r(1) = q*sqrt((1.0+stateData->m_e)/fac)*u*c1;
      _r(2) = 0.0;

      _v(0) = -k*_r(1)/R;
      _v(1) = k*(_r(0)/R+stateData->m_e);
      _v(2) = 0.0;

    } else if (IsElliptic()) 
    {
      // elliptic or circular state
      // calculation of eccentric anomaly

      stateData->m_E=stateData->m_M; // for circular orbits
      if (stateData->m_e!=0.0) 
        stateData->m_E=EccentricAnomaly(stateData->m_M, stateData->m_e);

      const double cosE=cos(stateData->m_E);
      const double sinE=sin(stateData->m_E);

      // perifocal coordinates
      const double fac=sqrt(1.0 - stateData->m_e*stateData->m_e);
      const double _r_=stateData->m_a*(1.0-stateData->m_e*cosE); // distance
      const double _v_=sqrt(lGM*stateData->m_a)/_r_; // velocity

      _r(0) = stateData->m_a*(cosE-stateData->m_e);
      _r(1) = stateData->m_a*fac*sinE;
      _r(2) = 0.0;

      _v(0) = -_v_*sinE;
      _v(1) = _v_*fac*cosE;
      _v(2) = 0.0;

    } else if (IsHyperbolic())
    {
      // hyperbolic state
      const double _a_=fabs(stateData->m_a);

      // calculation of hyperbolic anomaly
      const double H=HyperbolicAnomaly(stateData->m_M, stateData->m_e);

      const double coshH=cosh(H);
      const double sinhH=sinh(H);

      // perifocal coordinates
      const double k  =sqrt(lGM/_a_);
      const double fac=sqrt(stateData->m_e*stateData->m_e-1.0);
      const double rho=stateData->m_e*coshH-1.0;

      _r(0) = _a_*(stateData->m_e-coshH);
      _r(1) = _a_*fac*sinhH;
      _r(2) = 0.0;

      _v(0) = -k*sinhH/rho;
      _v(1) = k*fac*coshH/rho;
      _v(2) = 0.0;
    } else 
    {
      FatalError("Kepler2Cart(): state must be elliptic, parabolic \
      or hyperbolic");
    }

    R_z(-stateData->m_omega, m_Rz);
    R_x(-stateData->m_i, m_Rx);
    MultiplyMatrices(m_Rx, m_Rz, m_Ry);
    R_z(-stateData->m_Omega, m_Rz);
    MultiplyMatrices(m_Rz, m_Ry, _PQW);

    stateData->cartesian.position = _PQW*_r;
    stateData->cartesian.velocity = _PQW*_v;
  }; // end CState::Kepler2Cart


  // --------------------------------------------------------------------------
  double CState::GM_AUd(void)const {

    double lGM=0.0;

    if (stateData->m_central_state!=0) {
      lGM=stateData->m_central_state->getGM(AUD);
    } else {
      switch (stateData->m_central_body) {

        // the values below are taken from
        // http://ssd.jpl.nasa.gov/astro_constants.html

        case SUN:
          lGM=0.0002959122082855911;
          break;
        case MERCURY:
          lGM=4.912547451450812e-11;
          break;
        case VENUS:
          lGM=7.243452486162703e-10;
          break;
        case EARTH:
          lGM=8.88769239011351e-10;
          break;
        case MARS:
          lGM=9.549535105779258e-11;
          break;
        case JUPITER:
          lGM=2.825345909524226e-07;
          break;
        case SATURN:
          lGM=8.459715185680659e-08;
          break;
        case URANUS:
          lGM=1.292024916781969e-08;
          break;
        case NEPTUNE:
          lGM=1.524358900784276e-08;
          break;
        case PLUTO:
          lGM=2.18869976542597e-12;
          break;
        case LUNA:
          lGM=1.093189565989891e-11;
          break;
        default:
          lGM=0.0;
      }
    }

    return lGM;
  };


  // --------------------------------------------------------------------------
  double CState::GM_SI(void)const { return GM_AUd()*GM_AUD2GM_SI; };


  // --------------------------------------------------------------------------
  double CState::GM_old(void)const {
    double gm=0.0;
    switch (stateData->unitType) {
      case AUD:
        gm=GM_AUd();
        break;
      case SI:
        gm=GM_SI();
        break;
      default:
        FatalError("double CState::GM_old(void) => unknown unit type!");
    }
    return gm;
  }


  // --------------------------------------------------------------------------
  double CState::E_PARABOLIC1(void)const {
    double limit;

    if (stateData->unitType==AUD)
      limit=1e6; else if (stateData->unitType==SI)
      limit=1e17; else
      FatalError("CState::E_PARABOLIC1(): unknown unit type");

    return limit;
  };


  // --------------------------------------------------------------------------
  double CState::E_PARABOLIC2(void)const {
    double limit;
    if (stateData->unitType==AUD) limit=1e9; else if (stateData->unitType==SI)
      limit=1e20; else
      FatalError("CState::E_PARABOLIC2(): unknown unit type");
    return limit;
  };


  // --------------------------------------------------------------------------
  double CState::H_PARABOLIC1(void)const {
    double limit;
    if (stateData->unitType==AUD) limit=-1e6; else if (stateData->unitType==SI)
      limit=-1e17; else
      FatalError("CState::H_PARABOLIC1(): unknown unit type");
    return limit;
  };



  // --------------------------------------------------------------------------
  double CState::H_PARABOLIC2(void)const {
    double limit;
    if (stateData->unitType==AUD) limit=-1e9; else if (stateData->unitType==SI)
      limit=-1e20; else
      FatalError("CState::H_PARABOLIC2(): unknown unit type");
    return limit;
  };

  // --------------------------------------------------------------------------
  CState::~CState(void)
  {
    SafeDelete(stateData);
  }

  // --------------------------------------------------------------------------
  CState::CState(
   const double a_mjd, // = MJD_J2000
   const BodyID cb, // = SUN
   const UnitType ut) // = AUD
  : stateData(0)
  {
    stateData=new CStateData;
    if (stateData==0)
      FatalError("CState:CState error with memory allocation");

    //stateData->useDE405=false;
    //stateData->indexDE405=DE405_NO_BODY;
    //stateData->type=BODY;
    //stateData->refSys = 0;
    //stateData->GM_AUD=0;
    //stateData->GM_SI=0;
    //stateData->m_central_state=0;
    stateData->m_central_body=cb;
    stateData->unitType=ut;
    stateData->m_mjd=a_mjd;

    Set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, CARTESIAN, a_mjd, cb, ut);
  };


  // --------------------------------------------------------------------------
  CState::CState(const CState& to_copy)
    :stateData(new CStateData(*to_copy.stateData))
  {
    if (stateData==0)
      FatalError("CState:CState error with memory allocation");
  }


  // --------------------------------------------------------------------------
  CState::CState(
  const Vector&V,       // position/velocity wrt center body
  const CoordType st,   // = CARTESIAN
  const double a_mjd,   // = MJD_J2000
  const BodyID cb, // = SUN
  const UnitType ut):   // = AUD
   stateData(new CStateData)
   {

    if (V.getDimension()!=6)
      FatalError("CState::CState(): wrong vector size");

    if (stateData==0)
      FatalError("CState:CState error with memory allocation");

    //stateData->useDE405=false;
    //stateData->indexDE405=DE405_NO_BODY;
    //stateData->type=BODY;
    //stateData->refSys = 0;
    //stateData->GM_AUD=0;
    //stateData->GM_SI=0;
    //stateData->m_central_state=0;
    stateData->m_central_body=cb;
    stateData->unitType=ut;
    stateData->m_mjd=a_mjd;

    Set(V(0), V(1), V(2), V(3), V(4), V(5), st, a_mjd, cb, ut);
  };


  // --------------------------------------------------------------------------
  CState::CState(
   const double*p,
   const CoordType st, // = CARTESIAN
   const double a_mjd, // = MJD_J2000
   const BodyID cb, // = SUN
   const UnitType ut): // = AUD
    stateData(new CStateData)
  {
    if (stateData==0)
      FatalError("CState:CState error with memory allocation");

    //stateData->useDE405=false;
    //stateData->indexDE405=DE405_NO_BODY;
    //stateData->type=BODY;
    //stateData->refSys = 0;
    //stateData->GM_AUD=0;
    //stateData->GM_SI=0;
    //stateData->m_central_state=0;
    stateData->m_central_body=cb;
    stateData->unitType=ut;
    stateData->m_mjd=a_mjd;

    if (p!=0)
      Set(p[0], p[1], p[2], p[3], p[4], p[5], st, a_mjd, cb, ut);
  };


  // --------------------------------------------------------------------------
    CState::CState(const double v0, const double v1, const double v2,
      const double v3, const double v4, const double v5, const CoordType st, // = CARTESIAN
      const double a_mjd, // = MJD_J2000
      const BodyID cb, // = SUN
      const UnitType ut): // = AUD
    stateData(new CStateData)
    {
      if (stateData==0)
        FatalError("CState:CState error with memory allocation");
      //stateData->useDE405=false;
      //stateData->indexDE405=DE405_NO_BODY;
      //stateData->type=BODY;
      //stateData->refSys = 0;
      //stateData->GM_AUD=0;
      //stateData->GM_SI=0;
      //stateData->m_central_state=0;
      stateData->m_central_body=cb;
      stateData->unitType=ut;
      stateData->m_mjd=a_mjd;

      Set(v0, v1, v2, v3, v4, v5, st, a_mjd, cb, ut);
    };



  // --------------------------------------------------------------------------
  const CState& CState::operator=(const CState& toCopy) 
  {
    if ((&toCopy) != this) 
    {
      if (stateData==0)
        FatalError("CState::operator= data pointer not allocated!");
      *stateData=*toCopy.stateData;
    }
    return (*this);
  }


  // --------------------------------------------------------------------------
  bool CState::operator==(const CState& rhs) 
  {
    return (*stateData) == (*(rhs.stateData));
  }


  // --------------------------------------------------------------------------
  const CAttitude& CState::getAttitude(void) const
  {
    return stateData->attitude;
  }


  // --------------------------------------------------------------------------
  const MathLib::Vector& CState::getAngularVelocity(void) const
  {
    return stateData->angular_velocity;
  }


  // --------------------------------------------------------------------------
  void CState::SetState(
   const Vector& V,      // Vector of data to set
   const CoordType st,   // = CARTESIAN
   const double a_mjd,   // = MJD_J2000
   const BodyID cb, // = SUN
   const UnitType ut,    // = AUD
   const CStateType aST,  // = STATETYPE_BODY
   const bool useDE405,  // = false
   const DE405StateIndex indexDE405, // = 0
   CState* aCS           // center state
  ) {

    if (V.getDimension() != 6) FatalError("CState::SetState(): wrong vector size");

    stateData->type=aST;
    stateData->useDE405=useDE405;
    stateData->indexDE405=indexDE405;
    if (aCS != 0) stateData->m_central_state = aCS;

    Set(V(0), V(1), V(2), V(3), V(4), V(5), st, a_mjd, cb, ut);

  };


  // --------------------------------------------------------------------------
  void CState::SetState(const double*p, const CoordType st, /* = CARTESIAN */
  const double a_mjd, /* = MJD_J2000 */
  const BodyID cb, /* = SUN */
  const UnitType ut) {
    /* = AUD */

    Set(p[0], p[1], p[2], p[3], p[4], p[5], st, a_mjd, cb, ut);
  };


  // --------------------------------------------------------------------------
  void CState::SetState(const double v0, const double v1, const double v2,
       const double v3, const double v4, const double v5, const CoordType st,
       /* = CARTESIAN */
       const double a_mjd, /* = MJD_J2000 */
       const BodyID cb, /* = SUN */
       const UnitType ut) /* = AUD */
  { Set(v0, v1, v2, v3, v4, v5, st, a_mjd, cb, ut); };


  // ---------------------------------------------------------------------------
  void CState::useDE405(bool _val) { stateData->useDE405=_val; }


  // --------------------------------------------------------------------------
  bool CState::useDE405(void) { return stateData->useDE405; }


  // --------------------------------------------------------------------------
  DE405StateIndex CState::getIndexDE405(void) { return stateData->indexDE405; }


  // --------------------------------------------------------------------------
  void CState::setType(CStateType aType) { stateData->type=aType; }


  // --------------------------------------------------------------------------
  CStateType CState::getType(void)const { return stateData->type; }


  // --------------------------------------------------------------------------
  void CState::setCenterState(CState*sCS) { stateData->m_central_state=sCS; }


  // --------------------------------------------------------------------------
  CState*CState::getCenterState(void) { return stateData->m_central_state; }


  // --------------------------------------------------------------------------
  void CState::setCenterBodyId(const BodyID aPI) 
  { 
    stateData->m_central_body=aPI; 
  }


  // --------------------------------------------------------------------------
  BodyID CState::getCenterBodyId(void)const {
    return stateData->m_central_body;
  };


  // --------------------------------------------------------------------------
  UnitType CState::Unit(void)const { return stateData->unitType; }


  // --------------------------------------------------------------------------
  double CState::MJD(void)const { return stateData->m_mjd; }


  // --------------------------------------------------------------------------
  void CState::setGM(double aGM, UnitType aType) {
    switch (aType) {
      case SI:
        stateData->GM_AUD=aGM*GM_SI2GM_AUD;
        stateData->GM_SI=aGM;
        break;
      case AUD:
        stateData->GM_AUD=aGM;
        stateData->GM_SI=aGM*GM_AUD2GM_SI;
        break;
      default:
        FatalError("CState::setGM(double, UnitType) => Unknown unit type!");
    }
  }


  // --------------------------------------------------------------------------
  double CState::getGM(UnitType aType) const {
    double gm=0.0;
    switch (stateData->unitType) {
      case SI:
        gm=(aType==SI?stateData->GM_SI:stateData->GM_AUD);
        break;
      case AUD:
        gm=(aType==AUD?stateData->GM_AUD:stateData->GM_SI);
        break;
      default:
        FatalError("double CState::getGM(UnitType) => Unknown unit type!");
    }
    return gm;
  }


  double CState::x(void)const { return stateData->cartesian.position(0); }


  //-----------------------------------------------------------------------------
  double CState::dx(void)const { return stateData->cartesian.velocity(0); }


  //-----------------------------------------------------------------------------
  double CState::y(void)const { return stateData->cartesian.position(1); }


  //-----------------------------------------------------------------------------
  double CState::dy(void)const { return stateData->cartesian.velocity(1); }


  //-----------------------------------------------------------------------------
  double CState::z(void)const { return stateData->cartesian.position(2); }


  //-----------------------------------------------------------------------------
  double CState::dz(void)const { return stateData->cartesian.velocity(2); }


  // ----------------------------------------------------------------------------
  double CState::r(void)const { return stateData->sphere.position(0); }


  // ----------------------------------------------------------------------------
  double CState::phi(void)const { return stateData->sphere.position(1); }


  // ----------------------------------------------------------------------------
  double CState::theta(void)const { return stateData->sphere.position(2); }


  // ----------------------------------------------------------------------------
  double CState::dr(void)const { return stateData->sphere.velocity(0); }


  // ----------------------------------------------------------------------------
  double CState::dphi(void)const { return stateData->sphere.velocity(1); }


  // ----------------------------------------------------------------------------
  double CState::dtheta(void)const { return stateData->sphere.velocity(2); }


  // ----------------------------------------------------------------------------
  double CState::a(void)const { return stateData->m_a; }


  // ----------------------------------------------------------------------------
  double CState::e(void)const { return stateData->m_e; }


  // ----------------------------------------------------------------------------
  double CState::i(void)const { return stateData->m_i; }


  // ----------------------------------------------------------------------------
  double CState::Omega(void)const { return stateData->m_Omega; }


  // ----------------------------------------------------------------------------
  double CState::omega(void)const { return stateData->m_omega; }


  // ----------------------------------------------------------------------------
  double CState::curlypi(void)const { return stateData->m_curlypi; }


  // ----------------------------------------------------------------------------
  double CState::M(void)const { return stateData->m_M; }


  // ----------------------------------------------------------------------------
  double CState::E(void)const { return stateData->m_E; }


  // ----------------------------------------------------------------------------
  double CState::f(void)const { return stateData->m_f; }


  // ----------------------------------------------------------------------------
  double CState::L(void)const { return stateData->m_L; }


  // ----------------------------------------------------------------------------
  double CState::vr(void)const {
    return stateData->sphere.velocity(0); // same as CState::dr()
  }


  // ----------------------------------------------------------------------------
  double CState::vphi(void)const {
    return stateData->sphere.velocity(1) // dphi
    *stateData->sphere.position(0) // r
    *cos(stateData->sphere.position(2)); // cos(theta)
  }


  // ----------------------------------------------------------------------------
  double CState::vtheta(void)const {
    return stateData->sphere.velocity(2) // dtheta
    *stateData->sphere.position(0); // r;
  }


  // ----------------------------------------------------------------------------
  double CState::h(void)const {
    double result=stateData->m_e*sin(stateData->m_Omega+stateData->m_omega);
    AssertAngularRange(result);
    return result;
  };


  // ----------------------------------------------------------------------------
  double CState::k(void)const {
    double result=stateData->m_e*cos(stateData->m_Omega+stateData->m_omega);
    AssertAngularRange(result);
    return result;
  };


  // ----------------------------------------------------------------------------
  double CState::p(void)const {
    double result=sin(0.5*stateData->m_i)*sin(stateData->m_Omega);
    AssertAngularRange(result);
    return result;
  };


  // ----------------------------------------------------------------------------
  double CState::q(void)const {
    double result=sin(0.5*stateData->m_i)*cos(stateData->m_Omega);
    AssertAngularRange(result);
    return result;
  };


  // ----------------------------------------------------------------------------
  double CState::l(void)const {
    double result=stateData->m_Omega+stateData->m_omega+stateData->m_M;
    AssertAngularRange(result);
    return result;
  };


  // ----------------------------------------------------------------------------
  double CState::Momentum(void)const {

    return Norm(Cross(stateData->cartesian.position,
         stateData->cartesian.velocity));
  };


  // ----------------------------------------------------------------------------
  double CState::SemiMinorAxis(void)const {
    if (IsElliptic())
      return stateData->m_a*sqrt(1.0-stateData->m_e*stateData->m_e);
    else if (IsHyperbolic())
      return stateData->m_a*sqrt(stateData->m_e*stateData->m_e-1.0); else
      return 0.0;
  };


  // ----------------------------------------------------------------------------
  double CState::PotentialEnergy(void)const {
    if (stateData->sphere.position(0)!=0.0)
      return -GM_old()/stateData->sphere.position(0); else
      return 0.0;
  }


  // ----------------------------------------------------------------------------
  double CState::Velocity(void)const {
    return Norm(stateData->cartesian.velocity);
  };


  // ----------------------------------------------------------------------------
  double CState::KineticEnergy(void)const {
    return 0.5*(stateData->cartesian.velocity(0)*
         stateData->cartesian.velocity(0)+stateData->cartesian.velocity(1)*
         stateData->cartesian.velocity(1)+stateData->cartesian.velocity(2)*
         stateData->cartesian.velocity(2));
  };


  // ----------------------------------------------------------------------------
  double CState::OrbitPeriod(void)const {
    double _GM=GM_old();
    double local_a=stateData->m_a<0.0?-1.0*stateData->m_a:stateData->m_a;

    if (_GM!=0.0) return 2.0*PI*sqrt(local_a*local_a*local_a/_GM); else
      return 0.0;
  }


  // ----------------------------------------------------------------------------
  double CState::MeanMotion(void)const {
    double _OP=OrbitPeriod();
    if (_OP==0.0)
      return 0.0; else
      return TWOPI/_OP;
  }


  // ----------------------------------------------------------------------------
  double CState::SemilatusRectum(void)const {
    double _GM=GM_old();
    if (stateData->type==BODY&&_GM!=0.0)
      return pow(Momentum(), 2)/_GM; else
      return 0.0;
  }


  // ----------------------------------------------------------------------------
  double CState::Zeta(void)const {
    return atan2(stateData->sphere.velocity(2),
         stateData->sphere.velocity(1)*cos(stateData->sphere.position(2)));
  }


  // ----------------------------------------------------------------------------
  Vector CState::Cart(const UnitType unit_type) const {

    switch (unit_type) {
      case SI:
        if (stateData->unitType==AUD)
          return Stack(stateData->cartesian.position*AU2M,
               stateData->cartesian.velocity*AUD2MS);
        break;
      case AUD:
        if (stateData->unitType==SI)
          return Stack(stateData->cartesian.position*M2AU,
               stateData->cartesian.velocity*MS2AUD);
        break;
      case STORED:;
    }
    // no conversion is needed
    return Stack(stateData->cartesian.position, stateData->cartesian.velocity);
  };


  // ----------------------------------------------------------------------------
  Vector CState::Polar(void)const {

    return Stack(stateData->sphere.position, stateData->sphere.velocity);
  };


  // ----------------------------------------------------------------------------
  Vector CState::Polarv(void)const {
    return Vector(stateData->sphere.position(0), stateData->sphere.position(1),
         stateData->sphere.position(2), vr(), vphi(), vtheta());
  };


  // ----------------------------------------------------------------------------
  Vector CState::KeplerM(void)const {
    return Vector(stateData->m_a, stateData->m_e, stateData->m_i,
         stateData->m_Omega, stateData->m_omega, stateData->m_M);
  };


  // ----------------------------------------------------------------------------
  Vector CState::Keplerf(void)const {
    return Vector(stateData->m_a, stateData->m_e, stateData->m_i,
         stateData->m_Omega, stateData->m_omega, stateData->m_f);
  };


  // ----------------------------------------------------------------------------
  Vector CState::rVector(void)const {

    return stateData->cartesian.position/Norm(stateData->cartesian.position);

  };


  // ----------------------------------------------------------------------------
  Vector CState::tVector(void)const {
    const Vector t=Cross(hVector(), rVector());
    const double _t_=Norm(t);
    return t/_t_;
  };


  // ----------------------------------------------------------------------------
  Vector CState::hVector(void)const {

    const Vector h=Cross(stateData->cartesian.position, stateData->cartesian.velocity); //Cross(r,v);

    const double _h_=Norm(h);
    return h/_h_;
  };


  // ----------------------------------------------------------------------------
  Vector CState::PVector(void)const {
    const double sini    =sin(stateData->m_i);
    const double cosi    =cos(stateData->m_i);
    const double sinOmega=sin(stateData->m_Omega);
    const double cosOmega=cos(stateData->m_Omega);
    const double sinomega=sin(stateData->m_omega);
    const double cosomega=cos(stateData->m_omega);

    return Vector(cosomega*cosOmega-sinomega*cosi*sinOmega,
         cosomega*sinOmega+sinomega*cosi*cosOmega, sinomega*sini);
  };


  // ----------------------------------------------------------------------------
  Vector CState::QVector(void)const {
    const double sini    =sin(stateData->m_i);
    const double cosi    =cos(stateData->m_i);
    const double sinOmega=sin(stateData->m_Omega);
    const double cosOmega=cos(stateData->m_Omega);
    const double sinomega=sin(stateData->m_omega);
    const double cosomega=cos(stateData->m_omega);

    return Vector(-sinomega*cosOmega-cosomega*cosi*sinOmega,
         -sinomega*sinOmega+cosomega*cosi*cosOmega, cosomega*sini);
  };


  // ----------------------------------------------------------------------------
  Vector CState::RVector(void)const {
    const double sini    =sin(stateData->m_i);
    const double cosi    =cos(stateData->m_i);
    const double sinOmega=sin(stateData->m_Omega);
    const double cosOmega=cos(stateData->m_Omega);

    return Vector(sini*sinOmega, -sini*cosOmega, cosi);
  };


  // ----------------------------------------------------------------------------
  double CState::GM(void)const {

    double result=0.0;

    if (Unit()==AUD)
      result=stateData->GM_AUD; else if (Unit()==SI)
      result=stateData->GM_AUD*GM_AUD2GM_SI; else
      FatalError("CState::GM(): unknown unit type");

    return result;
  };


  // ----------------------------------------------------------------------------
  bool CState::IsElliptic(void)const {
    const bool small_M=
         (fabs(stateData->m_M)<=E_M_SMALL||
         fabs(stateData->m_M-TWOPI)<=E_M_SMALL)?true:false;

    if (stateData->m_a<0.0)
      return false; // hyperbolic
    if (stateData->m_a>=E_PARABOLIC2()) // parabolic
      return false;
    if (stateData->m_a>=E_PARABOLIC1()&&small_M) // parabolic
      return false;

    return true;
  };


  // ----------------------------------------------------------------------------
  bool CState::IsParabolic(void)const {
    const bool small_ME=
         (fabs(stateData->m_M)<=E_M_SMALL||
         fabs(stateData->m_M-TWOPI)<=E_M_SMALL)?true:false;

    const bool small_MH=
         (fabs(stateData->m_M)<=H_M_SMALL||
         fabs(stateData->m_M-TWOPI)<=H_M_SMALL)?true:false;

    if (stateData->m_a>=E_PARABOLIC2()||stateData->m_a<=H_PARABOLIC2())
      return true;

    if ((stateData->m_a>=E_PARABOLIC1()&&small_ME)||
         (stateData->m_a<=H_PARABOLIC1()&&small_MH))
           return true;

    if (stateData->m_e==1.0)
      return true;
    return false;
  };


  // ----------------------------------------------------------------------------
  bool CState::IsHyperbolic(void)const {
    const bool small_M=
         (fabs(stateData->m_M)<=H_M_SMALL||
         fabs(stateData->m_M-TWOPI)<=H_M_SMALL)?true:false;

    if (stateData->m_a>0.0)
      return false; // elliptic
    if (stateData->m_a<=H_PARABOLIC2()) // parabolic
      return false;
    if (stateData->m_a<=H_PARABOLIC1()&&small_M) // parabolic
      return false;

    return true;
  };

  // ---------------------------------------------------------------------------
  bool CState::isBodycentric(void) const
  {
    return (stateData->m_central_body != SUN);
  } // bool isBodycentric(void) const

  // ---------------------------------------------------------------------------
  bool CState::isHeliocentric(void) const
  {
    return (stateData->m_central_body == SUN);
  } // bool isHeliocentric(void) const

  // ---------------------------------------------------------------------------
  void CState::PrintHeader(ofstream&out, const char*comment /* = 0 */) const {
    if (!comment) {
      out<<"MJD " <<"[d], ";
      out<<"x " <<"[AU], ";
      out<<"y " <<"[AU], ";
      out<<"z " <<"[AU], ";
      out<<"dx " <<"[AU/d], ";
      out<<"dy " <<"[AU/d], ";
      out<<"dz " <<"[AU/d], ";
      out<<"r " <<"[AU], ";
      out<<"phi " <<"[rad], ";
      out<<"theta " <<"[rad], ";
      out<<"dr " <<"[AU/d], ";
      out<<"dphi " <<"[rad/d], ";
      out<<"dtheta " <<"[rad/d], ";
      out<<"a " <<"[AU], ";
      out<<"e, ";
      out<<"i " <<"[rad], ";
      out<<"Omega " <<"[rad], ";
      out<<"omega " <<"[rad], ";
      out<<"curlypi " <<"[rad], ";
      out<<"M " <<"[rad], ";
      out<<"E " <<"[rad], ";
      out<<"f " <<"[rad], ";
      out<<"L " <<"[rad], ";
      out<<"v " <<"[AU/d], ";
      out<<"Ekin " <<"[AU^2/d^2], ";
      out<<"Epot " <<"[AU^2/d^2], ";
      out<<"Etot " <<"[AU^2/d^2]";
    } else {
      out<<"MJD (" <<comment<<") [d], ";
      out<<"x (" <<comment<<") [AU], ";
      out<<"y (" <<comment<<") [AU], ";
      out<<"z (" <<comment<<") [AU], ";
      out<<"dx (" <<comment<<") [AU/d], ";
      out<<"dy (" <<comment<<") [AU/d], ";
      out<<"dz (" <<comment<<") [AU/d], ";
      out<<"r (" <<comment<<") [AU], ";
      out<<"phi (" <<comment<<") [rad], ";
      out<<"theta (" <<comment<<") [rad], ";
      out<<"dr (" <<comment<<") [AU/d], ";
      out<<"dphi (" <<comment<<") [rad/d], ";
      out<<"dtheta (" <<comment<<") [rad/d], ";
      out<<"a (" <<comment<<") [AU], ";
      out<<"e (" <<comment<<"), ";
      out<<"i (" <<comment<<") [rad], ";
      out<<"Omega (" <<comment<<") [rad], ";
      out<<"omega (" <<comment<<") [rad], ";
      out<<"curlypi (" <<comment<<") [rad], ";
      out<<"M (" <<comment<<") [rad], ";
      out<<"E (" <<comment<<") [rad], ";
      out<<"f (" <<comment<<") [rad], ";
      out<<"L (" <<comment<<") [rad], ";
      out<<"v (" <<comment<<") [AU/d], ";
      out<<"Ekin (" <<comment<<") [AU^2/d^2], ";
      out<<"Epot (" <<comment<<") [AU^2/d^2], ";
      out<<"Etot (" <<comment<<") [AU^2/d^2]";
    }
  } // end CState::PrintHeader


  // ----------------------------------------------------------------------------
  void CState::Print(ofstream&out) const {
    out<<setprecision(12);
    out<<stateData->m_mjd<<", ";
    out<<stateData->cartesian.position(0)<<", ";
    out<<stateData->cartesian.position(1)<<", ";
    out<<stateData->cartesian.position(2)<<", ";
    out<<stateData->cartesian.velocity(0)<<", ";
    out<<stateData->cartesian.velocity(1)<<", ";
    out<<stateData->cartesian.velocity(2)<<", ";
    out<<stateData->sphere.position(0)<<", ";
    out<<stateData->sphere.position(1)<<", ";
    out<<stateData->sphere.position(2)<<", ";
    out<<stateData->sphere.velocity(0)<<", ";
    out<<stateData->sphere.velocity(1)<<", ";
    out<<stateData->sphere.velocity(2)<<", ";
    out<<stateData->m_a<<", ";
    out<<stateData->m_e<<", ";
    out<<stateData->m_i<<", ";
    out<<stateData->m_Omega<<", ";
    out<<stateData->m_omega<<", ";
    out<<stateData->m_curlypi<<", ";
    out<<stateData->m_M<<", ";
    out<<stateData->m_E<<", ";
    out<<stateData->m_f<<", ";
    out<<stateData->m_L<<", ";
    out<<Velocity()<<", ";
    out<<KineticEnergy()<<", ";
    out<<PotentialEnergy()<<", ";
    out<<Energy();

  } // end CState::Print

} // of namespace SpaceLib

// end of file
