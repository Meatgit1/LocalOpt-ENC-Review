/*!
\file Kepler.cpp
\version 1.7
\brief Miscellaneous functions for Keplerian elements.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Tanja D. Schmidt <tanja.schmidt@dlr.de>c (TS)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 09.05.2006 (AO)

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
   This code is based on the code from:
   Montenbruck, Gill: Satellite Orbits. Springer, 2000 and
   Montenbruck, Pfleger: Astronomie mit dem Personal Computer (3. Aufl.). Springer, 1999

\b Copyright:
   (c) 1999-2006 Oliver Montenbruck, Eberhard Gill, Thomas Pfleger
*/

//SpaceLib
#include <configSpaceLib.h>

// standard C library
#include <math.h>

//SpaceLib
#include <Global.h>
#include <MLConst.h>
#include <MLGlobal.h>
#include <Kepler.h>

using namespace std;
using namespace Global;
using namespace MathLib;

namespace SpaceLib {

  const int MAXIT = 50; ///< tbd

  //--------------------------------------------------------------------------
  void Stumpff(double E2, double &c1, double &c2, double &c3) {

    double n, add;

    c1 = c2 = c3 = 0.0;
    add = n = 1.0;

    do {
      c1 += add;
      add /= (2.0 * n);
      c2 += add;
      add /= (2.0 * n + 1.0);
      c3 += add;
      add *= -E2;
      n += 1.0;
    } while (fabs(add) >= EPS);
  };



  //--------------------------------------------------------------------------
  double PeriapsisDistance(double orbitParameter, double e) {
    return orbitParameter / (1.0 + e);
  }


  //--------------------------------------------------------------------------
  double OrbitParameter(double q, double e) {
    return q * (1.0 + e);
  }


  //--------------------------------------------------------------------------
  double EccentricAnomaly(double M, const double e) {

    double E, f;

    AssertAngularRange2pi(M);
    E = (e < 0.8) ? M : PI;     // start value of iteration
    int loop = 0;

    do {
      f = E - e * sin(E) - M;
      E = E - f / (1.0 - e * cos(E));
      if (++loop == MAXIT) {
        cerr << "convergence problems in EccentricAnomaly()" << endl;
        break;
      }
    } while (fabs(f) > EPS);

    return E;
  };


  //--------------------------------------------------------------------------
  double HyperbolicAnomaly(double Mh, double e) {

    double H = pow(6.0 * Mh, 1.0 / 3.0);  // start value for iteration

    if (Mh < 0.0) H = -H;

    int loop = 0;
    double H_old;

    do {
      H_old = H;
      H = H - (e * sinh(H) - H - Mh) / (e * cosh(H) - 1.0);
      if (++loop == MAXIT) {
        cerr << "convergence problems in HyperbolicAnomaly()" << endl;
        break;
      }
    } while (fabs(H - H_old) > EPS);

    return H;
  };


  //--------------------------------------------------------------------------
  void CalcAnomaliesM(const double a,
                      const double M,
                      const double e,
                      double &EH,
                      double &f) {
    if (a > 0.0) {
      EH = EccentricAnomaly(M, e);
      AssertAngularRange2pi(EH);
      f = 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(0.5 * EH));
      AssertAngularRange2pi(f);
    }
    else if (a < 0.0) {
      EH = HyperbolicAnomaly(M, e);
      f = 2.0 * atan(sqrt((e + 1.0) / (e - 1.0)) * tanh(0.5 * EH));
      AssertAngularRange2pi(f);
    }
    else
      FatalError("CalcAnomaliesM(): semi-major axis is exactly zero");

  };


  //--------------------------------------------------------------------------
  void CalcAnomaliesf(const double a,
                      const double f,
                      const double e,
                      double &EH,
                      double &M) {

    EH = 2.0 * atan(sqrt((1.0 - e) / (1.0 + e)) * tan(0.5 * f));
    AssertAngularRange2pi(EH);

    if (a > 0.0) {
      M = EH - e * sin(EH);
    }
    else if (a < 0.0) {
      M = e * sinh(EH) - EH;
    }
    else
      FatalError("CalcAnomaliesf(): semi-major axis is exactly zero");
  };


  //--------------------------------------------------------------------------
  double HyperbolicMeanAnomaly(double t0, double t, double gm, double a) {
    return sqrt(gm * (t - t0) / pow(fabs(a), 3.0));
  }

  //--------------------------------------------------------------------------
  bool AssertElevationRange(const double theta) {
    if (!(theta >= -HALFPI && theta <= HALFPI))
    {
      //cout << "\ntheta:" << theta <<endl;
      return false;
      FatalError("AssertElevationRange() => elevation out of range");
    }
    return true;
  };


  //--------------------------------------------------------------------------
  void AssertEccentricityRange(const double e) {
    if (!(e >= 0.0))
      FatalError("AssertEccentricityRange() => eccentricity out of range");
  };


  //--------------------------------------------------------------------------
  void AssertInclinationRange(const double i) {
    if (!(i >= 0.0 && i <= PI))
      FatalError("AssertInclinationRange() => inclination out of range");
  };

}                               // namespace SpaceLib

// end of file

