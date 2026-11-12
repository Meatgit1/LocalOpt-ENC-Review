/*!
\file MLGJ4P.cpp
\version 1.7
\brief A class for the numerical integration if 2nd order DES using the Gauss-Jackson-algorithm (GJ4P).
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\date 06.01.2005 (AO)

\b Purpose:
   This module of SpaceLib provides a class that is able to integrate
   a second order DES with the algorithm of Gauss-Jackson by using a
   4th order predictor Runge-Kutta scheme. DES of the first order can not
   by solved with this class. The handling is similar to the classes
   declarated in 'MLDE.h'.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.
   This code is based on the code from:
   Montenbruck, Gill: Satellite Orbits. Springer, 2000

\b Copyright:
   (c) 1999-2004 Oliver Montenbruck, Eberhard Gill
*/

//SpaceLib
#include <configSpaceLib.h>
#include <MLMatrix.h>
#include <MLGJ4P.h>

using namespace std;

namespace MathLib {


  // 4th order Runge-Kutta step for 2nd order differential equation
  //------------------------------------------------------------------------------
  void GJ4P::RK4(double &t, Vector & r, Vector & v, double h)
  {
    Vector v_1, v_2, v_3, v_4;
    Vector a_1, a_2, a_3, a_4;

      v_1 = v;
      f(t, r, v_1, a_1, pAux);
      v_2 = v + (h / 2.0) * a_1;
      f(t + h / 2.0, r + (h / 2.0) * v_1, v_2, a_2, pAux);
      v_3 = v + (h / 2.0) * a_2;
      f(t + h / 2.0, r + (h / 2.0) * v_2, v_3, a_3, pAux);
      v_4 = v + h * a_3;
      f(t + h, r + h * v_3, v_4, a_4, pAux);

      t = t + h;
      r = r + (h / 6.0) * (v_1 + 2.0 * v_2 + 2.0 * v_3 + v_4);
      v = v + (h / 6.0) * (a_1 + 2.0 * a_2 + 2.0 * a_3 + a_4);

  };                            // end of GJ4P::RK4



  // Initialization of backwards differences from initial conditions
  //------------------------------------------------------------------------------
  void GJ4P::Init(double t_0, const Vector & r_0, const Vector & v_0, double h_)
  {

    const int m = 4;            // Order of method

    // Coefficients gamma/delta of 1st/2nd order Moulton/Cowell corrector method

    const double gc[m + 1] = { +1.0, -1 / 2.0, -1 / 12.0, -1 / 24.0, -19 / 720.0 };
    const double dc[m + 2] = { +1.0, -1.0, +1 / 12.0, 0.0, -1 / 240.0, -1 / 240.0 };

    int i, j;
    double t = t_0;
    Vector r = r_0;
    Vector v = v_0;

    // Save step size

    h = h_;

    // Create table of accelerations at past times t-3h, t-2h, and t-h using
    // RK4 steps

    f(t, r, v, D[0], pAux);     // D[i]=a(t-ih)
    for (i = 1; i <= m - 1; i++) {
      RK4(t, r, v, -h);
      f(t, r, v, D[i], pAux);
    };

    // Compute backwards differences

    for (i = 1; i <= m - 1; i++)
      for (j = m - 1; j >= i; j--)
        D[j] = D[j - 1] - D[j];

    // Initialize backwards sums using 4th order GJ corrector

    S1 = v_0 / h;
    for (i = 1; i <= m; i++)
      S1 -= gc[i] * D[i - 1];
    S2 = r_0 / (h * h) - dc[1] * S1;
    for (i = 2; i <= m + 1; i++)
      S2 -= dc[i] * D[i - 2];

  }                             // end of GJ4P::Init


  // Step from t to t+h
  //------------------------------------------------------------------------------
  void GJ4P::Step(double &t, Vector & r, Vector & v)
  {

    const int m = 4;            // Order of method

    // Coefficients gamma/delta of 1st/2nd order Bashforth/Stoermr predictor

    const double gp[m + 1] = { +1.0, +1 / 2.0, +5 / 12.0, +3 / 8.0, +251 / 720.0 };
    const double dp[m + 2] = { +1.0, 0.0, +1 / 12.0, +1 / 12.0, +19 / 240.0, +3 / 40.0 };

    int i;

    // 4th order predictor

    r_p = dp[0] * S2;
    for (i = 2; i <= m + 1; i++)
      r_p += dp[i] * D[i - 2];
    r_p = (h * h) * r_p;
    v_p = gp[0] * S1;
    for (i = 1; i <= m; i++)
      v_p += gp[i] * D[i - 1];
    v_p = h * v_p;

    // Update backwards difference table

    f(t + h, r_p, v_p, d[0], pAux); // Acceleration at t+h
    for (i = 1; i <= m - 1; i++)
      d[i] = d[i - 1] - D[i - 1]; // New differences at t+h
    for (i = 0; i <= m - 1; i++)
      D[i] = d[i];              // Update differences
    S1 += d[0];
    S2 += S1;                   // Update sums

    // Update independent variable and solution

    t = t + h;
    r = r_p;
    v = v_p;

  }                             // end of GJ4P::Step

}                               // end of namespace MathLib

// end of implementation

