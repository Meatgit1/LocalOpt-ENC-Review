/*!
\file MLGJ4P.h
\version 1.7
\brief A class for the numerical integration if 2nd order DES using the Gauss-Jackson-algorithm (GJ4P).
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\date 27.10.2004 (BD)

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

#ifndef GJ4P_H
#define GJ4P_H

#include <MLVector.h>

namespace MathLib
{

  class Matrix;
  class RK4;


  /// Function prototype for second order differential equations
  /// void f (double t, const Vector& r, const Vector& v, Vector& a)

  /**
     \param t The independent variable.
     \param r The variable vector (position).
     \param v The first derivation of r (velocity).
     \param a The second derivation of r = (first derivation of v = acceleration).
     \param pAux A pointer to auxiliary data.
   */
  typedef void (*GJ4Pfunct) (double t, const Vector & r, const Vector & v, Vector & a, void *pAux);


  /// The integrator class of Gauss-Jackson algorithm.

  /**
     Although the inteface of GJ4P looks quite the same as the integrator classes
     defined in 'MLDE.h' there is one essential difference. It can only handle
     differential equation system (DES) of the order two where the others can only
     handle DES of the order one. This difference shows up in the declaration of the
     used DES-function (GJ4Pfunct).
   */
  class GJ4P
  {

  private:

    // 4th order Runge-Kutta step
    void RK4(double &t,         // Independent variable; updated by t+h
      Vector & r,               // Value of r(t); updated by r(t+h)
      Vector & v,               // Value of v(t)=dr/dt(t); updated by v(t+h)
      double h                  // Step size
      );

    // Elements
    int n_eqn;                  // Dimension
    GJ4Pfunct f;                // Differential equation
    double h;                   // Step size
    void *pAux;                 // Pointer to auxiliary data requird by f
    Vector S2, S1;              // First and second sum of acceleration
    Vector D[4];                // Backward differences of acceleration at t
    Vector d[4];                // Backward differences of acceleration at t+h
    Vector r_p, v_p;            // Predictor

  public:

    /**
    * \brief Constructor
    * \param f_     The actual DES.
    * \param n_eqn_ The DES' dimension = the number of equations.
    * \param pAux_   A pointer to auxiliary data.
    *
    *  The interface of the function that incorporates the DES to integrate.
    *  This function (in fact a pointer to such a function) has to be provided
    *  by the user of this integrator class, because he is the only one who
    *  knows, or at least should know about the mathematical nature of his
    *  problem.
    */
    GJ4P(GJ4Pfunct f_, int n_eqn_, void *pAux_):n_eqn(n_eqn_), f(f_), pAux(pAux_)
    {};


    /// Initialization of integrator's start values.

    /**
       \param t_0 The initial value of independent variable.
       \param r_0 The initial value r_0=r(t_0).
       \param v_0 The initial value v_0=dr/dt(t_0).
       \param h_  The step size.

       Except the step size h_ all arguments of this method are the initial
       conditions of a so called 'Anfangswert' - problem and like the DES-function
       have to be provided by the user.
     */
    void Init(double t_0, const Vector & r_0, const Vector & v_0, double h_);


    /// One integration step (of given step size).

    /**
       \param t The independent variable; updated by t+h.
       \param r The value of r(t); updated by r(t+h).
       \param v The value of v(t)=dr/dt(t); updated by v(t+h).

       This method actually does the integration of the DES with the given
       step size. After the integration the parameters are set with the new
       values, that means t, r and v will be updated with reference to t+h.
     */
    void Step(double &t, Vector & r, Vector & v);
  };

}                               // of namespace MathLib

#endif // of #ifndef GJ4P_H

// end of file

