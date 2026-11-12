/*!
\file MLGlobal.h
\version 1.7
\brief General mathematical functions.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 09.05.2006 (AO)

\b Purpose:
   This file provides the interface for some widely used general mathematical functions.

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

#ifndef MATHGLOBAL_H
#define MATHGLOBAL_H

#include <Global.h>

namespace MathLib
{
  class Vector;

  /** \brief Assures that the given angle is in the interval ]-pi,pi]
      Modifies argument angle!
      \param  angle angle that should be in the desired interval
      \return given angle but between ]-pi,pi] */
  template<class T>
  inline 
  T AssertAngularRange(T& angle)
  {
    while (angle > PI)   angle -= TWOPI;
    while (angle <= -PI) angle += TWOPI;
    return angle;
  }

  /** \brief Assures that the given angle is in the range ]0,2*pi]
      Modifies argument angle!
      \param angle angle that should be in the desired interval
      \return given angle but in between ]0,2+pi] */
  double AssertAngularRange2pi(double& angle);

  /** \brief Calculates the difference of two angles.
      Subtracts angle two from angle one and assures that the difference is in the range ]-pi,pi].
      \param left_angle first angle
      \param right_angle second angle
      \return difference angle
      \sa AssertAngularRange */
  double DifferenceAngle(const double left_angle, const double right_angle);

  /** \brief Returns the sign of the given value.
      Tests, whether the given expression is negative or not.
      \param expression value for which the sign has do be determined
      \return -1.0 if the expression is negative and 1.0 if it is positive */
  double signum(const double expression);

  /** \brief Returns the given number bounded by the given constraints.
      The function returns the lower boundary if the value to test
      is smaller than the lower boundary and the upper boundary if the
      value is greater than the upper boundary. Otherwise the function
      returns the given argument.
      \param what value to test against the boundaries
      \param l lower bound
      \param u upper bound
      \return what, upper bound, or lower bound */
  template<class T>
    inline const T& bounded(const T& what, const T& l, const T& u)
  {
    if (l < u)
    {
      if (what < l) return l;
      else if (what > u) return u;
      else return what;
    }
    else
    {
      if (what < u) return u;
      else if (what > l) return l;
      else return what;
    }
  } // const T& bounded(const T&, const T&, const T&)

  /// Returns a (pseudo-) random number from the interval [0,1]
  inline double rand01(void)
  {
    return (static_cast < double >(rand()) / static_cast < double >(RAND_MAX));
  }

  /// Returns a (pseudo-) random number from the intervall [0,1[
  inline double rand01ex(void)
  {
    double result;
    do {} while ((result = rand01()) == 1.0);
    return result;
  }

  /** \brief Returns a (pseudo-)random number from the interval [-interval, +interval].
      \param interval half the size of the 0-centered interval
      \return random number from the interval [-interval, +interval] */
  double randSymInterval(const double interval);

  /** \brief Returns a (pseudo-)random number from the interval [min_val, max_val].
      \param min_val lower boundary
      \param max_val upper boundary
      \return random number from the interval [min_val, max_val] */
  double randInterval(const double min_val, const double max_val);

  /** \brief Returns a (pseudo-)random integer number from the interval [min_val,max_val].
      \param min_val lower boundary
      \param max_val upper boundary
      \return random integer number from the interval [min_val, max_val] */
  int randIntInterval(const int min_val, const int max_val);

  /** \brief Returns a (pseudo-) random vector of dimension dim, where all components are in the interval [0,1].
      \param dim dimension of the target vector
      \return random vector (all components are in the interval [0,1]) */
  Vector randVec01(const int dim);

  /** \brief Returns a (pseudo-) random vector of dimension dim, where all components are in the interval [-interval, +interval].
      \param dim dimension of the vector
      \param interval half the size of the 0-centered interval
      \return random vector (all components are in the interval [-interval, +interval])
      \sa randVec01, randVecInterval, randVecNorm */
  Vector randVecSymInterval(const int dim, const double interval);

  /** \brief Returns a (pseudo-)random vector of dimension dim, where all components are in the interval [min_val, max_val].
      \param dim dimension of the vector
      \param min_val lower boundary for vector components
      \param max_val upper boundary for vector components
      \return random vector (all components are in the interval [min_val, max_val])
      \sa randVec01, randVecSymInterval, randVecNorm */
  Vector randVecInterval(const int dim, const double min_val, const double max_val);

  /** \brief Returns a (pseudo-)random vector of dimension dim, where all components are in the interval [0,1] and the sum of all components is 1.
      \param dim dimension of the target vector
      \return random vector (all components are in the interval [0,1] and the sum of all components is 1)
      \sa randVec01, randVecSymInterval, randVecIntervall */
  Vector randVecNorm(const int dim);

  /** \brief Returns a (pseudo-)random number from the normal distribution.
      Implements the polar form of the Box-Muller transformation.
      (c) Copyright 1994, Everett F. Carter Jr.
      Permission is granted by the author to use this software for any
      application provided this copyright notice is preserved.
      \param mean Mean
      \param sdev Standard deviation
      \return random number from the normal distribution with mean and standard deviation sdev */
  double randNormal(const double mean, const double sdev);

  /** \brief Simulates the flipping of an unbiased coin.
      \return true or false with a 50/50 probability */
  bool FlipCoin(void);

  /** \brief Simulates the flipping of a biased coin.
      \param true_probability The probabilty for the true to be be returned
      \return true with (true_probability)/(1-true_probability) probability (or false) */
  template<class T>
  bool Flip(const T true_probability)
  {
    if (rand01ex() < true_probability) return true;else return false;
  }

  /// Computes the arcus cosinus hyperbolicus (Microsoft C++ 6.0).
  double acosh(const double x);

  /** t.b.d. */
  void AxisAngle2Quaternion(
    const double  v1
  , const double  v2
  , const double  v3
  , const double  phi
  ,       double& q1
  ,       double& q2
  ,       double& q3
  ,       double& q4);

  /** t.b.d. */
  void AxisAngle2Quaternion(
    const double  v1
  , const double  v2
  , const double  v3
  , const double  phi
  , double q[4]
  );

  /** t.b.d. */
  void AxisAngle2Quaternion(
    const double axis[3],
    const double angle
  ,       double& q1
  ,       double& q2
  ,       double& q3
  ,       double& q4);

  /** t.b.d. */
  void AxisAngle2Quaternion(
    const double axis[3]
   ,const double angle
   ,      double q[4]
  );

  /**
    \brief Correct attitude quaternion to norm of 1 (references to doubles).
    \param q1 attitude quaternion component 1 (first vector component)
    \param q2 attitude quaternion component 2 (second vector component)
    \param q3 attitude quaternion component 3 (third vector component)
    \param q4 attitude quaternion component 4 (fourth vector component)
    \param eps Desired maximal deviation from unit length.
    \return Corrected attitude quaternion parameter.
  */
  double NormalizeQuaternion(
    double& q1
   ,double& q2
   ,double& q3
   ,double& q4
   ,const double eps = 1.0E-6);


  /**
    \brief Correct attitude quaternion to norm of 1 (array of doubles).
    \param q attitude quaternion
    \param eps Desired maximal deviation from unit length.
    \return Corrected attitude quaternion parameter.
  */
  inline double NormalizeQuaternion(double q[4],const double eps = 1.0E-6)
  {
    return NormalizeQuaternion(q[0], q[1], q[2], q[3], eps);
  }


} // end of namespace Mathlib

#endif // of MATHGLOBAL_H

