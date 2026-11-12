/*!
\file MLGlobal.cpp
\version 1.7
\brief General mathematical functions.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 09.05.2006 (AO)

\b Purpose:
   This file implements some widely used general mathematical functions.

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

//SpaceLib
#include <configSpaceLib.h>
#include <Global.h>
#include <MLConst.h>
#include <MLVector.h>

//standard C library
#include <math.h>

using namespace std;
using namespace Global;

namespace MathLib
{
  //---------------------------------------------------------------------------
  double AssertAngularRange2pi(double& angle)
  {
    while (angle >= TWOPI)
      angle -= TWOPI;
    while (angle < 0.0)
      angle += TWOPI;
    return angle;
  };
  //---------------------------------------------------------------------------
  double DifferenceAngle(const double left_angle, const double right_angle)
  {
    double diff = left_angle - right_angle;
    return AssertAngularRange(diff);
  }
  //---------------------------------------------------------------------------
  double signum(const double expression)
  {
    return (expression >= 0.0) ? 1.0 : -1.0;
  }

  //---------------------------------------------------------------------------
  double randSymInterval(const double interval)
  {
    return interval * 2.0 * (rand01() - 0.5);
  }
  //---------------------------------------------------------------------------
  double randInterval(const double min_val, const double max_val)
  {
    double interval_size = max_val - min_val;
    return min_val + interval_size * rand01();
  }
  //---------------------------------------------------------------------------
  int randIntInterval(const int min_val, const int max_val)
  {
    return min_val + static_cast<int>(rand01ex() * static_cast<float>(max_val - min_val + 1.0));

  }
  //---------------------------------------------------------------------------
  Vector randVec01(const int dim)
  {
    double *component = new double[dim];

    for (int i = 0; i < dim; i++)
      component[i] = rand01();

    Vector result(component, dim);

    delete[]component;

    return result;
  }
  //---------------------------------------------------------------------------
  Vector randVecSymInterval(const int dim, const double interval)
  {
    double *component = new double[dim];

    for (int i = 0; i < dim; i++)
      component[i] = randSymInterval(interval);

    Vector result(component, dim);

    delete[]component;

    return result;
  }
  //---------------------------------------------------------------------------
  Vector randVecInterval(const int dim, const double min_val, const double max_val)
  {
    double *component = new double[dim];

    for (int i = 0; i < dim; ++i)
      component[i] = randInterval(min_val, max_val);

    Vector result(component, dim);

    delete[]component;

    return result;
  }
  //---------------------------------------------------------------------------
  Vector randVecNorm(const int dim)
  {
    Vector vec(dim);

    vec = randVec01(dim);

    double vec_sum = 0.0;

    for (int i = 0; i < dim; i++)
      vec_sum += vec(i);

    return vec / vec_sum;
  }
  //---------------------------------------------------------------------------
  double randNormal(const double mean, const double sdev)
  {
    double x1, x2, w, y1;
    static double y2;
    static bool use_last = false;

    if (use_last) // use value from previous call
    {
      y1 = y2;
      use_last = false;
    }
    else
    {
      do
      {
        x1 = 2.0 * rand01() - 1.0;
        x2 = 2.0 * rand01() - 1.0;
        w = x1 * x1 + x2 * x2;
      } while (w >= 1.0);

      w = sqrt((-2.0 * log(w)) / w);
      y1 = x1 * w;
      y2 = x2 * w;
      use_last = true;
    }
    return mean + y1 * sdev;
  }
  //---------------------------------------------------------------------------
  bool FlipCoin(void)
  {
    if (rand01ex() >= 0.5)
      return true;
    else
      return false;
  }

  //---------------------------------------------------------------------------
  double acosh(const double x)
  {
    return log(x + sqrt(x * x - 1));
  }

  //---------------------------------------------------------------------------
  void AxisAngle2Quaternion(
    const double  v1
  , const double  v2
  , const double  v3
  , const double  angle
  ,       double& q1
  ,       double& q2
  ,       double& q3
  ,       double& q4)
  {
    const double sinhp( sin( angle / 2.0 ) );
    q1 = v1 * sinhp;
    q2 = v2 * sinhp;
    q3 = v3 * sinhp;
    q4 = cos( angle / 2.0 );
  }

  //---------------------------------------------------------------------------
  void AxisAngle2Quaternion(
    const double  v1
  , const double  v2
  , const double  v3
  , const double  angle
  , double q[4]
    )
  {
    AxisAngle2Quaternion(v1, v2, v3, angle, q[0], q[1], q[2], q[3]);
  };

  //---------------------------------------------------------------------------
  void AxisAngle2Quaternion(
    const double axis[3],
    const double angle
  ,       double& q1
  ,       double& q2
  ,       double& q3
  ,       double& q4)
  {
    AxisAngle2Quaternion(axis[0], axis[1], axis[2], angle,q1,q2,q3,q4);
  }

  //---------------------------------------------------------------------------
  void AxisAngle2Quaternion(
    const double axis[3]
   ,const double angle
   ,      double q[4]
  )
  {
    AxisAngle2Quaternion(axis[0],axis[1],axis[2], angle, q[0], q[1],
      q[2], q[3]);
  }


  //---------------------------------------------------------------------------
  double NormalizeQuaternion(double& q1, double& q2,
    double& q3, double& q4, const double eps)
  {
    const double norm_q(sqrt(q1*q1+q2*q2+q3*q3+q4*q4));
    if (fabs(1.0-norm_q) > eps)
    {
      q1/=norm_q;
      q2/=norm_q;
      q3/=norm_q;
      q4/=norm_q;
      return norm_q;
    }
    else return 1.0;
  }





} // end of namespace MathLib
