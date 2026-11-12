/*!
\file MLVector.cpp
\version 1.7
\brief A class for mathematical vectors and functions for vectors (Vector).
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 09.05.2006 (AO)

\b Purpose:
   This file implements the class Vector.
   Vector implements 'double' vectors.

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
   (c) 1999-2006 Oliver Montenbruck, Eberhard Gill, Andreas Ohndorf, Bernd Dachwald
*/

//SpaceLib
#include <configSpaceLib.h>

//STL
#include <iomanip>

//standard C library
#include <math.h>
#include <cassert>

//SpaceLib
#include <Global.h>
#include <MLGlobal.h>
#include <MLVector.h>

using namespace std;
using namespace Global;

namespace MathLib
{
//===== member functions ===================================================

  Vector::Vector()
  : n(0)
  {
    v = NULL;
  }

//--------------------------------------------------------------------------
  Vector::Vector(int dim)
  : n(dim)
  {
    v = new double[dim];
    for (int i = 0; i < dim; ++i)
      v[i] = 0.0;
  }

//--------------------------------------------------------------------------
  Vector::Vector(const Vector & vec)
  : n(vec.n)
  {
    v = new double[vec.n];
    for (int i = 0; i < vec.n; ++i)
      v[i] = vec.v[i];
  }

//--------------------------------------------------------------------------
  Vector::Vector(const double *p, int dim)
  : n(dim)
  {
    v = new double[dim];
    for (int i = 0; i < dim; ++i)
      v[i] = p[i];
  }

//--------------------------------------------------------------------------
  Vector::Vector(double v0, double v1)
  : n(2)
  {
    v = new double[2];
    v[0] = v0;
    v[1] = v1;
  }

//--------------------------------------------------------------------------
  Vector::Vector(double v0, double v1, double v2)
  : n(3)
  {
    v = new double[3];
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
  }

//--------------------------------------------------------------------------
  Vector::Vector(double v0, double v1, double v2, double v3)
  : n(4)
  {
    v = new double[4];
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;
  }

//--------------------------------------------------------------------------
  Vector::Vector(double v0, double v1, double v2, double v3, double v4)
  : n(5)
  {
    v = new double[5];
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;
    v[4] = v4;
  }

//--------------------------------------------------------------------------
  Vector::Vector(double v0, double v1, double v2, double v3, double v4, double v5)
  : n(6)
  {
    v = new double[6];
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;
    v[4] = v4;
    v[5] = v5;
  }

//--------------------------------------------------------------------------
  Vector::Vector(double v0, double v1, double v2, double v3, double v4, double v5, double v6)
  : n(7)
  {
    v = new double[7];
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;
    v[4] = v4;
    v[5] = v5;
    v[6] = v6;
  }

//--------------------------------------------------------------------------
  Vector::Vector(double v0, double v1, double v2, double v3, double v4, double v5, double v6, double v7)
  : n(8)
  {
    v = new double[8];
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;
    v[4] = v4;
    v[5] = v5;
    v[6] = v6;
    v[7] = v7;
  }

//--------------------------------------------------------------------------
  Vector::Vector(double v0, double v1, double v2, double v3, double v4, double v5, double v6, double v7, double v8)
  : n(9)
  {
    v = new double[9];
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;
    v[4] = v4;
    v[5] = v5;
    v[6] = v6;
    v[7] = v7;
    v[8] = v8;
  }

  //--------------------------------------------------------------------------
  Vector::Vector(double v0, double v1, double v2, double v3, double v4, double v5, double v6, double v7, double v8, double v9)
    : n(10)
  {
    v = new double[10];
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;
    v[4] = v4;
    v[5] = v5;
    v[6] = v6;
    v[7] = v7;
    v[8] = v8;
    v[9] = v9;
  }


//-----------------------------------------------------------------------------
  const MathLib::Vector& Vector::Set(const MathLib::Vector& vec)
  {
    if (this == &vec) return (*this);        // check identity
      
    if (v == NULL)                           // allocate memory if necessary
    {
      n = vec.n;
      v = new double[n];
    }
    else if (vec.n != Vector::n)             // in case of not equal dimensions
    {
      delete[] v;                              // free old data vec
      v = new double[n=vec.n];               // and allocate new space
    }

    for (int i = 0; i < n; ++i) v[i] = vec.v[i]; // copy values

    return *this;
  }



//--------------------------------------------------------------------------
  Vector Vector::Slice(int first, int last) const
  {
//#define development_version
//#ifdef development_version
    assert(first >= 0 && last >= 0 && first < n && last < n);
    if (!(first >= 0 && last >= 0 && first < n && last < n))
    {
      flush(cout << "first:"<<first<<endl);
      flush(cout << "last :"<<last <<endl);
      flush(cout << "n    :"<<n<<endl);
      //FatalError("Vector::Slice() => invalid index");
    }
//#endif

    int dim = last-first+1;

#ifdef development_version
    assert(dim >= 1);
//      FatalError("Vector::Slice() => invalid dimension");
#endif

    Vector result(dim);

    for (int i = first; i <= last; ++i)
      result.v[i-first] = v[i];

    return result;
  }

//--------------------------------------------------------------------------
  const Vector& Vector::operator=(const double value)
  {
    for (int i = 0; i < n; ++i)
      v[i] = value;

    return *this;
  }

//--------------------------------------------------------------------------
  const Vector& Vector::operator=(const Vector & vec)
  {
    if (this == &vec) // check identity
      return (*this);

    if (v == NULL) // allocate memory if necessary
    {
      n = vec.n;
      v = new double[n];
    };

#ifdef development_version
    if (n != vec.n)
      FatalError("Vector::operator=() => incompatible sizes");
#endif

    for (int i = 0; i < n; ++i)
      v[i] = vec.v[i];

    return *this;
  }

//--------------------------------------------------------------------------
  bool Vector::operator==(const Vector & vec) const
  {
    if (n != vec.n)
      return false;

    for (int i = 0; i < n; ++i)
      if (v[i] != vec.v[i])
        return false;

    return true;
  }

//--------------------------------------------------------------------------
  bool Vector::operator!=(const Vector & vec)
  {
    if (n != vec.n)
      return true;

    for(int i = 0; i < n; ++i)
      if (v[i] != vec.v[i])
        return true;

    return false;
  }

//--------------------------------------------------------------------------
  double Vector::operator()(int i) const
  {
#ifdef development_version
    assert((i<n) && (i>=0));
#endif
    return v[i];
  }

//--------------------------------------------------------------------------
  double &Vector::operator()(int i)
  {
#ifdef development_version
    assert((i<n) && (i>=0));
#endif

    return v[i];
  }

//--------------------------------------------------------------------------
  double Vector::operator[](int i) const
  {
#ifdef development_version
    assert((i<n) && (i>=0));
#endif

    return v[i];
  }

//--------------------------------------------------------------------------
  double &Vector::operator[](int i)
  {
#ifdef development_version
    assert((i<n) && (i>=0));
#endif

    return v[i];
  }

//--------------------------------------------------------------------------
  double Vector::getComponent(const int i) const
  {
#ifdef development_version
    assert((i<n) && (i>=0));
#endif

    return v[i];
  }

//--------------------------------------------------------------------------
  void Vector::setComponent(const int i, const double value)
  {
#ifdef development_version
    assert((i<n) && (i>=0));
#endif

    v[i] = value;
  }

//--------------------------------------------------------------------------
  Vector & Vector::operator+=(const Vector & vec)
  {
#ifdef development_version
    assert(n == vec.n);
#endif

    for (int i = 0; i < n; ++i)
      v[i] += vec.v[i];

    return *this;
  }

//--------------------------------------------------------------------------
  Vector & Vector::operator-=(const Vector & vec)
  {
#ifdef development_version
    assert(n==vec.n);
#endif

    for (int i = 0; i < n; ++i)
      v[i] -= vec.v[i];

    return *this;
  }

//--------------------------------------------------------------------------
  Vector & Vector::operator*=(const double factor)
  {
    for (int i = 0; i < n; ++i)
      v[i] *= factor;

    return *this;
  }

//--------------------------------------------------------------------------
  Vector & Vector::operator/=(const double divisor)
  {
    for (int i = 0; i < n; ++i)
      v[i] /= divisor;

    return *this;
  }

//===== friend functions ===================================================

  //---------------------------------------------------------------------------
  void AxisAngle2Quaternion(
    const Vector& axis,
    const double  angle,
    Vector& q
  )
  {
    Vector v_q(4);
    AxisAngle2Quaternion(axis.v[0],axis.v[1],axis.v[2], angle, q.v);
  }


  //---------------------------------------------------------------------------
  void AxisAngle2Quaternion(
      const Vector& axis,
      const double  angle,
      double q[4]
  )
  {
    AxisAngle2Quaternion(axis.v, angle, q);
  }


  //---------------------------------------------------------------------------
  Vector Stack(const Vector & a, const Vector & b)
  {
    Vector result(a.n + b.n);

    int i;
    for (i = 0; i < a.n; ++i)
      result.v[i] = a.v[i];
    for (i = 0; i < b.n; ++i)
      result.v[a.n + i] = b.v[i];

    return result;
  }

//--------------------------------------------------------------------------
  Vector VecPolar(double azim, double elev, double r)
  {
    return Vector(r * cos(azim) * cos(elev), r * sin(azim) * cos(elev), r * sin(elev));
  }

//--------------------------------------------------------------------------
  double Dot(const Vector & a, const Vector & b)
  {
#ifdef development_version
    if (a.n != b.n)
      FatalError("Dot() => incompatible sizes");
#endif

    double result = 0.0;

    for (int i = 0; i < a.n; ++i)
      result += a.v[i] * b.v[i];

    return result;
  }

//--------------------------------------------------------------------------
  Vector Cross(const Vector & left, const Vector & right)
  {
#ifdef development_version
    if (left.n != 3 || right.n != 3)
      FatalError("Cross() => invalid dimension");
#endif

    Vector result(3);

    result.v[0] = left.v[1] * right.v[2] - left.v[2] * right.v[1];
    result.v[1] = left.v[2] * right.v[0] - left.v[0] * right.v[2];
    result.v[2] = left.v[0] * right.v[1] - left.v[1] * right.v[0];

    return result;
  }

//--------------------------------------------------------------------------
  double Norm(const Vector & vec)
  {
    return sqrt(Dot(vec, vec));
  }


  // Calculate axis/angle representation from quaternion
  //---------------------------------------------------------------------------
  void Quaternion2AxisAngle(const double  q1  // = e1*sin(phi/2)
                          , const double  q2  // = e2*sin(phi/2)
                          , const double  q3  // = e3*sin(phi/2)
                          , const double  q4  // = cos(phi/2)
                          ,       Vector& axis
                          ,       double& angle) 
  {
    const double sinhp( sqrt( 1.0 - q4 * q4 ) );
    if ( sinhp > 0.0 )
    {
      angle = 2.0 * acos( q4 );                      // return rotation angle
      axis[0] = q1/sinhp;
      if (axis[0]>1.0)
      {
        axis[0] = 1.0;axis[1] = 0.0;axis[2] = 0.0;
        return;
      }
      axis[1] = q2/sinhp;
      if (axis[1]>1.0)
      {
        axis[0] = 0.0;axis[1] = 1.0;axis[2] = 0.0;
        return;
      }
      axis[2] = q3/sinhp;
      if (axis[2]>1.0)
      {
        axis[0] = 0.0;axis[1] = 0.0;axis[2] = 1.0;
        return;
      }     
    }
    else
    {
      axis.Set(Vector(3)); // initialize axis with 3-dimensinal zero values vector
      angle = 0.0; // set rotation angle to zero
    }
  }



//--------------------------------------------------------------------------
  Vector operator*(double value, const Vector & vec)
  {
    Vector result(vec.n);

    for (int i = 0; i < vec.n; ++i)
      result.v[i] = value * vec.v[i];

    return result;
  }

//--------------------------------------------------------------------------
  Vector operator*(const Vector & vec, double value)
  {
    return value * vec;
  }

//--------------------------------------------------------------------------
  Vector operator/(const Vector & vec, double value)
  {
    Vector result(vec.n);

    for (int i = 0; i < vec.n; ++i)
      result.v[i] = vec.v[i] / value;

    return result;
  }

//--------------------------------------------------------------------------
  Vector operator-(const Vector & vec)
  {
    Vector result(vec.n);

    for (int i = 0; i < vec.n; ++i)
      result.v[i] = -vec.v[i];

    return result;
  }

//--------------------------------------------------------------------------
  Vector operator+(const Vector & a, const Vector & b)
  {
#ifdef development_version
    if (a.n != b.n)
      FatalError("Vector operator+(Vector&) => incompatible sizes");
#endif

    Vector result(a.n);

    for (int i = 0; i < a.n; ++i)
      result.v[i] = a.v[i] + b.v[i];

    return result;
  }

//--------------------------------------------------------------------------
  Vector operator-(const Vector & left, const Vector & right)
  {
#ifdef development_version
    if (left.n != right.n)
      FatalError("Vector operator-(Vector&): incompatible sizes");
#endif

    Vector result(left.n);

    for (int i = 0; i < left.n; ++i)
      result.v[i] = left.v[i] - right.v[i];

    return result;
  }

//--------------------------------------------------------------------------
  std::ostream & operator<<(std::ostream & os, const Vector & vec)
  {
    for (int i = 0; i < vec.n; ++i)
    {
      os << setw(os.width()) << setprecision(os.precision())<< vec.v[i];
      if (i < vec.n - 1)
        os << ", ";
    }
    return os;
  }

  //----------------------------------------------------------------------------
  void MultiplyQuaternion(
    const Vector& q1,
    const int     i1,
    const Vector& q2,
    const int     i2,
    Vector& r,
    const int     ir
    )
  {
    if (r.n>(ir+3) && q1.n>(i1+3) && q2.n>(i2+3)) 
    {
      r.v[ir  ] = q2.v[i2+3]*q1.v[i1  ]+q1.v[i1+3]*q2.v[i2  ]+q1.v[i1+1]*q2.v[i2+2]-q1.v[i1+2]*q2.v[i2+1];
      r.v[ir+1] = q2.v[i2+3]*q1.v[i1+1]+q1.v[i1+3]*q2.v[i2+1]+q1.v[i1+2]*q2.v[i2  ]-q1.v[i1  ]*q2.v[i2+2];
      r.v[ir+2] = q2.v[i2+3]*q1.v[i1+2]+q1.v[i1+3]*q2.v[i2+2]+q1.v[i1  ]*q2.v[i2+1]-q1.v[i1+1]*q2.v[i2  ];
      r.v[ir+3] = q2.v[i2+3]*q1.v[i1+3]-q1.v[i1  ]*q2.v[i2  ]-q1.v[i1+1]*q2.v[i2+1]-q1.v[i1+2]*q2.v[i2+2];
    }
  }

  //----------------------------------------------------------------------------
  void CalcRelativeQuaternion(
    const Vector& q1,
    const int     i1,
    const Vector& q2,
    const int     i2,
          Vector& dq,
    const int     idq
  )
  {
    if (dq.n>(idq+3) && q1.n>(i1+3) && q2.n>(i2+3)) 
    {
      dq.v[idq  ] = q2.v[i2+3]*q1.v[i1  ]-q1.v[i1+3]*q2.v[i2  ]+q1.v[i1+1]*q2.v[i2+2]-q1.v[i1+2]*q2.v[i2+1];
      dq.v[idq+1] = q2.v[i2+3]*q1.v[i1+1]-q1.v[i1+3]*q2.v[i2+1]+q1.v[i1+2]*q2.v[i2  ]-q1.v[i1  ]*q2.v[i2+2];
      dq.v[idq+2] = q2.v[i2+3]*q1.v[i1+2]-q1.v[i1+3]*q2.v[i2+2]+q1.v[i1  ]*q2.v[i2+1]-q1.v[i1+1]*q2.v[i2  ];
      dq.v[idq+3] = q2.v[i2+3]*q1.v[i1+3]+q1.v[i1  ]*q2.v[i2  ]+q1.v[i1+1]*q2.v[i2+1]+q1.v[i1+2]*q2.v[i2+2];
    }
  }

}
