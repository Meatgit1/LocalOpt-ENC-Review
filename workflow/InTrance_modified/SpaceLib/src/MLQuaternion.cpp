/*!
\file MLQuaternion.cpp
\version 1.7
\brief Functions for quaternions.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 13.03.2005 (AO)

\b Purpose:
   This file provides functions for quaternions.

\b Note:
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.
   It is based on the code available at http://hiwaay.net/~crispen/vrmlworks/eai-faq/
*/

#include <math.h>

#include <MLQuaternion.h>

namespace MathLib
{
  //--------------------------------------------------------------------------
  void ConvertToQuaternion(const double axis[3], double angle, double result[4])
  {
    register double sin_a(sin(0.5*angle));
    for (int i = 0; i < 3; i++) result[i] = axis[i] * sin(0.5 * angle);
    result[3] = cos(0.5 * angle);
  }


  //--------------------------------------------------------------------------
  void ConvertToQuaternion(const Vector& axis, double angle, double result[4])
  {
    if (axis.getDimension() == 3)
    {
      register double sin_a(sin(0.5*angle));
      for (int i = 0; i < 3; i++) result[i] = axis(i) * sin_a;
      result[3] = cos(0.5 * angle);
    }
    else
    {
      result[0] = result[1] = result[2] = 0.0; result[3] = 1.0;
    }
  }


//--------------------------------------------------------------------------
  void ConvertFromQuaternion(double q[4], double axis[3], double* angle)
  {
    *angle = 2.0 * acos(q[3]);

    for (int i = 0; i < 3; i++)
    {
      if (*angle != 0.0)
        axis[i] = q[i] / sin(*angle / 2.0);
      else
        axis[i] = 0.0;
    }
  }


  //--------------------------------------------------------------------------
  void MultiplyQuaternion(
    const double q1[4]
   ,const double q2[4]
   ,double result[4])
  {
    // correct code for quaternion multiplication (used in member function of
    // CAttitudeQuaternion)
    //result[0] = q2[3] * q1[0] + q2[0] * q1[3] + q1[1] * q2[2] - q1[2] * q2[1];
    //result[1] = q2[3] * q1[1] + q2[1] * q1[3] + q1[2] * q2[0] - q1[0] * q2[2];
    //result[2] = q2[3] * q1[2] + q2[2] * q1[3] + q1[0] * q2[1] - q1[1] * q2[0];
    //result[3] = q2[3] * q1[3] - q2[0] * q1[0] - q2[1] * q1[1] - q2[2] * q1[2];

    // Wrong code for quaternion multiplication though it involves an incorrect
    // formula for calculating the vector part of resulting quaternion.
    // Nevertheless the wrong code has to be used because in 
    // PrintVRML(const char*) the wrong behaviour is needed because the rest of
    // this code corrects it implicitely. So, until no major code workthough
    // is done in the mentioned function, it's useless to correct the formula
    // only here.
    result[0] = q2[3] * q1[0] + q2[0] * q1[3] - q1[1] * q2[2] + q1[2] * q2[1];
    result[1] = q2[3] * q1[1] + q2[1] * q1[3] - q1[2] * q2[0] + q1[0] * q2[2];
    result[2] = q2[3] * q1[2] + q2[2] * q1[3] - q1[0] * q2[1] + q1[1] * q2[0];
    result[3] = q2[3] * q1[3] - q2[0] * q1[0] - q2[1] * q1[1] - q2[2] * q1[2];
  }


  //--------------------------------------------------------------------------
  Matrix Quaternion2Matrix(const Vector& v)
  {

    Vector v_v(v.Slice(0,2)); // get vector component of quaternion

    Matrix R(3,3);
    R(0,0) = v[0]*v[0]-v[1]*v[1]-v[2]*v[2]+v[3]*v[3];
    R(0,1) =2*(v[0]*v[1]+v[2]*v[3]);
    R(0,2) =2*(v[0]*v[2]-v[1]*v[3]);
    R(1,0) =2*(v[0]*v[1]-v[2]*v[3]);
    R(1,1) = -v[0]*v[0]+v[1]*v[1]-v[2]*v[2]+v[3]*v[3];
    R(1,2) =2*(v[1]*v[2]+v[0]*v[3]);
    R(2,0) =2*(v[0]*v[2]+v[1]*v[3]);
    R(2,1) =2*(v[1]*v[2]-v[0]*v[3]);
    R(2,2) = -v[0]*v[0]-v[1]*v[1]+v[2]*v[2]+v[3]*v[3];
//    std::cout <<"R:\n"<<R<<std::endl;
    // cross product of quaternion's vector component as matrix
  
//    Matrix Q(CrossMatrix(v_v));
/*    Matrix R2(3,3);
    R2 = (v(3)*v(3) - Dot(v_v,v_v) )*Id(3) + 2.0*Dyadic(v_v,v_v) - 2.0*v(3)*Q;
    std::cout <<"R2:\n"<<R2<<std::endl;

*/
    return R;//(v(3)*v(3) - Dot(v_v,v_v) )*Id(3) + 2.0*Dyadic(v_v,v_v) - 2.0*v(3)*Q;
  }


  //---------------------------------------------------------------------------
  CAttitudeQuaternion operator*(const double               lhs
                               ,const CAttitudeQuaternion& rhs)
  {
    return CAttitudeQuaternion(Vector(
                              lhs*rhs.m_dat[0]
                             ,lhs*rhs.m_dat[1]
                             ,lhs*rhs.m_dat[2]
                             ,lhs*rhs.m_dat[3]));
  }





} // end of namespace MathLib
