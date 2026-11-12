// file     : CAttitudeData.h
// author   : Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
// compiler : GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++ 
// date     : 19.05.2005 (AO)
// purpose  : Definition of the data structure of CAttitude class.

#ifndef CATTITUDEDATA_H
#define CATTITUDEDATA_H

// global functions of SpaceLib
#include <Global.h>

// MathLib
#include <MLVector.h>
#include <MLMatrix.h>
#include <MLQuaternion.h>

#define NEW_CODE 1
/*!
\todo AO:remove old code and conditional compiling of CAttitudeData (#define NEW_CODE)
*/

namespace SpaceLib 
{

  /// CAttitude internal data structure 
  struct CAttitudeData 
  {

    /// attitude matrix representation (or matrix of direction cosinus)
    MathLib::Matrix m_matrix;

#if defined NEW_CODE
    /// rotation quaternion representation
    MathLib::CAttitudeQuaternion m_q;
#else
    double m_q[4];
#endif

    // axis/angle representation
    MathLib::Vector m_axis;
    double m_angle;


    /// helper method
    void updateFromMatrix(void)
    {
      Matrix2Quaternion(m_matrix, m_q);
      Matrix2AxisAngle(m_matrix, m_axis, m_angle);
    }

    /// helper method
    void updateFromQuaternion(void)
    {
      m_matrix = Quaternion2Matrix(m_q);
      Quaternion2AxisAngle(m_q, m_axis, m_angle);
    }

    /// helper method
    void updateFromAxisAngle(void)
    {
      AxisAngle2Matrix(m_axis, m_angle, m_matrix);
      AxisAngle2Quaternion(m_axis, m_angle, m_q);
    }


    /// standard constructor
    CAttitudeData(void)
      :m_matrix(MathLib::Id(3))
#if defined NEW_CODE
      ,m_q()
#endif
      ,m_axis(3)
      ,m_angle(0.0) 
    {
#ifndef NEW_CODE
      m_q[0] = m_q[1] = m_q[2] = 0.0; m_q[3] = 1.0;
#endif
    }

#if defined NEW_CODE
    /// attitude quaternion constructor
    CAttitudeData(const MathLib::CAttitudeQuaternion& q)
      :m_matrix(MathLib::Id(3))
      ,m_q(q)
      ,m_axis(3)
      ,m_angle(0.0)
    {}
#endif


    /// matrix constructor
    CAttitudeData(const MathLib::Matrix& matrix):
      m_matrix(matrix)
#if defined NEW_CODE
      ,m_q()
#endif
      ,m_axis(3)
      ,m_angle(0.0)
    {
#ifndef NEW_CODE        
       q[0] = q[1] = q[2] = 0.0; q[3] = 1.0;
#endif
    }


    /// quaternion constructor
    CAttitudeData(const double q1
                , const double q2
                , const double q3
                , const double q4)
    : m_matrix(MathLib::Id(3))
#if defined NEW_CODE
    , m_q(q1, q2, q3, q4)
#endif
    , m_axis(3)
    , m_angle(0.0)
    {
#ifndef NEW_CODE
       m_q[0] = m_q1;q[1] = m_q2; m_q[2] = m_q3; m_q[3] = q4;
#endif
    }


    /// axis/angle constructor
    CAttitudeData(const MathLib::Vector& _axis, const double _angle)
    : m_matrix(MathLib::Id(3))
#if defined NEW_CODE
    , m_q()
#endif
    , m_axis(_axis)
    , m_angle(_angle) 
    {
#ifndef NEW_CODE
      m_q[0] = m_q[1] = m_q[2] = 0.0; m_q[3] = 1.0;
#endif
    }


    /// unit vectors constructor
    CAttitudeData(
       const MathLib::Vector& e1
     , const MathLib::Vector& e2
     , const MathLib::Vector& e3)
    : m_matrix(MathLib::Id(3))
#if defined NEW_CODE
    , m_q()
#endif
    , m_axis(3)
    , m_angle(0.0)
    {
      if ((e1.getDimension() == 3) && (e2.getDimension() == 3) && (e3.getDimension() == 3)) {
        m_matrix.SetRow(0, e1);
        m_matrix.SetRow(1, e2);
        m_matrix.SetRow(2, e3);
      } else FatalError("CAttitudeData::CAttitudeData(const Vector&, const Vector&, const Vector&) => Wrong index of unit vector(s)");
#ifndef NEW_CODE
      m_q[0] = m_q[1] = m_q[2] = 0.0; m_q[3] = 1.0;
#endif

    }



    /// copy constructor
    CAttitudeData(const struct CAttitudeData& to_copy) 
      : m_matrix(to_copy.m_matrix)
#if defined NEW_CODE
      , m_q(to_copy.m_q)
#endif
      , m_axis(to_copy.m_axis)
      , m_angle(to_copy.m_angle)
    {
#ifndef NEW_CODE
       m_q[0] = to_copy.m_q[0];q[1] = to_copy.m_q[1];
       m_q[2] = to_copy.m_q[2];q[3] = to_copy.m_q[3];
#endif
    };



    /// assign operator
    const struct CAttitudeData& operator=(const struct CAttitudeData& to_copy) 
    {
      m_matrix = to_copy.m_matrix;
#if defined NEW_CODE
      m_q = to_copy.m_q;
#endif
      m_axis = to_copy.m_axis;
      m_angle = to_copy.m_angle;

#ifndef NEW_CODE
      m_q[0] = to_copy.m_q[0];q[1] = to_copy.m_q[1];
      m_q[2] = to_copy.m_q[2];q[3] = to_copy.m_q[3];
#endif

      return *this;
    };

  };

} // namespace SpaceLib

#endif
