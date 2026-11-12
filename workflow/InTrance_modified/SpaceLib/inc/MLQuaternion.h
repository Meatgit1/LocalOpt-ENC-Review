/*!
\file MLQuaternion.h
\version 1.7
\brief Functions for quaternions.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 11.05.2005 (AO)

\b Purpose:
   This file provides functions for quaternions.

\b Note:
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.
   It is based on the code available at http://hiwaay.net/~crispen/vrmlworks/eai-faq/
*/

#ifndef QUATERNION_H
#define QUATERNION_H

#include <math.h>

// STL
#include <ostream>
#include <iostream>

// MathLib
#include <MLVector.h>
#include <MLMatrix.h>

namespace MathLib
{
  /** \brief Converts axis/angle representation to quaternion form.
      Assumes axis is unit-length, angle is in radians.
      \param axis axis unit vector (3-dim array of doubles)
      \param angle rotation angle [rad]
      \param result resulting quaternion (pointer to a 4-dim array of doubles) */
  void ConvertToQuaternion(const double axis[3], double angle, double result[4]);


  /** \brief Converts axis/angle representation to quaternion form.
      Assumes axis is unit-length, angle is in radians.
      \param axis axis unit vector (Vector class object)
      \param angle rotation angle [rad]
      \param result resulting quaternion (pointer to a 4-dim array of doubles) */
  void ConvertToQuaternion(const Vector& axis, double angle, double result[4]);


  /** \brief Converts unit-length quaternion form to axis/angle representation.
      \param q quaternion that is to be converted (4-dim array of doubles)
      \param axis resulting axis vector (3-dim array of doubles)
      \param angle angle [rad] (pointer to a double) */
  void ConvertFromQuaternion(double q[4], double axis[3], double * angle);

  /** \brief Multiplication of two quaternions.
      \param q1 left quaternion (4-dim array of doubles)
      \param q2 right quaternion (4-dim array of doubles)
      \param result resulting quaternion (4-dim array of doubles) */
  void MultiplyQuaternion(const double q1[4], const double q2[4], double result[4]);

  class CAttitudeQuaternion;
  /** \brief Converts quaternion representation to matrix form.

      \param q resulting quaternion (pointer to a 4-dim array of doubles)
      \return AttitudeMatrix corresponding to given attitude quaternion.
  */

  Matrix Quaternion2Matrix(const Vector&);

  inline Matrix Quaternion2Matrix(double q[4]) 
  {
    return Quaternion2Matrix(Vector(q,4));
  }
/*
  inline Matrix Quaternion2Matrix(const CAttitudeQuaternion& q)
  {
    return Quaternion2Matrix(q.getDataVector());
  }
*/
  inline Matrix Quaternion2Matrix(const double q1, const double q2
                                , const double q3, const double q4)
  {
    return Quaternion2Matrix(Vector(q1, q2, q3, q4));
  }


  /*!
    \brief Converts quaternion to axis angle representation

  \param q1 component 1 of quaternion to convert
  \param q2 component 2 of quaternion to convert
  \param q3 component 3 of quaternion to convert  
  \param q4 component 4 of quaternion to convert
  \param axis return value of rotation axis vector
  \return rotation angle
  */
  double Quaternion2AxisAngle(const double q1
                            , const double q2
                            , const double q3
                            , const double q4
                            , Vector& axis);


  /*!
    \brief Converts quaternion to axis angle representation

  \param q quaternion to convert (pointer to 4-dimensional array)
  \param axis return value of rotation axis vector
  \return rotation angle
  */
  inline void Quaternion2AxisAngle(const double q[4], Vector& axis, double& angle)
  {
    Quaternion2AxisAngle(q[0], q[1], q[2], q[3], axis, angle);
  }



  /*!
  \brief Expression of attitudes in 3D-space with quaternions.

  */
  class CAttitudeQuaternion
  {

  private:

    double m_dat[4];


    /// Sets internal data to represent a unit quaternion (no rotation).
    inline void initialize(void)
    {
        for (int i=0; i<3; ++i) m_dat[i] = 0.0; m_dat[3]=1.0;
    };


  public:

    /*!
    \brief Possible input types for quaternion initialization with vector/value
 type parameters

      The vector/value - type constructor takes a vector argument and a value
      argument. A quaternion may be initialized with a rotation axis and
      rotation angle or a vector component in combination with a value
      component. So, if you use the vector/value - constructor of this class,
      with the third parameter (be aware of default parameters) you are able to
      determine the type of data you are using for initialization.
    */
    enum QUAT_INPUT_TYPE
    {
       VECTOR_VALUE = 0
      ,AXIS_ANGLE   = 1
    };


      /// standard constructor
      CAttitudeQuaternion(void) 
      { 
        initialize();
      }



      /// Axis/angle constructor
      CAttitudeQuaternion(const double axis[3], const double angle)
      {
        ConvertToQuaternion(axis, angle, m_dat);
      };



      /// Vector/Value constructor.
      CAttitudeQuaternion(
        const Vector& vec
       ,const double val
       ,const QUAT_INPUT_TYPE type = AXIS_ANGLE)
      {
        switch (type)
        {
          case VECTOR_VALUE:
            if (vec.getDimension() == 3)
            {
              for (int i=0;i<3;++i) m_dat[i] = vec(i);
              m_dat[3] = val;
            }
          break;

          case AXIS_ANGLE:
            AxisAngle2Quaternion(vec, val, m_dat);
          break;
          default:
            {
              std::cout << "CAttitudeQuaternion(const Vector&, const double):"
                        << "Trying to initialize quaternion with " 
                        << vec.getDimension() << "-dimensional"
                        << " vector! It will be initialized with unit"
                        << "quaternion.\n";
 
              initialize();          // init with default values
            }
        }
      } // CAttitudeQuaternion



      /// Quaternion elements constructor.
      CAttitudeQuaternion(const double v1
                        ,const double v2
                        ,const double v3
                        ,const double v4) 
      {
        m_dat[0] = v1;
        m_dat[1] = v2;
        m_dat[2] = v3;
        m_dat[3] = v4;
      };



      /// Quaternion elements constructor (vector of quat. elements parameter)
      CAttitudeQuaternion(const Vector& data)
      {
        m_dat[0] = data(0);
        m_dat[1] = data(1);
        m_dat[2] = data(2);
        m_dat[3] = data(3);
      };

      

      /// Vector part of quaternion (first three elements of quaternion)
      Vector getVector(void) const
      {
        return Vector(m_dat[0], m_dat[1], m_dat[2]);
      }


      /// Vector filled with all four quaternion elements
      Vector getDataVector(void) const
      {
        return Vector(m_dat[0], m_dat[1], m_dat[2], m_dat[3]);
      }


      /// Replaces internal four quaternion components with given values
      /// but only if argument vector has a dimension of four.
      void setDataVector(const Vector& vec)
      {
        if (vec.getDimension() == 4)
        {
          m_dat[0] = vec(0);m_dat[1] = vec(1);
          m_dat[2] = vec(2);m_dat[3] = vec(3);
        }
      }


      /// Value part of quaternion
      double getValue(void) const
      {
        return m_dat[3];
      }


      /// Sets value part of quaternion.
      void setValue(const double val)
      {
        m_dat[3] = val;
      }


      /// Sets quaternion components with vector/value arguments depending on
      /// given type of inputparameters.
      void set(const Vector& vec
             , const double val
             , const QUAT_INPUT_TYPE type = AXIS_ANGLE)
      {
        switch (type)
        {
         case VECTOR_VALUE:
          if (vec.getDimension() == 3)
          {
            m_dat[0] = vec(0);
            m_dat[1] = vec(1);
            m_dat[2] = vec(2);
            m_dat[3] = val;
          }
          break;
         case AXIS_ANGLE:
           ConvertToQuaternion(vec, val, m_dat);
          break;
        default:
          break;
        } // switch (type)
      }


      /// Return inverse of quaternion
      CAttitudeQuaternion inverse(void) const
      {
        return CAttitudeQuaternion(-m_dat[0], -m_dat[1], -m_dat[2], m_dat[3]);
      }


      /// Invert quaternion
      void invert(void)
      {
        for (int i = 0;i<3;++i) m_dat[i] = - m_dat[i];
      }


      /// Normalizes quaternion
      void normalize(void)
      {
        const double norm(sqrt(m_dat[0]*m_dat[0]+
                               m_dat[1]*m_dat[1]+
                               m_dat[2]*m_dat[2]+
                               m_dat[3]*m_dat[3]));

        for (int i=0; i<4; ++i) m_dat[i]/=norm;
      }


      /// Multiplication operator
      inline CAttitudeQuaternion operator*(const CAttitudeQuaternion& rhs) const
      {
        return CAttitudeQuaternion(
          rhs[3] * m_dat[0] + rhs[0] * m_dat[3] + m_dat[1] * rhs[2] - m_dat[2] * rhs[1]
         ,rhs[3] * m_dat[1] + rhs[1] * m_dat[3] + m_dat[2] * rhs[0] - m_dat[0] * rhs[2]
         ,rhs[3] * m_dat[2] + rhs[2] * m_dat[3] + m_dat[0] * rhs[1] - m_dat[1] * rhs[0]
         ,rhs[3] * m_dat[3] - rhs[0] * m_dat[0] - rhs[1] * m_dat[1] - rhs[2] * m_dat[2]
        );
      }

      /// multiplication assign operator
      inline const CAttitudeQuaternion& operator*=(const CAttitudeQuaternion& rhs)
      {
        double t[4];
        t[0] = rhs[3] * m_dat[0] + rhs[0] * m_dat[3] + m_dat[1] * rhs[2] - m_dat[2] * rhs[1];
        t[1] = rhs[3] * m_dat[1] + rhs[1] * m_dat[3] + m_dat[2] * rhs[0] - m_dat[0] * rhs[2];
        t[2] = rhs[3] * m_dat[2] + rhs[2] * m_dat[3] + m_dat[0] * rhs[1] - m_dat[1] * rhs[0];
        t[3] = rhs[3] * m_dat[3] - rhs[0] * m_dat[0] - rhs[1] * m_dat[1] - rhs[2] * m_dat[2];
        m_dat[0]=t[0];
        m_dat[1]=t[1];
        m_dat[2]=t[2];
        m_dat[3]=t[3];
        return (*this);
      }


      /// Destructor
      ~CAttitudeQuaternion(void) {};


      /// access quaternion elements by index (0,..,3)
      double operator() (const unsigned int index) const
      {
        if ( index < 4) return m_dat[index];
        std::cout << "double CAttitudeQuaternion::operator() => \
index not valid!\n";
        return 0;
      }

      /// reference access quaternion elements by index (0,..,3)
      double& operator() (const unsigned int index)
      {
        if ( index < 4) return m_dat[index];
        std::cout << "double& CAttitudeQuaternion::operator() => \
index not valid (returned reference to first element!\n";
        return m_dat[0];
      }

      /// fortran access quaternion elements by index (0,..,3)
      //
      double operator[] (const unsigned int index) const
      {
        if ( index < 4) return m_dat[index];
        else 
        {
          std::cout << "double CAttitudeQuaternion::operator[]: index out of\
 bound (0..3)!\n";
          return 0.0;
        }
      }

      /// fortran reference access quaternion elements by index (0,..,3)
      double& operator[] (const unsigned int index)
      {
        if ( index < 4) return m_dat[index];
        std::cout << "double& CAttitudeQuaternion::operator[] => \
index not valid (returned reference to first element!\n";
        return m_dat[0];
      }


      /// Multiplication/assign ooperator
      CAttitudeQuaternion& operator*=(const double rhs)
      {
        for (int i = 0; i<4;++i) m_dat[i]*=rhs;
        return *this;
      }


      /// Scalar multiplication operator
      CAttitudeQuaternion operator*(const double rhs)
      {
        CAttitudeQuaternion t(*this);
        t*= rhs;
        return t;
      }

      /// minus assign operator
      const CAttitudeQuaternion& operator-=(const CAttitudeQuaternion& rhs)
      {
        m_dat[0]-=rhs.m_dat[0];
        m_dat[1]-=rhs.m_dat[1];
        m_dat[2]-=rhs.m_dat[2];
        m_dat[3]-=rhs.m_dat[3];
        return (*this);
      }

      /** \brief Ángular velocity vector from quaternion and it's derivative .
      \param q_d Derivative of quaternion
      \return Angular velocity vector with in quaternion frame coordinates
      */
      inline Vector AngularVelocityVector(const CAttitudeQuaternion& q_d) const
      {
        return 2.0*(inverse()*q_d).getVector();
      }


      /** \brief first derivative of attitude quaternion due to angular velocity*/
      inline CAttitudeQuaternion calcFirstDerivative(const Vector& omega) const
      {
        return 0.5 * (*this * CAttitudeQuaternion(omega, 0.0,VECTOR_VALUE));
      }


      friend double NormalizeQuaternion(CAttitudeQuaternion& q,const double eps);
      friend CAttitudeQuaternion operator*(const double lhs, const CAttitudeQuaternion& rhs);

      friend inline void Matrix2Quaternion(const Matrix& M, CAttitudeQuaternion& q);

      friend inline void AxisAngle2Quaternion(
         const Vector&              axis
       , const double               angle
       ,       CAttitudeQuaternion& q);

  };


  

  /**
  \brief Calcluates a mesaure for deviation of given quaternion derivative from
 given quaternion.
   \param q attitude quaternion
   \param q_d assumed derivative of attitude quaternion
   \return Deviation of derivative quaternion from true derivative.
  */
  inline double QuaternionDerivativeDeviation(const CAttitudeQuaternion& q
                                             ,const CAttitudeQuaternion& q_d)
  {
    return -1.0/Sum(q.getDataVector())
      *Dot(q.getDataVector(), q_d.getDataVector());
  }


  /**
  \brief Corrects derivative of quaternion with respect to given quaternion
 (CAttitudeQuaternion)
  \param q attitude quaternion
  \param q_d derivative of attitude quaternion
  \param eps if deviation of quaternion derivative from true quaternion
  derivative is greater than eps it'll be corrected, otherwise not
  \return eps, if no correction had to be applied, the correction factor
 otherwise
  */
  inline double CorrectQuaternionDerivative(const CAttitudeQuaternion  q
                                           ,      CAttitudeQuaternion& q_d
                                           ,const double               eps = 1.0E-6)
  {
    const double corr_fac(QuaternionDerivativeDeviation(q, q_d));
    if (fabs(corr_fac) > eps)
    {
      q_d(0)+=corr_fac;
      q_d(1)+=corr_fac;
      q_d(2)+=corr_fac;
      q_d(3)+=corr_fac;
      return corr_fac;
    }
    return eps;
  }



  /**
  \brief Converts given derivative of also given quaternion to angular velocity
  vector.
  \param q attitude quaternion
  \param q_d first derivative of quaternion
  \return Angular velocity vector in coordinates of quaternion frame.
  */
  inline Vector QuaternionDerivative2AngularVelocityVector(
    const CAttitudeQuaternion& q
  , const CAttitudeQuaternion& q_d)
  {
    return q.AngularVelocityVector(q_d);
  }



  /**
    \brief Correct attitude quaternion to norm of 1 (CAttitudeQuaternion).
    \param q attitude quaternion
    \param eps Desired maximal deviation from unit length.
    \return Corrected attitude quaternion parameter.
  */
  inline double NormalizeQuaternion(CAttitudeQuaternion& q,const double eps = 1.0E-6)
  {
    return NormalizeQuaternion(q.m_dat[0], q.m_dat[1], q.m_dat[2], q.m_dat[3], eps);
  }

  /**
    \brief Correct attitude quaternion to norm of 1 (Vector).
    \param q attitude quaternion
    \param eps Desired maximal deviation from unit length.
    \return Corrected attitude quaternion parameter.
  */
/*  double CorrectQuaternion(Vector& q, const double eps = 1.0E-6)
  {
    if (q.getDimension() == 4)
    {
      CAttitudeQuaternion q_q(q);
      const double ret_val(CorrectQuaternion(q_q, eps));
      q = q_q.getDataVector();
      return ret_val;
    }
    return -1;
  }
*/

  /*!
    \brief Converts given matrix to quaternion 

    \param M rotation matrix
    \param q quaternion 
    \return modified parameter q
  */
  inline void Matrix2Quaternion(const Matrix& M, CAttitudeQuaternion& q)
  {
    if (M.isSquare() == true && M.getNumberOfColumns() == 3)
    {
      M.calcQuaternion(q.m_dat);
    }
    else
    {
      std::cout << "calcQuaternion(const Matrix& M, CAttitudeQuaternion& q):\
Given rotation marix is no square matrix or has dim > 3!\n";
    }
  }

  /**
  */
  inline void AxisAngle2Quaternion(
         const Vector&              axis
       , const double               angle
       ,       CAttitudeQuaternion& q)
  {
    AxisAngle2Quaternion(axis, angle, q.m_dat);
  }


  /*!
    \brief Converts quaternion to axis angle representation

  \param q quaternion to convert
  \param axis return value of rotation axis vector
  \return rotation angle
  */
  inline void Quaternion2AxisAngle(const CAttitudeQuaternion& q
                                  ,      Vector&              axis
                                  ,      double&              angle)
  {
    Quaternion2AxisAngle(q[0],q[1],q[2],q[3], axis, angle);
  };


  /// helper function to allow "Invert(q)"
  inline void Invert(CAttitudeQuaternion& q)
  {
    q.invert();
  }



  /// helper function to allow " AttitudeQuaternion q2 = Inverse(q1);"
  inline CAttitudeQuaternion Inverse(const CAttitudeQuaternion& q) 
  {
    CAttitudeQuaternion q2(q);
    q2.invert();
    return q2;
  }

  /// helper function to allow "Normalize(q)"
  inline void Normalize(CAttitudeQuaternion& q)
  {
    q.normalize();
  }


  double Quaternion2AxisAngle(const Vector& q, Vector& axis);

  double Quaternion2AxisAngle(const CAttitudeQuaternion& q, Vector& axis);


  inline Matrix Quaternion2Matrix(const CAttitudeQuaternion& q)
  {
    return Quaternion2Matrix(q.getDataVector());
  }
} // end of namespace MathLib

#endif // QUATERNION_H
