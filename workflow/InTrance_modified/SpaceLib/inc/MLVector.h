/*!
\file MLVector.h
\version 1.7
\brief A class for mathematical vectors and functions for vectors (Vector).
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 25.10.2006 (AO)

\b Purpose:
   This file provides the interface for the class Vector.
   Vector implements mathematical vectors.

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

#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cassert>
#include <MLGlobal.h>

namespace MathLib
{
  class Matrix;

  /** \brief A Class for mathematical vectors. */
  class Vector
  {
    friend void MultiplyQuaternion(
      const Vector& q1,
      const int     i1,
      const Vector& q2,
      const int     i2,
      Vector& r,
      const int     ir
      );

    friend void CalcRelativeQuaternion(
      const Vector& q1,
      const int     i1,
      const Vector& q2,
      const int     i2,
      Vector& dq,
      const int     idq
      );

    friend double NormalizeQuaternion(Vector& q,int first,const double eps);
    friend double Sum(const Vector& vec, const int first, const int last);
    friend Vector Stack(const Vector &, const Vector &);
    friend Vector VecPolar(double, double, double);
    friend double Dot(const Vector &, const Vector &);
    friend Vector Cross(const Vector &, const Vector &);
    friend Vector& Cross(const Vector&, const Vector&, Vector&);
    friend double Norm(const Vector &);
    friend double Sum(const Vector &);
    friend Vector operator*(double, const Vector &);
    friend Vector operator*(const Vector &, double);
    friend Vector operator/(const Vector &, double);
    friend Vector operator-(const Vector &);
    friend Vector operator+(const Vector &, const Vector &);
    friend Vector operator-(const Vector &, const Vector &);
    friend Matrix Diag(const Vector &);
    friend Vector operator*(const Matrix &, const Vector &);
    friend Vector operator*(const Vector &, const Matrix &);
    friend void MultiplyMatrixVector(const Matrix&, const Vector&, Vector&);
    friend Matrix Dyadic(const Vector &, const Vector &);
    friend std::ostream & operator<<(std::ostream &, const Vector &);

  protected:

    int      n; ///< Vector dimension.
    double * v; ///< Vector elements.

  public:

    friend class Matrix;

    /** \brief Default constructor. */
    Vector();

    /** \brief Constructs zero-initialized vector of given dimension. */
    Vector(int dim);

    /** \brief Copy constructor. */
    Vector(const Vector & vec);

    /** \brief Array constructor for vector with given dimension.
        Constructs a vector with the given dimension and initializes
        the vector with the values from the given array. It should
        be assured that the array contains size guilty double values.
        \param p    Pointer to an array of double values.
        \param dim  Dimension of the desired vector. */
    Vector(const double * p, int dim);

    /** \brief Constructor for a 2-dimensional vector. */
    Vector(double v0, double v1);

    /** \brief Constructor for a 3-dimensional vector. */
    Vector(double v0, double v1, double v2);

    /** \brief Constructor for a 4-dimensional vector. */
    Vector(double v0, double v1, double v2, double v3);

    /** \brief Constructor for a 5-dimensional vector. */
    Vector(double v0, double v1, double v2, double v3, double v4);

    /** \brief Constructor for a 6-dimensional vector. */
    Vector(double v0, double v1, double v2, double v3, double v4, double v5);

    /** \brief Constructor for a 7-dimensional vector. */
    Vector(double v0, double v1, double v2, double v3, double v4, double v5, double v6);

    /** \brief Constructor for a 8-dimensional vector. */
    Vector(double v0, double v1, double v2, double v3, double v4, double v5, double v6, double v7);

    /** \brief Constructor for a 9-dimensional vector. */
    Vector(double v0, double v1, double v2, double v3, double v4, double v5, double v6, double v7, double v8);

    /** \brief Constructor for a 10-dimensional vector. */
    Vector(double v0, double v1, double v2, double v3, double v4, double v5, double v6, double v7, double v8, double v9);

    /** Destructor. */
    inline ~Vector() {delete[] v;}

    /** \brief double assignment operator (all elements get the assigned value). */
    const Vector& operator=(const double value);

    /** \brief Vector assignment operator (equal dimensions if initialized). */
    const Vector& operator=(const Vector & vec);

    /** \brief Sets vector with argument (also with different dimensions). */
    const Vector& Set(const Vector& vec);

    /*! \brief sets vector with content of given data vector beginning from given index
     *
     * \param vec vector of new data values
     * \param idx start index in old vector
     */
    void Set(const Vector& vec, const int idx)
    {
      if (idx < n) for (int i=idx,j=0;i<n,j<vec.n;++i,++j) v[i]=vec.v[j];
    }

    /** \brief Sets all elements of vector with given value
     * \param val value to set all vector's elements with
     */
    inline void Set(const double val)
    {
      if (v!=0) for (int i = 0;i<n;++i) v[i] = val;
    }

    /** \brief Compare operator. */
    bool operator==(const Vector & vec) const;

    /** \brief Compare operator. */
    bool operator!=(const Vector & vec);

    /** \brief Component access method (Fortran notation). */
    double operator()(int i)const;

    /** \brief Component access method (Fortran notation) that returns reference. */
    double & operator()(int i);
    /** \brief Component access method (Fortran notation). */


    /// component access (copy)
    double operator[](int i) const;

    /// component access (reference)
    double& operator[](int i);


    /** \brief Component access method. */
    double getComponent(const int i) const;

    /** \brief Component modifier method.
        The modifier method sets the vector's value indicated by
        the first parameter (vector element index) to the value
        given by parameter two.
        \param i  The Index of vector element to set to new value.
        \param value The value the indicated vector element should be set to. */
    void setComponent(const int i, const double value);

    /** \brief Get sub-vector (Vector slicing).
        Returns a part of the vectors elements as new vector.
        Returns a vector that contains certain elements of the target
        vector. The desired elements are indicated by the inidices of
        the first and the last element of the source vector. The
        function does range checking.
        \param first Index of first element to extract.
        \param last  Index of last element to extract. */
    Vector Slice(int first, int last) const;

    /** \brief Vector addition. */
    void Add(const Vector & vec);

    /** \brief Vector subtraction. */
    void Subtract(const Vector & vec);

    /** \brief Scalar multiplication. */
    void Multiply(const double factor);

    /** \brief Scalar division. */
    void Divide(const double divisor);

    /** \brief Vector addition with assignment. */
    Vector & operator+=(const Vector & vec);

    /** \brief Vector subtraction with assignment. */
    Vector & operator-=(const Vector & vec);

    /** \brief Scalar multiplication with assignment. */
    Vector& operator*=(const double factor);

    /** \brief Scalar division with assignment. */
    Vector& operator/=(const double divisor);

    /** \brief Returns vector dimension. */
    int getDimension() const { return n; };

    /** \brief vector dimension (for STL like use) */
    int size(void) const {return n;};

    friend void AxisAngle2Quaternion(
      const Vector& axis,
      const double  angle,
      Vector& q
    );

    friend void AxisAngle2Quaternion(
      const Vector& axis,
      const double  angle,
      double q[4]
    );


  }; // end class Vector

  /// vector elements sum
  inline double Sum(const Vector & vec)
  {
    return Sum(vec, 0, vec.size()-1);
  }

  /// partial vector elements sum
  inline double Sum(const Vector& vec, const int first, const int last)
  {
    double sum(0.0);
    if ((vec.n>last) && (first<=last))
    {
      for (int i = first;i<=last;++i) sum+=vec[i];
    }
    return sum;
  }

    /**
  \brief Calcluates a mesaure for deviation of given quaternion derivative from
 given quaternion.
   \param q attitude quaternion
   \param q_d assumed derivative of attitude quaternion
   \return Deviation of derivative quaternion from true derivative.
  */
  inline double QuaternionDerivativeDeviation(const Vector& q
                                             ,const Vector& q_d)
  {
    return -1.0/Sum(q)*Dot(q, q_d);
  }

  inline double QuaternionDerivativeDeviation(const Vector& q, int q_first
                                             ,const Vector& q_d, int q_d_first)
  {
    const double dot(q[q_first]  *q_d[q_d_first]
                   + q[q_first+1]*q_d[q_d_first+1]
                   + q[q_first+2]*q_d[q_d_first+2]
                   + q[q_first+3]*q_d[q_d_first+3]);
    return -1.0/Sum(q, q_first, q_first+3)*dot;
  }


  /**
    \brief Correct attitude quaternion to norm of 1 (vector of doubles).
    \param q attitude quaternion
    \param eps Desired maximal deviation from unit length.
    \return Corrected attitude quaternion parameter.
  */
  inline double NormalizeQuaternion(Vector& q, const double eps = 1.0E-6)
  {
    return NormalizeQuaternion( q, 0, eps);
  }

  /// normalizes quaternion element of given data vector
  inline double NormalizeQuaternion(Vector& q,int first,const double eps=1.0E-6)
  {
    if (q.size() > (first+3)) 
      return NormalizeQuaternion(q.v[first],q.v[first+1]
      ,q.v[first+2],q.v[first+3], eps);
    else return -1.0;
  }


  /// corrects quaternion derivative element of rotational state vector
  inline double CorrectQuaternionDerivative(const Vector& q
                                    ,const int     q_first
                                    ,      Vector& q_d
                                    ,const int     q_d_first
                                    ,const double eps = 1.0E-6)
  {
    const double corr_fac(QuaternionDerivativeDeviation(q, q_first, q_d, q_d_first));
    if (fabs(corr_fac) > eps)
    {
      q_d(q_d_first)+=corr_fac;
      q_d(q_d_first+1)+=corr_fac;
      q_d(q_d_first+2)+=corr_fac;
      q_d(q_d_first+3)+=corr_fac;
      return corr_fac;
    }
    return eps;
  }
  

  /*!
   * \brief assures derivative condition for given rotational state
   *
   * \param s rotational state vector (quaternion/quaternion derivative)
   * \param eps accuracy boundary (preset to 1.0e-6)
   */
  inline void CorrectQuaternionState(Vector& s, const double eps = 1.0E-15)
  {
    switch(s.size())
    {
    case 7:
      NormalizeQuaternion(s, 0, eps);
      break;
    case 8:
      NormalizeQuaternion(s, 0, eps);    
      CorrectQuaternionDerivative(s,0,s,4,eps);
      break;
    }
  }


  /**
  \brief Corrects derivative of quaternion with respect to given quaternion (vector)
  \param q attitude quaternion
  \param q_d derivative of attitude quaternion
  \param eps if deviation of quaternion derivative from true quaternion
  derivative is greater than eps it'll be corrected, otherwise not
  \return eps, if no correction had to be applied, the correction factor otherwise
  */
  inline double CorrectQuaternionDerivative(const Vector& q
                                           ,      Vector& q_d
                                           ,const double eps = 1.0E-6)
  {
    return CorrectQuaternionDerivative(q, 0, q_d, 0, eps);
  }




  /**
  \brief Converts attitude quaternion to axis angle representation
  \param q1 quaternion component 1
  \param q2 quaternion component 2
  \param q3 quaternion component 3
  \param q4 quaternion component 4
  \param axis reference to Vector class object which is intended to hold
 the calculated rotation axis
 \param angle reference to double value which is intended to hold the calculated
 rotation angle
  */
  void Quaternion2AxisAngle(const double  q1
                          , const double  q2
                          , const double  q3
                          , const double  q4
                          ,       Vector& axis
                          ,       double& angle);

  /**
  \brief Converts axis/angle attitude representation to quaternion representation
  \param axis Vector object which holds the axis of rotation.
  \param angle The rotation Angle.
  \param q Vector object which will be filled with the four quaternion components.
  */
  void AxisAngle2Quaternion(
    const Vector& axis,
    const double  angle,
    Vector& q
  );


  /**
  \brief Converts axis/angle attitude representation to quaternion representation
  \param axis Vector object which holds the axis of rotation.
  \param angle The rotation Angle.
  \param q pointer to double array which will be filled with the four quaternion components.
  */
  void AxisAngle2Quaternion(
      const Vector& axis,
      const double  angle,
      double q[4]
  );


  void MultiplyQuaternion(
    const Vector& q1,
    const int     i1,
    const Vector& q2,
    const int     i2,
          Vector& r,
    const int     ir
  );



  /**
  \brief Calculates q1 attitude relative to q2 (expressed as quaternion)
  \param q1 attitude 1 (quaternion in vector form)
  \param i1 start index in vector q1 (where quaternion data starts)
  \param q2 attitude 2 (quaternion in vector form)
  \param i2 start index in vector q2 (where quaternion data starts)
  \param dq reference to vector that will be filled with relative attitude quaternion data
  \param i start index for modifying dq vector content
  */
  void CalcRelativeQuaternion(
    const Vector& q1,
    const int     i1,
    const Vector& q2,
    const int     i2,
    Vector& dq,
    const int     idq
    );

  /// calculates angular velocity vector from attitude quaternion and derivative
  inline void CalcAngularVelocityVec(
    const MathLib::Vector& q,    ///< attitude quaternion
    const int              iq,   ///< start index of quaternion data in q
    const MathLib::Vector& qd,   ///< attitude quaternion derivative
    const int              iqd,  ///< start index of quaternion data in qd
          MathLib::Vector& o,    ///< angular velocity vector
    const int              io    ///< start index of data in ang. vel. vector
  )
  {
   o[io  ]=2.0*(q[iq+3]*qd[iqd  ]-qd[iqd+3]*q[iq  ]+qd[iqd+1]*q[iq+2]-qd[iqd+2]*q[iq+1]);
   o[io+1]=2.0*(q[iq+3]*qd[iqd+1]-qd[iqd+3]*q[iq+1]+qd[iqd+2]*q[iq  ]-qd[iqd  ]*q[iq+2]);
   o[io+2]=2.0*(q[iq+3]*qd[iqd+2]-qd[iqd+3]*q[iq+2]+qd[iqd  ]*q[iq+1]-qd[iqd+1]*q[iq  ]);
  }

  /// calculates angular velocity vector from attitude quaternion and derivative
  inline void CalcAngularVelocityVec(
    const MathLib::Vector& q,    ///< attitude quaternion
    const MathLib::Vector& qd,   ///< attitude quaternion derivative
          MathLib::Vector& o)    ///< angular velocity vector
  {
    CalcAngularVelocityVec(q,0,qd,0,o,0);
  }

  /** \brief Vector concatenation. */
  Vector Stack(const Vector & a, const Vector & b);

  /** \brief Creates vector object from given spherical coordinates. */
  Vector VecPolar(double azim, double elev, double r = 1.0);

  /** \brief Dot product.
      Returns the dot product of the given vector arguments. */
  double Dot(const Vector & a, const Vector & b);

  /** \brief Cross product.
      Returns the cross product of the given vector arguments. */
  Vector Cross(const Vector & left, const Vector & right);

  /** \brief cross product of two vectors*/
  inline Vector& Cross(const Vector& l, const Vector& r, Vector& result) throw()
  {
#if defined development_version
    assert(l.size()==r.size() && r.size()==result.size() && result.size()==3);
#endif
    result(0) = l(1)*r(2) - l(2)*r(1);
    result(1) = l(2)*r(0) - l(0)*r(2);
    result(2) = l(0)*r(1) - l(1)*r(0);
    return result;
  } // Vector& Cross(const Vector& l,const Vector& right,Vector& result) throw()

  /** \brief Cartesian norm (modulus).
      Returns the cartesian norm (modulus) of the given vector argument. */
  double Norm(const Vector & vec);

  /** \brief Sum norm.
      Returns the sum of the argument vectors' elements. */
  double Sum(const Vector & vec);

  /** \brief Multiplication operator (scalar times vector). */
  Vector operator*(double value, const Vector & vec);

  /** \brief Multiplication operator (vector times scalar). */
  Vector operator*(const Vector & vec, double value);

  /** \brief Division operator (vector divided by double value). */
  Vector operator/(const Vector & vec, double value);

  /** \brief Negation operator (unary minus). */
  Vector operator-(const Vector & vec);

  /** \brief Vector addition operator (vector + vector). */
  Vector operator+(const Vector & a, const Vector & b);

  /** \brief Vector subtraction operator (vector - vector). */
  Vector operator-(const Vector & left, const Vector & right);

  /** \brief Vector output.
      Stream operator for output of the given vector argument's elements. */
  std::ostream & operator<<(std::ostream & os, const Vector & vec);

} // namespace Mathlib

#endif // VECTOR_H
