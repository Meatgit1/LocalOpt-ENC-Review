/*!
\file MLMatrix.h
\version 1.7
\brief A class for mathematical matrices and functions for matrices (Matrix).
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 21.09.2006 (AO)

\b Purpose:
   Like the class vector the class declarated in this
   file (Matrix) is heavily used in the library
   SpaceLib to perform different transformations of
   vectors (distance and velocity). MathLib's member
   class Matrix provides some useful functions for
   this. Of course the basic Matrix operations are
   possible too.

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
   (c) 1999-2005 Oliver Montenbruck, Eberhard Gill, Andreas Ohndorf, Bernd Dachwald
*/

#ifndef MATRIX_H
#define MATRIX_H

// STL
#include <iostream>


// MathLib
#include <MLVector.h>


namespace MathLib 
{

  class Vector;

  /** \brief A class for mathematical matrices and functions for matrices.
      The class Matrix provides an encapsulation of data representation and
      mathematical operations for matrices. */
  class Matrix
  {
    friend Matrix CrossMatrix(const Vector&);
    friend Matrix Id(int);
    friend Matrix Diag(const Vector &);
    friend Matrix R_x(double);
    friend void R_x(double, Matrix&);
    friend Matrix R_y(double);
    friend void R_y(double, Matrix&);
    friend Matrix R_z(double);
    friend void R_z(double, Matrix&);
    friend Matrix Transp(const Matrix &);
    friend Matrix Inv(const Matrix &);
    friend Matrix operator*(double, const Matrix &);
    friend Matrix operator*(const Matrix& , double);
    friend Matrix operator/(const Matrix &, double);
    friend Matrix operator-(const Matrix &);
    friend Matrix operator+(const Matrix &, const Matrix &);
    friend Matrix operator-(const Matrix &, const Matrix &);
    friend Matrix operator*(const Matrix &, const Matrix &);
    friend void MultiplyMatrices(const Matrix&, const Matrix&, Matrix&);
    friend void MultiplyMatrixVector(const Matrix&, const Vector&, Vector&);
    friend Vector operator*(const Matrix &, const Vector &);
    friend Vector operator*(const Vector &, const Matrix &);

    friend Matrix Dyadic(const Vector&, const Vector &);
    friend std::ostream &  operator<<(std::ostream &, const Matrix &);
  
  protected:

    int       n; ///< Number of rows.
    int       m; ///< Number of columns.
    double ** M; ///< Matrix elements.

  public:

    /** \brief Default constructor.
        Constructs a Matrix object of the dimension (0 x 0). */
    Matrix(void);

    /** \brief Constructor for a (dim_row x dim_col) matrix.
        Constructs a matrix with the dimension (dim_row x dim_col).
        All elements are initialized with the value 0.
        \param dim_row number of rows
        \param dim_col number of columns */
    Matrix(int dim_row, int dim_col);

    /** \brief Copy constructor. */
    Matrix(const Matrix & mat);

    /** \brief Array copy constructor.

        This constructor initializes a matrix with
        data values from a given array of doubles. The
        argument pointer to a double array should consist of
        the matrix row vectors. That means the first
        n_col double values should be the ones of the
        first row of the matrix to construct.
        \param array pointer to double data array
        \param n_row number of matrix rows
        \param n_col number of matrix columns
    */
    Matrix(const double *array, int n_row, int n_col);

    /** \brief Destructor. */
    ~Matrix();

    /** \brief Assignment operator (double).
        Assigns the given double value to all elements of the matrix. */
    Matrix & operator=(const double value);

    /** \brief Assignment operator (Matrix).
        Copies the content of the matrix given as argument into the matrix data fields. */
    Matrix & operator=(const Matrix & mat);

    /** \brief Returns the number of the rows. */
    int getNumberOfRows() const
    {
      return n;
    };

    /** \brief Returns the number of columns. */
    int getNumberOfColumns() const
    {
      return m;
    };

    /** \brief Component access (Fortran notation). */
    double operator()(int row, int col)const
    {
      return M[row][col];
    }

    /** \brief Component access.
        \param row row
        \param col column */
    double &operator()(int row, int col)
    {
      return M[row][col];
    }

    /** \brief Returns a matrix column as vector.
        \param col column */
    Vector Col(int col) const;

    /** \brief Returns a matrix row as vector.
        \param i row */
    Vector Row(int i) const;

    /** \brief Returns the diagonal of a quadratic matrix as vector. */
    Vector Diag(void) const;

    /** \brief Returns the trace of a quadratic matrix.
        The trace of a quadratic matrix is the sum of all elements along the matrix diagonal. */
    double Trace(void) const;

    /** \brief True if it is a square matrix.*/
    bool isSquare(void) const { return m==n;};


    /** \brief Converts rotation matrix to quaternion (double values)
     If the matrix is of dimension (3x3) and a rotation transformation, this
     transformation may also be described by a quaternion. Because a quaternion
     can be thought of as combination of a vector and a single value, the 
     returned quaternion values will be of the form vector/value with 
     the following order:\n
     q=(q1 = vec1, q2 = vec2, q3 = vec3, q4 = val)\n

    \param q Vector that will be filled with the 4 quaternion components(vec/val).
    \return modified parameter q
    */
    void calcQuaternion(double& q1, double& q2, double& q3, double& q4) const;
 
    
    /** \brief Converts rotation matrix to quaternion(double array)
     If the matrix is of dimension (3x3) and a rotation transformation, this
     transformation may also be described by a rotation quaternion. Because a
     quaternion can be thought of as combination of a vector and a single value,
     the returned quaternion values will be of the form vector/value with 
     the following order:\n
     q=(q[0] = vec1, q[1] = vec2, q[2] = vec3, q[3] = val)\n
    \param q double array where the resulting quaternion components will be
 stored in
    \return Updated parameter q
    */
    void calcQuaternion(double q[4]) const
    {
      calcQuaternion(q[0], q[1], q[2], q[3]);
    }



    /** \brief Converts rotation matrix to quaternion(vector of double)
     If the matrix is of dimension (3x3) and a rotation transformation, this
     transformation may also be described by a rotation quaternion. Because a
     quaternion can be thought of as combination of a vector and a single value,
     the returned quaternion values will be of the form vector/value with 
     the following order:\n
     q=(q(0) = vec1, q(1) = vec2, q(2) = vec3, q(3) = val)\n
    \param q double array where the resulting quaternion components will be
 stored in
    \return Updated parameter q
    */
    void calcQuaternion(Vector& q) const;



    /** \brief Converts rotation matrix to axis/angle (in refrences to double)
     If the matrix is of dimension (3x3) and a rotation transformation, this
     transformation may also be described by a rotation about a certain axis
     with a certain angle. The returned axis/angle values are organized as
     follows:\n
     axis =(v1 = vec_1, v2 = vec_2, v3 = vec_3, angle = val)\n
    \param v_1 axis vector component 1
    \param v_2 axis vector component 2
    \param v_3 axis vector component 3
    \param angle rotation angle
    \return Updated parameters
    */
    void calcAxisAngle(double& v1, double& v2, double& v3, double& angle) const;


    /** \brief Converts rotation matrix to axis/angle (double array, double)
     If the matrix is of dimension (3x3) and a rotation transformation, this
     transformation may also be described by a rotation about a certain axis
     with a certain angle. The returned axis/angle values are organized as
     follows:\n
     axis =(axis[0] = vec_1, axis[1] = vec_2, axis[3] = vec_3, angle = val)\n
    \param axis rotation axis
    \param angle rotation angle
    \return Updated parameters
    */
    void calcAxisAngle(double axis[3], double& angle) const;


    /** \brief Converts rotation matrix to axis/angle (vector, double)
     If the matrix is of dimension (3x3) and a rotation transformation, this
     transformation may also be described by a rotation about a certain axis
     with a certain angle. The returned axis/angle values are organized as
     follows:\n
     axis =(vec[0] = vec_1, vec[1] = vec_2, vec[3] = vec_3, angle = val)\n
    \param vec rotation axis
    \param angle rotation angle
    \return Updated parameters
    */
    void calcAxisAngle( Vector& vec, double& angle) const;


    /** \brief Returns a submatrix indicated by row and column indices.
        The returned matrix is the submatrix whose upper left
        corner is indicated by the first_*-indices and the
        lower right corner by the last_*-indices. */
    Matrix slice(int first_row, int last_row, int first_col, int last_col);

    /** \brief Sets the desired column with a vector's content.
        The elements of the column indicated by the argument
        col will be replaced by the vectors elements.
        \param col column index in matrix of the column to change
        \param colvec vector containing the new values of the column to change */
    void SetCol(int col, const Vector & colvec);

    /** \brief Sets the desired row with a vector's content.
        The elements of the row indicated by the argument
        row will be replaced by the vectors elements.
        \param row row index in matrix of the column to change
        \param rowvec vector containing the new values of the row to change */
   void SetRow(int row, const Vector & rowvec);

    /** \brief Sets the elements of desired submatrix.
        The elements of the matrix's submatrix, which upper left
        corner is indicated by the given row and column index will
        be set with the elements of the given data matrix. Thus the
        element of the upper left corner will be set with the content
        of the upper left corner's element of the given data matrix.
        If there aren't enough rows and/or columns in the matrix
        to set, to transfer all the data matrix' content, no element
        will be altered.
        \param row row of upper left corner of submatrix to set
        \param col column of upper left corner of submatrix to set
        \param mat matrix with the new elements of the submatrix */
    void SetSubMatrix(int row, int col, const Matrix & mat);

  }; // end of class Matrix







  /*!\brief Rotation matrix -> roational quaternion (array of double)
  \param matrix 3x3 rotation transformation matrix
  \param q array of at least 4 double values to store quaternion elements
  \return updated array of double values (parameter q)
  */
  inline void Matrix2Quaternion(const Matrix& matrix, double q[4])
  {
    matrix.calcQuaternion(q);
  }

  /*!\brief Rotation matrix -> roational quaternion (4x1 vector)
  \param matrix 3x3 rotation transformation matrix
  \param v Vector to store quaternion elements (will be resized if necessary)
  \return updated data vector (parameter v)
  */
  inline void Matrix2Quaternion(const Matrix& matrix, Vector& v)
  {
    matrix.calcQuaternion(v);
  }


  /*!\brief Rotation matrix -> roational quaternion (4 references to doubles)
  \param matrix 3x3 rotation transformation matrix
  \param q1 quaternion component 1 (vector component 1)
  \param q2 quaternion component 2 (vector component 2)
  \param q3 quaternion component 3 (vector component 3)
  \return updated double values
  */
  inline void Matrix2Quaternion(
    const Matrix& matrix
   ,double& q1
   ,double& q2
   ,double& q3
   ,double& q4)
  {
    matrix.calcQuaternion(q1,q2,q3,q4);
  }


  /*!\brief Rotation matrix -> axis/angle representation (four double values)
  \param matrix 3x3 rotation transformation matrix
  \param vec1 (axis vector component 1)
  \param vec2 (axis vector component 2)
  \param vec3 (axis vector component 2)
  \param angle (rotation angle)
  \return updated double values
  */
  inline void Matrix2AxisAngle(
    const Matrix& matrix
   ,double& vec1
   ,double& vec2
   ,double& vec3
   ,double& angle
  )
  {
    matrix.calcAxisAngle(vec1, vec2, vec3, angle);
  }


  /*!\brief Rotation matrix -> axis/angle representation (array of doubles/double)
  \param matrix 3x3 rotation transformation matrix
  \param axis array of double values holding rotation axis
  \param angle (rotation angle)
  \return updated axis/angle parameters
  */
  inline void Matrix2AxisAngle(
    const Matrix& matrix
   ,double axis[3]
   ,double& angle
  )
  {
    matrix.calcAxisAngle(axis[0], axis[1], axis[2], angle);
  }


  /*!\brief Rotation matrix -> axis/angle representation (vector/double)
  \param matrix 3x3 rotation transformation matrix
  \param axis vector holding rotation axis
  \param angle (rotation angle)
  \return updated axis/angle parameters
  */
  inline void Matrix2AxisAngle(
    const Matrix& matrix
   ,Vector& axis
   ,double& angle
  )
  {
    axis.Set(Vector(3));
    matrix.calcAxisAngle(axis[0], axis[1], axis[2], angle);
  }


  /*!\brief Rotation matrix -> axis/angle representation (four double values)
  \param vec1 (axis vector component 1)
  \param vec2 (axis vector component 2)
  \param vec3 (axis vector component 2)
  \param angle (rotation angle)
  \param matrix 3x3 rotation transformation matrix
  \return updated rotation matrix parameter
  */
  void AxisAngle2Matrix( 
    const double vec1 
   ,const double vec2
   ,const double vec3
   ,const double angle
   ,Matrix& matrix);

  /*!\brief Rotation matrix -> axis/angle representation (array of doubles
/double)
  \param axis rotation axis
  \param angle (rotation angle)
  \param matrix rotation transformation matrix
  \return updated rotation matrix parameter
  */
  inline void AxisAngle2Matrix( 
    const double axis[3]
   ,const double angle
   ,Matrix& matrix)
  {
    AxisAngle2Matrix(axis[0],axis[1],axis[2],angle, matrix);
  }


  /*!\brief Rotation matrix -> axis/angle representation (Vector/double)
  \param axis rotation axis
  \param angle (rotation angle)
  \param matrix rotation transformation matrix
  \return updated rotation matrix parameter
  */
  inline void AxisAngle2Matrix( 
    const Vector& axis
   ,const double angle
   ,Matrix& matrix)
  {
    AxisAngle2Matrix(axis(0), axis(1), axis(2), angle, matrix);
  }



  /** \brief LU-Decomposition.
      Given an n x n Matrix A, this routine replaces it by the LU decomposition
      of a rowwise permutation of itself. A is output, arranged as in
      equation (2.3.14) of Press et al. (1986); Indx is an ouput Vector which
      records the row permutation effected by partial pivoting. This routine is
      used in combination with LU_BackSub to solve linear equations or invert
      a Matrix. [Adapted from LUDCMP of Press et al. (1986).]
      \param A square matrix; replaced by LU decomposition of permutation of
               A on output
      \param Indx permutation index Vector */
  void LU_Decomp (Matrix & A, Vector & Indx, bool mod = false);

  /** \brief LU Backsubstitution.
      Solves the set of n linear equations Ax=b. Here A is input, not as the
      Matrix A but rather as its LU decomposition, determined by the function
      LU_Decomp. b is input as the right-hand side Vector b, and returns with
      the solution Vector x. A and Indx are not modified by this function and
      can be left in place for successive calls with different right-hand
      sides b. This routine takes into account the posssibility that B will
      begin with many zero elements, so it is efficient for use in Matrix
      inversions.
      \param A LU decomposition of permutation of A
      \param Indx permutation index vector
      \param b right-hand side vector b; replaced by solution x of Ax=b
               on output */
  void LU_BackSub(Matrix & A, Vector & Indx, Vector & b, bool mod = false);



  /** 
    \brief Cross product operator.
  
    If a cross product of the argument vector is conducted with another
    vector, the result ís also obtained if the returned matrix is
    multiplied with the right hand side vector.
  */
  Matrix CrossMatrix(const Vector& v);

  /** \brief Returns a unity matrix.
      The elements of the returned matrix are zero except
      the ones on the main diagonal which are 1.0.
      \param dim matrix dimension */
  Matrix Id(int dim);

  /** Returns a diagonal matrix initialized with the given vector.
      The elements of the matrix diagonal will be set with the content of the
      given vector. Therefore, the dimension of the resulting quadratic matrix
      is equal to the size of the argument vector.
      \param vec Vector that contains the new diagonal elements. */
  Matrix Diag(const Vector& vec);

  /** \brief Rotation matrix for rotation about axis 1 (clockwise vector 
  rotation or counterclockwise coordinate frame rotation).
      The resulting rotation matrix for a rotation about
      the axis 1 (usually x-axis in a righthanded cartesian
      coordinate system) depends on the given angle.
      \param angle rotation angle. */
  Matrix R_x(double angle);

  /** \brief transformation matrix for rotation about axis 1
   \param angle rotation angle of counterclockwise frame rotation
   \param m transformation matrix
   \return changed transformation matrix parameter*/
  void R_x(double, Matrix&);

  /** \brief Rotation matrix for rotation about axis 2 (clockwise vector 
  rotation or counterclockwise coordinate frame rotation).
      The resulting rotation matrix for a rotation about
      the axis 2 (usually y-axis in a righthanded cartesian
      coordinate system) depends on the given angle.
      \param angle rotation angle */
  Matrix R_y(double angle);

  /** \brief transformation matrix for rotation about axis 2
   \param angle rotation angle of counterclockwise frame rotation
   \param m transformation matrix
   \return changed transformation matrix parameter*/
  void R_y(double, Matrix&);

  /** \brief Rotation matrix for rotation about axis 3 (clockwise vector 
  rotation or counterclockwise coordinate frame rotation).
      The resulting rotation matrix for a rotation about
      the axis 3 (usually z-axis in a righthanded cartesian
      coordinate system) depends on the given angle.
      \param angle rotation angle */
  Matrix R_z(double angle);

  /** \brief transformation matrix for rotation about axis 3
    \param angle rotation angle of counterclockwise frame rotation
    \param m transformation matrix
    \return changed transformation matrix parameter*/
  void R_z(double angle, Matrix& m);

  /** Returns the transposition of the given matrix. */
  Matrix Transp(const Matrix & mat);

  /** Returns the inverse of the given matrix */
  Matrix Inv(const Matrix & mat);

  /** \brief Multiplication of scalar with matrix.
      \param value scalar
      \param mat matrix
      \returns result */
  Matrix operator*(double value, const Matrix & mat);

  /** \brief Multiplication of matrix with scalar.
      \param mat matrix
      \param value scalar
      \returns result */
  Matrix operator*(const Matrix & mat, double value);

  /** \brief Division of matrix by scalar.
      \param mat divident (matrix)
      \param value divisor (scalar)
      \returns result */
  Matrix operator/(const Matrix & mat, double value);

  /** \brief Matrix negation (all elements are multiplied with -1.0).
      \param mat matrix
      \returns result */
  Matrix operator-(const Matrix & mat);

  /** \brief Addition of two matrices.
      \param left left matrix
      \param right right matrix 
      \returns result */
  Matrix operator+(const Matrix & left, const Matrix & right);

  /** Returns the result of a matrix subtraction. */
  Matrix operator-(const Matrix & left, const Matrix & right);

  /** Returns the result of matrix multiplication operation. */
  Matrix operator*(const Matrix & left, const Matrix & right);

  void MultiplyMatrices(const Matrix& left, const Matrix& right, Matrix& final);

  /// matrix vector multiplication
  void MultiplyMatrixVector(const Matrix& A, const Vector& b, Vector& c);

  /** Returns the result of a matrix vector multplication operation. */
  Vector operator*(const Matrix & mat, const Vector & vec);

  /** Returns the result of a matrix vector multplication operation. */
  Vector operator*(const Vector & vec, const Matrix & mat);

  /** The returned matrix is the result of a dyadic product of two vectors. */
  Matrix Dyadic(const Vector & left, const Vector & right);

  /**
  \brief Sum of absolute element differences.

  Matrix one is compared with matrix two and the absoulte element differences
  are summed up and returned.

  \param A matrix one
  \param B matric two which should be compared to matrix one
  \return Sum of matrices difference values.
  */
  double Difference(const Matrix& A, const Matrix& B);

  /** Stream output operator. */
  std::ostream & operator<<(std::ostream & os, const Matrix & mat);

} // end of namespace MathLib

#endif // MATRIX_H

