/*!
\file MLMatrix.cpp
\version 1.7
\brief A class for mathematical matrices and functions for matrices (Matrix).
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 21.09.2006 (AO)

\b Purpose:
   Like the class vector the class declarated in this file (Matrix) is heavily
   used in the library SpaceLib to perform different transformations of vectors
   (distance and velocity). MathLib's member class Matrix provides some useful
   functions for this. Of course the basic Matrix operations are possible too.

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
#include <math.h>
#include <cassert>

//SpaceLib
#include <Global.h>
#include <MLVector.h>
#include <MLMatrix.h>

using namespace std;

namespace MathLib
{
//===== member functions ===================================================

  Matrix::Matrix(void)
  : n(0),
    m(0)
  {
    M = NULL;
  }

//--------------------------------------------------------------------------
  Matrix::Matrix(int dim_row, int dim_col)
  : n(dim_row),
    m(dim_col)
  {
    int row, col;

    M = new double*[dim_row];

    for (row = 0; row < n; ++row)
      M[row] = new double[dim_col];

    for (row = 0; row < n; ++row)
    {
      for (col = 0; col < m; ++col)
        M[row][col] = 0.0;
    }
  }

//--------------------------------------------------------------------------
  Matrix::Matrix(const Matrix & mat)
  {
    int row, col;

    n = mat.n;
    m = mat.m;

    M = new double *[n];

    for (row = 0; row < n; ++row)
      M[row] = new double[m];

    for (row = 0; row < n; ++row)
    {
      for (col = 0; col < m; ++col)
        M[row][col] = mat.M[row][col];
    }
  }

//--------------------------------------------------------------------------
  Matrix::Matrix(const double *array, int dim_row, int dim_col)
  : n(dim_row),
    m(dim_col)
  {
    int row, col;

    M = new double *[n];

    for (row = 0; row < n; ++row)
      M[row] = new double[m];

    for (row = 0; row < n; ++row)
    {
      for (col = 0; col < m; ++col)
        M[row][col] = array[row * dim_col + col];
    }
  }

//--------------------------------------------------------------------------
  Matrix::~Matrix()
  {
    for (int row = 0; row < n; ++row)
      delete[]M[row];
    delete[]M;
  }

//--------------------------------------------------------------------------
  Matrix & Matrix::operator=(const double value)
  {
    for (int row = 0; row < n; ++row)
      for (int col = 0; col < n; ++col)
        M[row][col] = value;

    return (*this);
  }

//--------------------------------------------------------------------------
  Matrix & Matrix::operator=(const Matrix & mat)
  {
    int row, col;

    if (this == &mat) // check identity
      return (*this);

    if (M == NULL) // allocate memory if necessary
    {
      n = mat.n;
      m = mat.m;
      M = new double*[n];

      for (row = 0; row < n; ++row)
        M[row] = new double[m];
    }

#ifdef development_version
    if (n != mat.n || m != mat.m)
      FatalError("Matrix::operator=() => incompatible sizes");
#endif

    for (row = 0; row < n; ++row)
      for (col = 0; col < m; ++col)
        M[row][col] = mat.M[row][col];

    return *this;
  }

//--------------------------------------------------------------------------
  Vector Matrix::Col(int col) const
  {
#ifdef development_version
    if ((col<0) || (col>=m))
      FatalError("Matrix::Col() => Index out of range");
#endif

    Vector result(n);
    for (int row = 0; row < n; ++row)
      result.v[row] = M[row][col];

    return result;
  }

//--------------------------------------------------------------------------
  Vector Matrix::Row(int row) const
  {
#ifdef development_version
    if ((row<0) || (row>=n))
      FatalError("Matrix::Row() => Index out of range");
#endif

    Vector result(m);
    for (int col = 0; col < m; ++col)
        result.v[col] = M[row][col];

    return result;
  }

//--------------------------------------------------------------------------
  Vector Matrix::Diag(void) const
  {
#ifdef development_version
    if (n != m)
      FatalError("Matrix::Diag() => Not a quadratic matrix");
#endif

    Vector result(n);
    for (int i = 0; i < n; ++i)
      result.v[i] = M[i][i];

    return result;
  }

//--------------------------------------------------------------------------
  double Matrix::Trace(void) const
  {
#ifdef development_version
    if (n != m)
      FatalError("Matrix::Trace() => Not a quadratic matrix");
#endif

    double result = 0.0;
    for (int i = 0; i < n; ++i)
      result += M[i][i];

    return result;

  }



  //--------------------------------------------------------------------------
  void Matrix::calcQuaternion(double& q1, double& q2, double& q3, double& q4)
    const
  {
    double val, t;
    int q( 1 );

    // begin with q1

    val = ( 1.0 + M[0][0] - M[1][1] - M[2][2] );

    // calculate a value for q2
    t = ( 1.0 - M[0][0] + M[1][1] - M[2][2] );

    // if q2 > q1 change index of element to calculate the other ones with
    if ( t > val )
    {
      val = t;
      q = 2; // now q2 is the primary element
    }

    // calculate a value for q3
    t = ( 1.0 - M[0][0] - M[1][1] + M[2][2] );

    // if q3 > old value, change index of element to calculate the other ones with
    if ( t > val )
    {
      val = t;
      q = 3; // now q3 is the primary element
    }

    // calculate a value for q4
    t = ( 1.0 + M[0][0] + M[1][1] + M[2][2] );

    // if q3 > old value, change index of element to calculate the other ones with
    if ( t > val )
    {
      val = t;
      q = 4; // now q4 is the primary element
    }

    // after estimating the biggest quaternion element, the remaining elements
    // can be computed
    switch ( q )
    {
      // compute q1 and from then the remaining
      case 1:
        val = 1.0 / ( 4.0 * ( q1 = 0.5 * sqrt( val ) ) );
        q2 = val * ( M[0][1] + M[1][0] );
        q3 = val * ( M[0][2] + M[2][0] );
        q4 = val * ( M[1][2] - M[2][1] );
      break;

        // the same][ if q2 should be used
      case 2:
        val = 1.0 / ( 4.0 * ( q2 = 0.5 * sqrt( val ) ) );
        q1 = val * ( M[0][1] + M[1][0] );
        q3 = val * ( M[1][2] + M[2][1] );
        q4 = val * ( M[2][0] - M[0][2] );
      break;

        // the same][ if q3 should be used
      case 3:
        val = 1.0 / ( 4.0 * ( q3 = 0.5 * sqrt( val ) ) );
        q1 = val * ( M[0][2] + M[2][0] );
        q2 = val * ( M[1][2] + M[2][1] );
        q4 = val * ( M[0][1] - M[1][0] );
      break;

        // the same][ if q4 should be used
      case 4:
        val = 1.0 / ( 4.0 * ( q4 = 0.5 * sqrt( val ) ) );
        q1 = val * ( M[1][2] - M[2][1] );
        q2 = val * ( M[2][0] - M[0][2] );
        q3 = val * ( M[0][1] - M[1][0] );
      break;

      default:
        FatalError( "void CAttitude::Matrix2Quaternion(void) => This part of the switch statement should have never been reached!" );
      break;
    }
  }


//--------------------------------------------------------------------------
  void Matrix::calcQuaternion(Vector& q) const
  {
    double q1,q2,q3,q4;
    calcQuaternion(q1,q2,q3,q4);
    q.Set(Vector(q1,q2,q3,q4));
  }


//--------------------------------------------------------------------------
  void Matrix::calcAxisAngle(double& v1, double& v2, double& v3, double& angle)
    const
  {
    if (m==3 && n==3)
    {
      const double cosp( 0.5 * ( Trace() - 1.0 ) );
      if ( abs( cosp ) <= 1.0 )
      {
        const double twosinp( 2.0 * sqrt( 1.0 - cosp * cosp ) );
        angle = acos( cosp );
        v1 = (M[1][2] - M[2][1] ) / twosinp;
        v2 = (M[2][0] - M[0][2] ) / twosinp;
        v3 = (M[0][1] - M[1][0] ) / twosinp;
      }
      else
      {
        v1 = v2 = v3 = angle = 0.0;
      }
    }
    else cout << "Matrix::calcAxisAngle(double& v1, double& v2, double& v3,\
 double& angle): rows or/and columns of matrix != 3!\n";
  }


//--------------------------------------------------------------------------
  void Matrix::calcAxisAngle( Vector& vec, double& angle) const
  {
    calcAxisAngle(vec(0),vec(1),vec(2),angle);
  }



//-----------------------------------------------------------------------------
  void Matrix::calcAxisAngle(double axis[3], double& angle) const 
  {
    calcAxisAngle(axis[0], axis[1], axis[3], angle);
  }



//--------------------------------------------------------------------------
  Matrix Matrix::slice(int first_row, int last_row, int first_col, int last_col)
  {
#ifdef development_version
    if (first_row < 0 || last_row < first_row || n - 1 < last_row ||
        first_col < 0 || last_col < first_col || m - 1 < last_col)
      FatalError("Matrix::slice() => Invalid argument(s)");
#endif

    Matrix result(last_row - first_row + 1, last_col - first_col + 1);
    for (int row = 0; row <= last_row - first_row; ++row)
      for (int col = 0; col <= last_col - first_col; ++col)
        result.M[row][col] = M[row + first_row][col + first_col];

    return result;
  }

//--------------------------------------------------------------------------
  void Matrix::SetCol(int col, const Vector & colvec)
  {
#ifdef development_version
    if (colvec.n != n)
      FatalError("Matrix::SetCol() => Incompatible sizes");
    if (col < 0 || m <= col)
      FatalError("Matrix::SetCol() => Column index out of range");
#endif

    for (int row = 0; row < n; ++row)
      M[row][col] = colvec.v[row];
  }

//--------------------------------------------------------------------------
  void Matrix::SetRow(int row, const Vector & rowvec)
  {
#ifdef development_version
    if (rowvec.n != m)
      FatalError("Matrix::Row() => Incompatible sizes");
    if (row < 0 || n <= row)
      FatalError("Matrix::Row() => Row index out of range");
#endif

    for (int col = 0; col < m; ++col)
      M[row][col] = rowvec.v[col];
  }

//--------------------------------------------------------------------------
  void Matrix::SetSubMatrix(int row, int col, const Matrix & mat)
  {
    if (n - row - mat.n >= 0 && m - col - mat.m >= 0)
    {
      for (int i = 0; i < mat.n; ++i)
        for (int j = 0; j < mat.m; ++j)
          M[row+i][col+j]=mat.M[i][j];
    }
  }

//===== matrix functions ===================================================

  void AxisAngle2Matrix( 
    const double vec1 
   ,const double vec2
   ,const double vec3
   ,const double angle
   ,Matrix&      matrix)
  {
    const double sinp( sin( angle ) );
    const double cosp( cos( angle ) );

    const double omcp(1.0 - cosp); // one minus cos(phi) 
    const double e1(vec1);
    const double e2(vec2); 
    const double e3(vec3);
    const double e1e2(e1*e2);
    const double e1e3(e1*e3); 
    const double e2e3(e2*e3);

    matrix(0,0) = cosp + e1*e1*omcp; 
    matrix(0,1) = e1e2*omcp+e3*sinp;
    matrix(0,2) = e1e3*omcp-e2*sinp;
    matrix(1,0) = matrix(0,1)-2.0*e3*sinp;
    matrix(1,1) = cosp + e2*e2*omcp;
    matrix(1,2) = e2e3*omcp+e1*sinp;
    matrix(2,0) = matrix(0,2)+2.0*e2*sinp;
    matrix(2,1) = matrix(1,2)-2.0*e1*sinp;
    matrix(2,2) = cosp + e3*e3*omcp;


    // much smaller in source code but probably slower
    /*!
    \todo check performance of both alternatives and drop slower one
    */
/*
    const Vector e( vec1,vec2,vec3 );
    matrix = cosp * Id( 3 ) + ( 1.0 - cosp ) * Dyadic( e, e ) 
      - sinp * CrossMatrix(2);
    */
  }

  //---------------------------------------------------------------------------
  void LU_Decomp(Matrix & A, Vector & Indx, bool mod)
  {
    // Constants

    const int n = A.getNumberOfRows();
    const double tiny = 1.0e-20;  // A small number

    // Variables

    int i, j, imax = 0, k;
    double aAmax, Sum, Dum;
    Vector V(n);

    // Loop over rows to get scaling information

    for (i = 0; i < n; ++i)
    {
      aAmax = 0.0;
      for (j = 0; j < n; ++j)
        if (fabs(A(i, j)) > aAmax)
          aAmax = fabs(A(i, j));
      if (aAmax == 0.0) {
        // No nonzero largest element
        cerr << "ERROR: Singular Matrix A in LU_Decomp";
        exit(1);
      };

      V(i) = 1.0 / aAmax;       // V stores the implicit scaling of each row
    };

    // Loop over columns of Crout's method

    for (j = 0; j < n; ++j) 
    {
      if (j > 0) {
        for (i = 0; i < j; ++i) { // This is equation 2.3.12 except for i=j
          Sum = A(i, j);
          if (i > 0) {
            for (k = 0; k < i; ++k)
              Sum -= A(i, k) * A(k, j);
            A(i, j) = Sum;
          };
        };
      };

      aAmax = 0.0;              // Initialize for the search of the largest
      // pivot element

      for (i = j; i < n; ++i) { // This is i=j of equation 2.3.12 and
        Sum = A(i, j);          // i=j+1..N of equation 2.3.13
        if (j > 0) {
          for (k = 0; k < j; ++k)
            Sum -= A(i, k) * A(k, j);
          A(i, j) = Sum;
        };
        Dum = V(i) * fabs(Sum); // Figure of merit for the pivot
        if (Dum >= aAmax) {     // Is it better than the best so far ?
          imax = i;
          aAmax = Dum;
        };
      };

      if (j != imax) {          // Do we need to interchange rows?
        for (k = 0; k < n; ++k) { // Yes, do so ...
          Dum = A(imax, k);
          A(imax, k) = A(j, k);
          A(j, k) = Dum;
        }
        V(imax) = V(j);         // Also interchange the scale factor
      };

      Indx(j) = imax;

      if (j != n - 1) {         // Now finally devide by the pivot element
        if (!mod && (A(j, j) == 0.0)) {   // If the pivot element is zero the Matrix
          A(j, j) = tiny;       // is singular (at least to the precision of
        };                      // the algorithm). For some applications on
        Dum = 1.0 / ((mod&&A(j,j)==0.0)?1.0:A(j, j));    // singular matrices, it is desirable to
        for (i = j + 1; i < n; ++i) 
        { // substitude tiny for zero.
          A(i, j) = A(i, j) * Dum;
        };
      };

    };                          // Go back for the next column in the reduction

    if (A(n - 1, n - 1) == 0.0)
      A(n - 1, n - 1) = tiny;
  }

//--------------------------------------------------------------------------
  void LU_BackSub(Matrix & A, Vector & Indx, Vector & b, bool mod)
  {
    // Constants

    const int n = A.getNumberOfRows();

    // Local variables

    int ii, i, ll, j;
    double Sum;

    //
    // Start
    //

    ii = -1;                    // When ii is set to a nonegative value, it will
    // become the first nonvanishing element of B.
    for (i = 0; i < n; ++i) 
    {   // We now do the forward substitution.
      ll = (int)Indx(i);        // The only wrinkle is to unscramble the
      Sum = b(ll);              // permutation as we go.
      b(ll) = b(i);
      if (ii != -1) {
        for (j = ii; j < i; ++j)
          Sum -= A(i, j) * b(j);
      }
      else {
        if (Sum != 0.0)
          ii = i;               // A nonzero element was encountered, so from
      };                        // now on we will have to do the sums in the
      b(i) = Sum;               // loop above.
    };

    for (i = n - 1; i >= 0; --i) 
    {  // Now we do the backsubstitution, eqn 2.3.7.
      Sum = b(i);
      if (i < n - 1) {
        for (j = i + 1; j < n; ++j) 
        {
          Sum = Sum - A(i, j) * b(j);
        };
      };
      b(i) = (mod&&(A(i,i)==0.0)?0.0:(Sum / A(i, i)));     // Store a component of the solution Vector X.
    };
  }





//===== friend functions ===================================================

//--------------------------------------------------------------------------
  Matrix MathLib::CrossMatrix(const Vector& v)
  {
    if (v.getDimension() == 3)
    {
      Matrix C(3,3);
      C.M[0][1] = -(C.M[1][0] = v(2));
      C.M[2][0] = -(C.M[0][2] = v(1));
      C.M[1][2] = -(C.M[2][1] = v(0));
      return C;
    }
    else return Id(3);
  }


//--------------------------------------------------------------------------
  Matrix MathLib::Id(int dim)
  {
    Matrix result(dim, dim);
    for (int i = 0; i < dim; ++i)
      result.M[i][i] = 1.0;
    
    return result;
  }

//--------------------------------------------------------------------------
  Matrix MathLib::Diag(const Vector & vec)
  {
    Matrix result(vec.n, vec.n);
    for (int i = 0; i < vec.n; ++i)
      result.M[i][i] = vec.v[i];

    return result;
  }

//--------------------------------------------------------------------------
  Matrix MathLib::R_x(double angle)
  {
    const double C = cos(angle);
    const double S = sin(angle);

    Matrix U(3, 3);
    U.M[0][0] = 1.0;
    U.M[0][1] = 0.0;
    U.M[0][2] = 0.0;
    U.M[1][0] = 0.0;
    U.M[1][1] = +C;
    U.M[1][2] = +S;
    U.M[2][0] = 0.0;
    U.M[2][1] = -S;
    U.M[2][2] = +C;

    return U;
  } // Matrix MathLib::R_x(double angle)

  //----------------------------------------------------------------------------
  void R_x(double angle, Matrix& M)
  {
#ifdef development_version
    assert( (M.getNumberOfRows()==3) && (M.isSquare() == true) );
#endif // development_version
    const double C = cos(angle);
    const double S = sin(angle);

    M.M[0][0] = 1.0;
    M.M[0][1] = 0.0;
    M.M[0][2] = 0.0;
    M.M[1][0] = 0.0;
    M.M[1][1] = +C;
    M.M[1][2] = +S;
    M.M[2][0] = 0.0;
    M.M[2][1] = -S;
    M.M[2][2] = +C;
    
  } // void R_x(double, Matrix&)

//--------------------------------------------------------------------------
  Matrix MathLib::R_y(double angle)
  {
    const double C = cos(angle);
    const double S = sin(angle);

    Matrix U(3, 3);
    U.M[0][0] = +C;
    U.M[0][1] = 0.0;
    U.M[0][2] = -S;
    U.M[1][0] = 0.0;
    U.M[1][1] = 1.0;
    U.M[1][2] = 0.0;
    U.M[2][0] = +S;
    U.M[2][1] = 0.0;
    U.M[2][2] = +C;

    return U;
  } // Matrix MathLib::R_y(double angle)

  //----------------------------------------------------------------------------
  void R_y(double angle, Matrix& M)
  {
#ifdef development_version
    assert( (M.getNumberOfRows()==3) && (M.isSquare() == true) );
#endif // development_version
    const double C = cos(angle);
    const double S = sin(angle);

    M.M[0][0] = +C;
    M.M[0][1] = 0.0;
    M.M[0][2] = -S;
    M.M[1][0] = 0.0;
    M.M[1][1] = 1.0;
    M.M[1][2] = 0.0;
    M.M[2][0] = +S;
    M.M[2][1] = 0.0;
    M.M[2][2] = +C;
  } // void R_y(double angle, Matrix& M)

//--------------------------------------------------------------------------
  Matrix MathLib::R_z(double angle)
  {
    const double C = cos(angle);
    const double S = sin(angle);

    Matrix U(3, 3);
    U.M[0][0] = +C;
    U.M[0][1] = +S;
    U.M[0][2] = 0.0;
    U.M[1][0] = -S;
    U.M[1][1] = +C;
    U.M[1][2] = 0.0;
    U.M[2][0] = 0.0;
    U.M[2][1] = 0.0;
    U.M[2][2] = 1.0;

    return U;
  }  // Matrix MathLib::R_z(double angle)

  //----------------------------------------------------------------------------
  void R_z(double angle, Matrix& M)
  {
#ifdef development_version
    assert( (M.getNumberOfRows()==3) && (M.isSquare() == true) );
#endif // development_version
    const double C = cos(angle);
    const double S = sin(angle);

    M.M[0][0] = +C;
    M.M[0][1] = +S;
    M.M[0][2] = 0.0;
    M.M[1][0] = -S;
    M.M[1][1] = +C;
    M.M[1][2] = 0.0;
    M.M[2][0] = 0.0;
    M.M[2][1] = 0.0;
    M.M[2][2] = 1.0;
  } // void R_z(double angle, Matrix& M)

//--------------------------------------------------------------------------
  Matrix MathLib::Transp(const Matrix & mat)
  {
    Matrix T(mat.m, mat.n);
    for (int i = 0; i < T.n; ++i)
      for (int j = 0; j < T.m; ++j)
        T.M[i][j] = mat.M[j][i];

    return T;
  }

//--------------------------------------------------------------------------
  Matrix MathLib::Inv(const Matrix & mat) {

    const int n = mat.n;        // number of rows

    Matrix LU(n, n), Inverse(n, n);
    Vector b(n), Indx(n);

    if (mat.m != mat.n) {
      cerr << "ERROR: Invalid shape in Inv(Matrix)" << endl;
      exit(1);
    };

    // LU decomposition
    LU = mat;
    LU_Decomp(LU, Indx);

    // Solve Ax=b for  unit Vectors b_1..b_n
    for (int j = 0; j < n; ++j) 
    {
      b = 0.0;
      b(j) = 1.0;               // Set b to j-th unit Vector
      LU_BackSub(LU, Indx, b);  // Solve Ax=b
      Inverse.SetCol(j, b);     // Copy result
    };

    return Inverse;
  }

//--------------------------------------------------------------------------
  Matrix MathLib::operator*(double value, const Matrix & mat)
  {
    Matrix result(mat.n, mat.m);
    for (int row = 0; row < mat.n; ++row)
      for (int col = 0; col < mat.m; ++col)
        result.M[row][col] = value * mat.M[row][col];

    return result;
  }

//--------------------------------------------------------------------------
  Matrix MathLib::operator*(const Matrix & mat, double value)
  {
    return value * mat;
  }

//--------------------------------------------------------------------------
  Matrix MathLib::operator/(const Matrix & mat, double value)
  {
    Matrix result(mat.n, mat.m);
    for (int row = 0; row < mat.n; ++row)
      for (int col = 0; col < mat.m; ++col)
        result.M[row][col] = mat.M[row][col] / value;

    return result;
  }


//--------------------------------------------------------------------------
  Matrix MathLib::operator-(const Matrix & mat)
  {
    Matrix result(mat.n, mat.m);
    for (int row = 0; row < mat.n; ++row)
      for (int col = 0; col < mat.m; ++col)
        result.M[row][col] = -mat.M[row][col];

    return result;
  }

//--------------------------------------------------------------------------
  Matrix MathLib::operator+(const Matrix & left, const Matrix & right)
  {
#ifdef development_version
    if ((left.n != right.n) || (left.m != right.m))
      FatalError("Matrix::operator+() => Incompatible sizes");
#endif

    Matrix result(left.n, left.m);
    for (int row = 0; row < left.n; ++row)
      for (int col = 0; col < left.m; ++col)
        result.M[row][col] = left.M[row][col] + right.M[row][col];

    return result;
  }

//--------------------------------------------------------------------------
  Matrix MathLib::operator-(const Matrix & left, const Matrix & right)
  {
#ifdef development_version
    if ((left.n != right.n) || (left.m != right.m))
      FatalError("Matrix::operator-() => Incompatible sizes");
#endif

    Matrix result(left.n, left.m);
    for (int row = 0; row < left.n; ++row)
      for (int col = 0; col < left.m; ++col)
        result.M[row][col] = left.M[row][col] - right.M[row][col];

    return result;
  }

//--------------------------------------------------------------------------
  Matrix MathLib::operator*(const Matrix & left, const Matrix & right)
  {
#ifdef development_version
    if ((left.m != right.n))
      FatalError("Matrix::operator*(Matrix,Matrix) => Incompatible sizes");
#endif

    Matrix result(left.n, right.m);
    double sum;
    for (int row = 0; row < left.n; ++row)
      for (int col = 0; col < right.m; ++col)
      {
        sum = 0.0;
        for (int k = 0; k < left.m; ++k)
          sum += left.M[row][k] * right.M[k][col];
        result.M[row][col] = sum;
      }
    return result;
  }


  //----------------------------------------------------------------------------
  void MultiplyMatrices(const Matrix& l, const Matrix& r, Matrix& f)
  {
#ifdef development_version
    assert(l.m==r.n);
    assert(l.getNumberOfRows()==f.getNumberOfRows());
    assert(r.getNumberOfColumns()==f.getNumberOfColumns());
#endif
    double sum;
    for (int row = 0; row < l.n; ++row)
      for (int col = 0; col < r.m; ++col)
      {
        sum = 0.0;
        for (int k = 0; k < l.m; ++k)
          sum += l.M[row][k] * r.M[k][col];
        f.M[row][col] = sum;
      }
  } // void MultiplyMatrices(const Matrix&, const Matrix&, Matrix&)

//--------------------------------------------------------------------------
  Vector MathLib::operator*(const Matrix & mat, const Vector & vec)
  {
#ifdef development_version
    assert(mat.m==vec.n);
    //if ((mat.m != vec.n))
    //  FatalError("Matrix::operator*(Matrix,Vector) => Incompatible sizes");
#endif

    Vector result(mat.n);
    double sum;
    for (int row = 0; row < mat.n; ++row)
    {
      sum = 0.0;
      for (int col = 0; col < mat.m; ++col)
        sum += mat.M[row][col] * vec.v[col];
      result.v[row] = sum;
    }
    return result;
  }


  //----------------------------------------------------------------------------
  void MultiplyMatrixVector(const Matrix& A, const Vector& b, Vector& c)
  {
#ifdef development_version
    assert(A.getNumberOfColumns() == b.size());
    assert(A.getNumberOfRows() == c.size());
#endif
    double sum;
    for (int row = 0; row < A.n; ++row)
    {
      sum = 0.0;
      for (int col = 0; col < A.m; ++col)
        sum += A.M[row][col] * b.v[col];
      c.v[row] = sum;
    }
  } // void MultiplyVector(const Matrix& A, const Vector& b, Vector& c)


//--------------------------------------------------------------------------
  Vector MathLib::operator*(const Vector & vec, const Matrix & mat)
  {
#ifdef development_version
    if ((mat.n != vec.n))
      FatalError("Matrix::operator*(Vector,Matrix) => Incompatible sizes");
#endif

    Vector result(mat.m);
    double sum;
    for (int col = 0; col < mat.m; ++col)
    {
      sum = 0.0;
      for (int row = 0; row < mat.n; ++row)
        sum += vec.v[row] * mat.M[row][col];
      result.v[col] = sum;
    }
    return result;
  }

//--------------------------------------------------------------------------
 Matrix MathLib::Dyadic(const Vector & left, const Vector & right)
  {

    Matrix result(left.n, right.n);
    for (int row = 0; row < left.n; ++row)
      for (int col = 0; col < right.n; ++col)
        result.M[row][col] = left.v[row] * right.v[col];

    return result;
  }

//--------------------------------------------------------------------------
  ostream & MathLib::operator<<(ostream & os, const Matrix & mat) {

    int w = os.width();
    int p = os.precision();
    for (int row = 0; row < mat.n; ++row)
    {
      for (int col = 0; col < mat.m; ++col)
        os << setw(p+7) << setprecision(p) << right << mat.M[row][col];// << setw(1) << " ";
      os << endl;
    }
    return os;
  }

//--------------------------------------------------------------------------
  // helper function
  double Difference(const Matrix& A, const Matrix& B)
  {
    double diff(0.0);
    if ( (A.getNumberOfRows()==B.getNumberOfRows()) && 
         (A.getNumberOfColumns()==B.getNumberOfColumns()) ) 
    {

      const int rows(A.getNumberOfRows());
      const int columns(A.getNumberOfColumns());

      for (int i = 0; i<rows; ++i)
        for (int j = 0; j<columns; ++j) diff+= fabs(A(i,j)-B(i,j));

    } else FatalError("MatrixDifference(const Matrix&, const Matrix&) => Dimensions not equal!");
    return diff;
  }




} // end of namespace MathLib
