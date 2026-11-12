/*!
\file MLPolynomial.h
\version 1.7
\brief A class for to handle real polynomials.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 28.08.2005 (AO)

\b Purpose:


\b Note:
This software is protected by national and international copyright.
Any unauthorized use, reproduction or modificaton is unlawful and
will be prosecuted. Commercial and non-private application of the
software in any form is strictly prohibited unless otherwise granted
by the authors.
The code is provided without any warranty; without even the implied
warranty of merchantibility or fitness for a particular purpose.\n
\n
References:\n
(1) Sedgewick Robert, "Algorithmen in C++", Addison-Wesley (Deutschland),
1995\n
(2) Press William H., "Numerical Recipes in C++, Second Edition",Cambridge 
University Press, 2003\n

\b Copyright:
(c) 2005 Andreas Ohndorf, Bernd Dachwald
*/
#ifndef MLPOLYNOMIAL_H
#define MLPOLYNOMIAL_H

#pragma warning (disable:4290) // ignored exception declaration

// STL
#include <stdexcept>
#include <vector>
#include <valarray>
#include <algorithm>


// SpaceLib
#include <MLVector.h>
#include <MLMatrix.h>

namespace MathLib {

  typedef enum {
    NON_ZERO_START = 0
    ,ZERO_START     = 1
  }  PolyApproxType;

  /*! \brief Maximum degree of polynomial usable for data approximation.
  \sa CalculateApproximationPolynomials
  */
  const unsigned int MAX_APPROX_POLYNOMIAL_DEGREE(50);

  /*!
  \brief Real value, double precision polynomials with maximum degree of 10.
  */
  class Polynomial
  {

    MathLib::Vector m_coeff;


  public:

    //---------------------------------------------------------------------------
    /*!
    \brief Standard constructor.
    */
    Polynomial(void):m_coeff() {};
    //---------------------------------------------------------------------------
    /*!
    \brief Copy constructor.

    \param to_copy The polynomial to copy.
    */
    Polynomial(const Polynomial& to_copy):m_coeff(to_copy.m_coeff) {};
    //---------------------------------------------------------------------------
    /*!
    \brief Initialize polynomial with given degree and coefficents.

    \param degree Polynomial should have this degree.
    \param coeff A Pointer to a double value array representing the polynomials
    coefficents.
    */
    Polynomial(int degree, double* coeff);
    //---------------------------------------------------------------------------
    /*!
    \brief Initialize polynomial with given coefficents (valarray).

    \param data data structure holding polynomial coefficents
    \param begin index of first parameter in data structure
    \param end index of last polynomial coefficent
    */
    template<class T>
      Polynomial(const T& data, const uint begin, const uint end)
      : m_coeff(begin<end?end-begin+1:0)
    {
      if (int(begin)<int(end))
      {
        for (uint i=0;i<(end-begin+1);++i) m_coeff[i] = data[begin+i];
      } else m_coeff.Set(Vector());
    }
    //---------------------------------------------------------------------------
    /*!
    \brief Assign operator.
    */
    Polynomial& operator=(const Polynomial& right);
    //---------------------------------------------------------------------------
    /*!
    \brief Destructor.
    */
    ~Polynomial(void) {};
    //---------------------------------------------------------------------------
    /*!
    \brief Returns polynomial value for given x.

    \param x Abscissa data we want the polynomials evaluated at.
    \return Value of polynomial at given input x.
    */
    double getValue(const double x) const;

    //---------------------------------------------------------------------------
    /*!
    \brief Returns value of polynomials first derivative for given x.

    \param x Abscissa data we want the polynomials first derivative evaluated
    at.
    \return Value of polynomials first derivative at given input x.
    */
    double firstDerivative(const double x) const;

    //---------------------------------------------------------------------------
    /*!
    \brief Returns value of polynomials second derivative for given x.

    \param x Abscissa data we want the polynomials second derivative evaluated
    at.
    \return Value of polynomials second derivative at given input x.
    */
    double secondDerivative(const double x) const;


  }; // class Polynomial


  /*!
  * \brief Initialization of approximation input data structures from given
  *segment size and polynomial degress for all segments.
  *
  * \param data original data values
  * \param segm output data structure of segmentation
  * \param degrees output data structure of polynomial degree
  * \param segment_size input segment size for all segments
  * \param poly_degree polynomial degree of all polynomials (all segments)
  * \sa SetUpApproximationArrays
  */
  void SetUpApproximationArrays(
    const std::valarray<double>& data
    ,      std::valarray<unsigned int>& segm
    ,      std::valarray<unsigned int>& degrees
    , const unsigned int segment_size
    , const unsigned int poly_degree);


  /*! \brief Initializes approximation segmentation data arrays from
  * pointers to segmentation information.
  *
  * The approximation functions expect their input i a certain type of
  * data structures (ie. std::valarray references). This helper function
  * creates these data structures and takes "casual" pointers to data
  * scopes as input.
  *
  *\param segs segment count (the number of segments).
  *\param seg_end_indices pointer to unsigned int array of indices
  *\param seg_poly_deg pointer to unsigned int array of segment polynomial degrees
  *\param segm new segmentation data structure
  *\param degrees new polynomial degree data structure
  */
  void SetUpApproximationArrays(
    unsigned int segs
    , unsigned int* seg_end_indices
    , unsigned int* seg_poly_deg
    ,      std::valarray<unsigned int>& segm
    ,      std::valarray<unsigned int>& degrees 
    );


  //---------------------------------------------------------------------------
  /*!
  * \brief Conversion of segmentation from "without start index" to 
  * "with start index 0"
  * \param i input data structure
  * \param o output data structure
  * \return modified output data structure
  */
  template<class T, class U>
    inline void CompModSegVec( const T& i, U& o)
  {
    if (i.size()>0)
    {
      o.resize(i.size()+1); // set size of new segmentation vector
      o[0]=0;               // start first segment at index 0
      for (uint idx=0;idx<i.size();++idx) o[idx+1]=i[idx]; // copy indices
    }
  };


  //-----------------------------------------------------------------------------
  /*!
  * \brief Copies desired data scope into new data struc
  * \param i input data structure
  * \param b first index of desired data scope
  * \param e last index of desired data scope
  * \param o output data structure
  * \return modified output data structure
  */
  template<class T, class U>
    inline void CalcDataVecs(const T& i, const uint b, const uint e, U& o)
  {
    if (i.size()>0 && int(e)<i.size() && int(b)<int(e))
    {
      o.resize(e-b+1);
      for (uint idx=0;idx<o.size();++idx) o[idx]=i[b+idx];
    }
  }


  //-----------------------------------------------------------------------------
  /*! \brief Sets components of a vector according to given input data structure.
  *
  * \param i input data structure with data to set vector with
  * \param y vector to modify
  * \param b start index in vector to modifiy components
  */
  template<class T>
    inline void SetVector(const T& i, Vector& y, const uint b)
  {
    if (i.size()>0 && y.getDimension()>0 && b+i.size()<=uint(y.getDimension()))
      for (uint idx=0;idx<i.size();++idx) y[b+idx]=i[idx];
  }


  //---------------------------------------------------------------------------
  /*!
  * \brief Computes system matrix for polynomial approximation
  *
  * \param input input data (independend variable)
  * \param val   data to approximate
  * \param b start index in input and value data structures
  * \param e end index in input and value data structures
  * \param degree degree of desired approximation polynomial
  * \param value_vec right hand side vector of system equation
  * \return System matrix and modified value data structure to solve for
  *polynomial parameters
  */
  template<class T>
    Matrix ComputeSystemMatrix(
    const T& input
    ,const T& val
    ,const uint b
    ,const uint e
    ,const uint degree
    ,T& value_vec);


  //---------------------------------------------------------------------------
  /*!
  * \brief Constructs functional vector for system matrix creation from input
  * data structure (possible sub scope)
  * 
  * The functional vector has the same dimension as the input data vector and
  * contains functional values of the desired base function which is in this
  * case a polynomial of given degree. For example, if parameter degree has a
  * value of 3, all members of given input vector are raised to the power of
  * 3 and stored in a data structure (of input data structure type) that is the
  * final return value.
  * 
  * \param input Input data structure (must provide index - operator and size() - method
  * \param begin first index of data scope to use
  * \param end end index of data scope to use
  * \param degree desired degree of functional polynomial
  * \return Data structure filled with base functionals evaluated with input data.
  */
  template<class T>
    T CalcFuncVec(
    const T& input
    ,const uint begin
    ,const uint end
    ,const int degree)
  {
    if ((int(end-begin)>0) && (input.size()>0))
    {
      const uint sz(std::min(end-begin+1, input.size()));
      T ret_val(sz);
      for (uint i=0;i<sz;++i) ret_val[i] = pow(input[begin+i], degree);
      return ret_val;
    }
    return T();
  }

  //---------------------------------------------------------------------------
  /*!
  * \brief Constructs functional vector for system matrix creation from input
  * data structure (all data)
  *
  * \param input Input data structure (must provide index - 
  *operator and size() - method)
  * \return Data structure filled with base functionals evaluated with input 
  *data.
  * \param degree desired degree of functional polynomial
  */
  template<class T>
    inline T CalcFuncVec(const T& input, const int degree)
  {
    return CalcFuncVec(input, 0, input.size()-1, degree);
  }


  //---------------------------------------------------------------------------
  /*! \brief Wrapper function for ComputeSystemMatrix for complete data structs
  *
  * \param input input data (independend variable)
  * \param val   data to approximate
  * \param degree degree of desired approximation polynomial
  * \param value_vec right hand side vector of system equation
  * \return System matrix and modified value data structure to solve for
  *polynomial parameters
  */
  template<class T>
    inline Matrix ComputeSystemMatrix(
    const T& input
    ,const T& val
    ,const uint degree
    ,T& value_vec)
  {
    return ComputeSystemMatrix(input, val, 0, input.size()-1, degree, value_vec);
  }

  //-----------------------------------------------------------------------------
  /*!
  * \brief Maps input data structure onto interval (-1,1).
  * \param i input data structure
  * \param b first index of desired data scope
  * \param e last index of desired data scope
  * \param o output data structure
  * \return modified output data structure
  */
  template<class T, class U>
    void CalcModInput(const T& i, const uint b, const uint e, U& o)
  {
    if (i.size()>0 && int(e)<i.size() && int(b)<int(e))
    {
      const uint sz(e-b+1);
      o.resize(sz);
      const T::value_type divisor(i[e]-i[b]);
      for (uint idx=0;idx<sz;++idx) o[idx]=(2.0*i[b+idx]-i[b]-i[e])/divisor;
    }
  }


  //-----------------------------------------------------------------------------
  /*! \brief Calculates a polynomials approximation of given data pair set.
  *
  */
  void CalcPolyApprox(
    const std::valarray<double>&   iv    ///< independend variable
    ,const std::valarray<double>&  y    ///< data
    ,const std::valarray<uint>&    s    ///< segmentation
    ,const std::valarray<uint>&    d    ///< degree of segment's polynomials
    ,std::vector<Polynomial>&      p     ///< polynomials
    ,std::valarray<double>&        ny    ///< new approximation of data 
    );


  //---------------------------------------------------------------------------
  /*! \brief Calculates polynomial approximation data at given independent
  *variable knots (here, time).
  */
  void CalcPolyApproxData(
     const std::vector<Polynomial>& polys  ///< vector of approximation polynomials
    ,const std::valarray<double>&   t      ///< original time resolution
    ,const std::valarray<uint>&     s      ///< segmentation array
    ,const double                   rf     ///< resolution factor
    ,      std::valarray<double>&   nt     ///< new time resolution
    ,      std::valarray<double>&   nd     ///< new data array
    );


  //-----------------------------------------------------------------------------
  /*! \brief Calculates first derivative of polynomial approximation data 
  *at given independent variable knots (here, time).
  */
  void CalcPolyApproxFirstDerivative(
    const std::vector<Polynomial>& polys  ///< vector of approximation polynomials
    ,const std::valarray<double>&  t      ///< original time resolution
    ,const std::valarray<uint>&    s      ///< segmentation array
    ,const double                  rf     ///< resolution factor
    ,      std::valarray<double>&  nt     ///< new time resolution
    ,      std::valarray<double>&  dd     ///< new derivatives array
    );


  //-----------------------------------------------------------------------------
  /*! \brief Calculates second derivative of polynomial approximation data 
  *at given independent variable knots (here, time).
  */
  void CalcPolyApproxSecondDerivative(
    const std::vector<Polynomial>& polys  ///< vector of approximation polynomials
    ,const std::valarray<double>&  t      ///< original time resolution
    ,const std::valarray<uint>&    s      ///< segmentation array
    ,const double                  rf         ///< resolution factor
    ,      std::valarray<double>&  nt     ///< new time resolution
    ,      std::valarray<double>&  dd     ///< new derivatives array
    );


  //-----------------------------------------------------------------------------
  /*!
  * \brief Calculate set of second derivatives at to satisfy spline interpolation
  *criteria for given set of data pairs (x,y)
  * 
  * \param x independent variable data structure
  * \param y ordinate value data structure of given data pairs (x,y)
  * \param p reference to data structure for calculated second derivatives
  */
  template<class T>
  void MakeSpline(const T& x, const T& y, T& p){};


  //-----------------------------------------------------------------------------
  // helper function for EvaluateSplinePoint
  template<class T>
    inline T f(T x)
  {
    return x*x*x-x;
  }

  //-----------------------------------------------------------------------------
  /*! \brief evaluates spline at given abscissa value
  *
  * \param x array of knot abscissa values
  * \param y array of knot ordinate values
  * \param p second derivatives at knots
  * \param v input abscissa value for which the ordinate value is wanted
  * \return ordinate value of spline at given abscissa
  * \sa SplineApproximation
  */
  template<class T>
    inline double EvaluateSplinePoint(const T& x, const T& y, const T& p, double v)
  {
    if (v>x[x.size()-1]) v = x[x.size()-1];
    if (v<x[0]) v = x[0];
    uint i(0);
    while (v>x[i+1] && i<x.size()-1) ++i;
    if (v>x[i] && i==(x.size()-1) )
    {
      v = x[i];
      --i;
    }
    T::value_type u(x[i+1]-x[i]);
    T::value_type t((v-x[i])/u);
    return t*y[i+1]+(1.0-t)*y[i]+u*u*(f(t)*p[i+1]+f(1-t)*p[i])/6.0;
  }


  //-----------------------------------------------------------------------------
  // helper function for EvaluateSplinePointDerivative
  template<class T>
    T f_d(T x)
  {
    return 3.0*x*x-1;
  }

  //-----------------------------------------------------------------------------
  /*! \brief evaluates spline at given abscissa value
  *
  * \param x array of knot abscissa values
  * \param y array of knot ordinate values
  * \param p second derivatives at knots
  * \param v input abscissa value for which the ordinate value is wanted
  * \return ordinate value of spline at given abscissa
  * \sa SplineApproximation
  */
  template<class T>
    double EvaluateSplinePointDerivative(const T& x    , const T& y, const T& p
    ,double v)
  {
    if (v>x[x.size()-1]) v = x[x.size()-1];
    if (v<x[0]) v = x[0];
    uint i(0);
    while (v>x[i+1] && i<x.size()-1) ++i;
    if (v>x[i] && i==(x.size()-1)) v = x[i];
    T::value_type u(x[i+1]-x[i]);
    T::value_type z( (y[i+1]-y[i]) /u);
    T::value_type t((v-x[i])/u);
    return (y[i+1]-y[i])/u+u*(f_d(t)*p[i+1]-f_d(1.0-t)*p[i])/6.0;
  }


  //-----------------------------------------------------------------------------
  /*! \brief Resamples given set of independent variable input data values.
  * \param x data structure containing independent variable values
  * \param new_x reference to data structure, where new data should be stored
  * \param res resolution factor (ratio of new data value count to old one)
  * \return modified data structure new_x
  */
  template<class T, class U>
  void ResampleArray(const T& x, U& new_x, const double res)
  {
    if (x.size()>0)
    {
      const double span(x[x.size()-1]-x[0]);
      new_x.resize(int(ceil(res*(x.size()-1)+1)));
      const double inc(span/(new_x.size()-1));
      uint idx(0);
      do
      {
        new_x[idx]=x[0]+idx*inc;
        ++idx;
      } while (idx<new_x.size());
    }
  }


  //-----------------------------------------------------------------------------
  /*! \brief Evaluates spline interpolation polynomials at given values of
  *independent variable (and first derivative).
  *
  * \param x old independent variable data
  * \param y old ordinate data values
  * \param p second derivatives from MakeSpline
  * \param new_x new independend variable input (a must give, ie. from ResampleArray)
  * \param new_y new (through interpolation) estimated values
  * \param new_y_d new first derivatives at given input knot data
  * \return modified data structures new_y and new_y_d
  * \sa MakeSpline, ResampleArray
  */
  template<class T, class U>
    inline void EvaluateSpline(const T& x      , const T& y     , const T& p
                              ,const T& new_x  ,       U& new_y ,       U& new_y_d)
  {
    const uint sz(new_x.size());
    new_y.resize(sz);new_y_d.resize(sz);
    double t(0.0),s(0.0);
    for (uint i=0;i<sz;++i)
    {
      new_y[i]  =s=EvaluateSplinePoint(x,y,p,new_x[i]);      
      new_y_d[i]=t=EvaluateSplinePointDerivative(x,y,p,new_x[i]);
    }
  }


  //-----------------------------------------------------------------------------
  /*! \brief Evaluates spline interpolation polynomials at given values of
  *independent variable.
  *
  * \param x old independent variable data
  * \param y old ordinate data values
  * \param p second derivatives from MakeSpline
  * \param new_x new independend variable input (a must give, ie. from ResampleArray)
  * \param new_y new (through interpolation) estimated values
  * \return modified data structure new_y
  * \sa MakeSpline, ResampleArray
  */
  template<class T, class U>
    inline void EvaluateSpline(const T& x      , const T& y     , const T& p
                              ,const T& new_x  ,       U& new_y )
  {
    const uint sz(new_x.size());
    new_y.resize(sz);new_y_d.resize(sz);
    for (uint i=0;i<sz;++i) new_y[i]=EvaluateSplinePoint(x,y,p,new_x[i]);
  }


  //-----------------------------------------------------------------------------
  /*! \brief Evaluates spline interpolation polynomials at values of independent
  * variable resulting from resampling original input data with given resolution
  * factor.
  *
  * \param x old independent variable (IV) data
  * \param y old ordinate data values
  * \param p second derivatives from MakeSpline
  * \param new_x new independend variable input (will be calculated)
  * \param new_y new (through interpolation) estimated values
  * \param res resolution factor (=number of new IV values/number of old IV values)
  * \return modified data structures new_x, new_y
  * \sa MakeSpline, ResampleArray
  */
  template<class T, class U>
    inline void EvaluateSpline(const T& x    , const T& y    , const T& p
                              ,      T& new_x,       U& new_y, const double res)
  {
    ResampleArray(x, new_x, res);
    EvaluateSpline(x, y, p, new_x, new_y);
  }


  //-----------------------------------------------------------------------------
  /*! \brief Evaluates spline interpolation polynomials (and first derivatives)
  * at values of independent variable resulting from resampling original input
  * data with given resolution factor.
  *
  * \param x old independent variable (IV) data
  * \param y old ordinate data values
  * \param p second derivatives from MakeSpline
  * \param new_x new independend variable input (will be calculated)
  * \param new_y new (through interpolation) estimated values
  * \param new_y_d new (through interpolation) estimated values if first derivatives
  * \param res resolution factor (=number of new IV values/number of old IV values)
  * \return modified data structures new_x, new_y, new_y_d
  * \sa MakeSpline, ResampleArray
  */
  template<class T, class U>
    inline void EvaluateSpline(const T& x      , const T& y     , const T& p
                       ,      T& new_x  ,       U& new_y ,       U& new_y_d
                       ,const double res)
  {
    ResampleArray(x, new_x, res);
    EvaluateSpline(x, y, p, new_x, new_y, new_y_d);
  }


  //-----------------------------------------------------------------------------
  /*!
  * \brief Spline interpolation of given set of data pairs (x,y)
  * 
  * for explanation of the used algorithm, look into:\n
  * 
  *
  * \param x abscissa value data structure (e.g. vector, valarray)
  * \param y ordinate value data structure
  * \param new_x new New abscissa data structure. It must contain the new desired 
  *values already at function calling.
  * \param new_y Interpolated ordinate values.
  * \param new_y_d first derivative of spline interpolation ordinate values
  * \return Updated data structures new_y, new_y_d
  */
  template<class T, class U>
    inline void SplineApproximation(const T& x      , const T& y     
                                   ,const T& new_x  ,       U& new_y 
                                   ,      U& new_y_d)
  {
    T p;
    MathLib::MakeSpline(x,y,p);
    EvaluateSpline(x,y,p,new_x,new_y,new_y_d);
  }


  //-----------------------------------------------------------------------------
  /*!
  * \brief Spline interpolation of given set of data pairs (x,y) (return 
  *second derivatives)
  * 
  * for explanation of the used algorithm, look into (1):\n
  * \n
  * References:\n
  * (1) Sedgewick Robert, "Algorithmen in C++", Addison-Wesley (Deutschland),
  *1995\n
  * 
  *
  * \param x abscissa value data structure (e.g. vector, valarray)
  * \param y ordinate value data structure
  * \param new_x new New abscissa data structure. It must contain the new desired 
  *values already at function calling.
  * \param new_y Interpolated ordinate values.
  * \param new_y_d first derivative of spline interpolation ordinate values
  * \param p second derivatives of interpolation splines at knots
  * \return Updated data structures new_y, new_y_d, p
  */
  template<class T, class U>
    inline void SplineApproximation(const T& x      , const T& y     
                                   ,const T& new_x  ,       U& new_y 
                                   ,      U& new_y_d, T& p)
  {
    MathLib::MakeSpline(x,y,p);
    EvaluateSpline(x,y,p,new_x,new_y,new_y_d);
  }


  //-----------------------------------------------------------------------------
  /*!
  * \brief Spline interpolation of given set of data pairs (x,y)
  * 
  *
  * \param x abscissa value data structure (e.g. vector, valarray)
  * \param y ordinate value data structure
  * \param new_x new New abscissa data structure. It must contain the new desired 
  *values already at function calling.
  * \param new_y Interpolated ordinate values.
  * \return Filled data structure new_y
  */
  template<class T, class U>
    inline void SplineApproximation(const T& x      , const T& y     
                                   ,const T& new_x  ,       U& new_y )
  {
    T p;
    MathLib::MakeSpline(x,y,p);
    EvaluateSpline(x,y,p,new_x,new_y);
  };

}; // namespace MathLib
#endif