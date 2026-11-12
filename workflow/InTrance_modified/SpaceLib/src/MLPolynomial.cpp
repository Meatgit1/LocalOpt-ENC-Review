/*!
\file MLPolynomial.cpp
\version 1.7
\brief A class for handling polynomials.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 10.05.2005 (AO)

\b Purpose:
   Class Polynomial is intended to simplifiy the use of
   polynomials, currently only double precision ones.
   It's main purpose has been the approximation of data
   series with cubic polynomials but it may also be used
   in other applications. This file also implements the
   approximation routines to fit polynomials to given
   sets of data. Right now, all the calculation is done
   with double precision.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.

\b Copyright:
   (c) 2005 Bernd Dachwald, Andreas Ohndorf
*/

// C standard library
#include <math.h>

// STL
#include <vector>
#include <sstream>

// MathLib
#include <MLGlobal.h>
#include <MLPolynomial.h>
#include <MLMatrix.h>


namespace MathLib {

  using namespace std;

  //---------------------------------------------------------------------------
  Polynomial::Polynomial(int degree, double* coeff): m_coeff(degree+1) {
    if (coeff != 0)
    {
      for (int i=0;i<=degree;++i) m_coeff(i) = coeff[i];
    }
    else
    {
      m_coeff.Set(Vector());
    }
  }

  //---------------------------------------------------------------------------
  Polynomial& Polynomial::operator=(const Polynomial& right)
  {
    m_coeff.Set(right.m_coeff);
    return *this;
  }

  //---------------------------------------------------------------------------
  double Polynomial::getValue(const double x) const
  {
    double return_val = 0;
    int n_coeff = m_coeff.getDimension()-1;
    if (n_coeff>0)
    {
      for (int i=0;i<=n_coeff;++i) return_val+=m_coeff(i)*pow(x,double(i));
    }
    return return_val;
  }

  //---------------------------------------------------------------------------
  double Polynomial::firstDerivative(const double x) const
  {
    double return_val = 0.0;
    int n_coeff = m_coeff.getDimension()-1;
    if (n_coeff >= 1)
    {
      for (int i=1;i <= n_coeff;++i) 
        return_val += (double(i) * m_coeff(i) * pow(x,double(i-1)));
    }
    return return_val;
  }


  //---------------------------------------------------------------------------
  double Polynomial::secondDerivative(const double x) const
  {
    double return_val = 0.0;
    int n_coeff = m_coeff.getDimension()-1;
    if (n_coeff >= 2)
    {
      for (int i=2;i <= n_coeff;++i) 
        return_val += (double(i) * double(i-1)* m_coeff(i) * pow(x,double(i-2)));
    }
    return return_val;
  } // Polynomial::secondDerivative(const double x)



  //---------------------------------------------------------------------------
  template<class T>
    void NormalizeSubArray(
           T& d               // data to normalize
    ,const T& iv              // independent variable
    ,const unsigned int begin // first index of scope to normalize
    ,const unsigned int end   // last index of scope to normalize
    ,double& C0               // normalizing factor 1 ( =(max d(i) )-(min d(i) ) )
    ,double& C1               // normalizing factor 2 ( = C0/(iv[end]-iv[begin]))
    )
  {
    if (begin > end) return;
    if (d.size()!=iv.size()) return;
    if (end >= d.size()) return;

    double max_val(0);
    double min_val(0);

    for (unsigned int i = begin; i<=end;++i)
    {
      max_val=max(max(d[i],max_val));
      min_val=min(max(d[i],max_val));
    }
    C0 = max_val - min_val;
    C1 = C0/(iv[end]-iv[begin]);
    for (unsigned int i = begin; i<=end;++i) d[i] = (d[i]-min_val)/C0;
  }




  //-----------------------------------------------------------------------------
  void SetUpApproximationArrays(
    const std::valarray<double>& data
    ,      std::valarray<unsigned int>& segm
    ,      std::valarray<unsigned int>& degrees
    , const unsigned int segment_size
    , const unsigned int poly_degree)
  {
    unsigned int curr_end(0);
    const uint end_idx(data.size()-1);
    uint segments(0);
    while (curr_end < end_idx)
    {
      ++segments;
      curr_end+=(segment_size-1);
    }
    segm.resize(segments);
    degrees.resize(segments);
    if (segment_size > end_idx) curr_end = end_idx;
    else curr_end = segment_size-1;
    for (uint i = 0; i < segments;++i)
    {
      segm[i] = curr_end>end_idx?end_idx:curr_end;
      degrees[i] = poly_degree;
      curr_end+=segment_size-1;
    }
  }



  //-----------------------------------------------------------------------------
  void SetUpApproximationArrays(
    unsigned int                      segs
    ,unsigned int*                    seg_end_indices
    ,unsigned int*                    seg_poly_deg
    ,         valarray<unsigned int>& segm
    ,         valarray<unsigned int>& degrees 
    )
  {
    segm.resize(segs);
    degrees.resize(segs);
    for (uint i=0; i< segs;++i)
    {
      segm[i] = seg_end_indices[i];
      degrees[i] = seg_poly_deg[i];
    }
  }


  //----------------------------------------------------------------------------
  // system matrix of one segment (helper function for CalcPolyApprox)
  template<class T>
    Matrix ComputeSystemMatrix(
    const T& input
    ,const T& val
    ,const uint b
    ,const uint e
    ,const uint degree
    ,T& value_vec)
  {
    if ( (input.size()>0) && (input.size()==val.size()) && (int(e-b)>0) && (degree>1) )
    {
      vector<T> f;f.resize(degree+1);
      for (uint i=0; i<=degree; ++i) f[i]=CalcFuncVec(input, b, e, i);
      Matrix M(degree+1, degree+1);
      T v(value_vec);value_vec.resize(degree+1);
      for (uint i=0; i<=degree; ++i)
      {
        M(i,i) = (f[i]*f[i]).sum();     // calculate "vector" - product
        value_vec[i]=(v*f[i]).sum();    // calculate right hand side element
        for (uint j=0; j<i; ++j) 
        {
          M(i,j) = M(j,i) = (f[i]*f[j]).sum();
          if (fabs(M(i,j))<1.0E-10) M(i,j)=M(j,i)=0.0;
        }
      }
      return M;
    }
    return Matrix();
  } // ComputeSystemMatrix


  //---------------------------------------------------------------------------
  /*! \brief Calculates transition condition between two segments.
   *
   * \param d_l  left side polynomial degree
   * \param d_r right side polynomial degree
   * \param l_bound independed variable at end of left side segment
   * \param mid independed variable at segment boundary
   * \param r_bound independed variable at end of right side segment
   */
  template<class T>
   inline Matrix ComputeTransitionCondition(
      const uint d_l
    , const uint d_r
    , const T    l_bound
    , const T    mid
    , const T    r_bound
    )
  {
    Matrix R(3, d_l+d_r+5);
    // continuity
    for (uint i=0;i<=d_l;++i) R(0,i)       = 1.0;
    for (uint i=0;i<=d_r;++i) R(0,d_l+4+i) = ((i%2)==0?-1.0:1.0);
    // smoothness (first derivative)
    const T f1(1.0/(mid-l_bound)), f2(1.0/(r_bound-mid));
    for (uint i=1;i<=d_l;++i) R(1,i) = 1.0*i*f1;
    for (uint i=1;i<=d_r;++i) R(1,d_l+4+i) = ((i%2)==0?1.0*i*f2:-1.0*i*f2);
    // second derivative
    const T f11(f1*f1), f22(f2*f2);
    for (uint i=2;i<=d_l;++i) R(2,i) = i*(i-1)*f11;
    for (uint i=2;i<=d_r;++i) R(2,d_l+4+i) = ((i%2)==0?i*(1.0-double(i))*f22:i*(i-1)*f22);
    return R;
  }


  //-----------------------------------------------------------------------------
  // calculation of polynomials for piecewise concatenation to obtain a
  // polynomial approximation of a given input data set (x,y pairs)
  void CalcPolyApprox(
    const valarray<double>& iv    ///< independend variable
   ,const valarray<double>&  y    ///< data
   ,const valarray<uint>&    s    ///< segmentation
   ,const valarray<uint>&    d    ///< degree of segment's polynomials
   ,vector<Polynomial>&     p     ///< polynomials
   ,valarray<double>&       ny    ///< new approximation of data 
   )
  {
    if (iv.size()>0 && iv.size()==y.size() && s.size()>0 && s.size()==d.size())
    {
      p.resize(s.size());
      valarray<uint> ss;            // new segmentation
      CompModSegVec(s, ss);         // alter segmentation vector
      // prepare set of new, input vectors whose elements are mapped
      // onto (-1,1) interval
      vector< valarray<double> > v_iv(ss.size()-1); // vecs of modified input data
      vector< valarray<double> > v_yv(ss.size()-1); // vecs of right sides
      uint row(0),col(0);
      const uint dim(d.sum()+d.size()+3*(ss.size()-2)); // system matrix dimension
      Matrix A(dim,dim);                                // system matrix (whole)
      Vector b(dim);
      for (uint i=0;i<v_iv.size();++i)                  // segment loop
      {
        CalcModInput(iv,ss[i],ss[i+1],v_iv[i]);
        CalcDataVecs( y,ss[i],ss[i+1],v_yv[i]);
        A.SetSubMatrix(row,col, ComputeSystemMatrix(v_iv[i],v_yv[i], d[i], v_yv[i]));
        SetVector(v_yv[i], b, col);
        if (i<v_iv.size()-1)
          A.SetSubMatrix(row+=(d[i]+1),col, 
            ComputeTransitionCondition(d[i], d[i+1], iv[ss[i]], iv[ss[i+1]], iv[ss[i+2]]));
        row+=3;col+=(d[i]+4);
      }
      Vector h(dim);
      LU_Decomp(A, h, true);
      LU_BackSub(A, h, b, true);

      ny.resize(y.size());
      uint j(0),m(0);
      for (uint i=0;i<v_iv.size();++i)
      {
        const double f1(2.0/(iv[ss[i+1]]-iv[ss[i]]));
        p[i]=Polynomial(b, j, j+d[i]);
        for (uint k=0;k<v_iv[i].size();++k)
        {
          ny[m+k]=p[i].getValue((v_iv[i])[k]);
        }
        m+=v_iv[i].size()-1;
        j+=(d[i]+4);
      }
    } else Warning("CalcPolyApprox(...): input invalid!\n");
  } // CalcPolyApprox



  //-----------------------------------------------------------------------------
  // data values approximation at given knot points
  void CalcPolyApproxData(
    const vector<Polynomial>&     polys  ///< vector of approximation polynomials
    ,const std::valarray<double>& t      ///< original time resolution
    ,const std::valarray<uint>&   s      ///< segmentation array
    ,const double             rf         ///< resolution factor
    ,      std::valarray<double>& nt     ///< new time resolution
    ,      std::valarray<double>& nd     ///< new data array
    )
  {
    if (s.size() != polys.size() || t.size() == 0) return;
    valarray<uint> ss;            // new segmentation
    CompModSegVec(s, ss);         // alter segmentation vector

    nd.resize(int(ceil(rf*(t.size()-1)+1)));
    nt.resize(nd.size());

    const double delta_t(t[t.size()-1]-t[0]); // input span
    const double t_inc(delta_t/(nt.size()-1));
    unsigned int idx(0);                      // index in ivar array
    unsigned int segment(0);                  // segment index
    Polynomial polynomial(polys[segment]);    // polynomial of current segment

    double start_ivar(t[0]);
    double last_start_ivar(start_ivar);
    double output_ivar(0.0);
    double calc_ivar(-1.0);
    double t_fac(2.0/(t[s[0]]-t[0]));

    do
    {
      // get approximation value
      output_ivar = start_ivar+idx*t_inc;

      calc_ivar = (output_ivar-last_start_ivar)*t_fac-1.0;
      nt[idx]   = output_ivar;
      nd[idx]   = polynomial.getValue(calc_ivar);
      ++idx;

      if (segment+1 < s.size() && (output_ivar+t_inc) > t[s[segment]])
      {
        last_start_ivar=t[s[segment]];
        t_fac = (2.0/(t[s[segment+1]]-t[s[segment]]));
        polynomial = polys[++segment];      // set new polynomial
      }
    }
    while (idx < nt.size());   // end loop after last point
  } // void CalcPolyApproxData(...)


  //-----------------------------------------------------------------------------
  // first derivatives of approximation at given knot points
  void CalcPolyApproxFirstDerivative(
    const vector<Polynomial>&     polys  ///< vector of approximation polynomials
    ,const std::valarray<double>& t      ///< original time resolution
    ,const std::valarray<uint>&   s      ///< segmentation array
    ,const double             rf         ///< resolution factor
    ,      std::valarray<double>& nt     ///< new time resolution
    ,      std::valarray<double>& dd     ///< new derivatives array
    )
  {
    if (s.size() != polys.size() || t.size() == 0) return;
    valarray<uint> ss;            // new segmentation
    CompModSegVec(s, ss);         // alter segmentation vector

    dd.resize(int(ceil(rf*(t.size()-1)+1)));
    nt.resize(dd.size());

    const double delta_t(t[t.size()-1]-t[0]); // input span
    const double t_inc(delta_t/(nt.size()-1));
    unsigned int idx(0);                      // index in ivar array
    unsigned int segment(0);                  // segment index
    Polynomial polynomial(polys[segment]);    // polynomial of current segment

    double start_ivar(t[0]);
    double last_start_ivar(start_ivar);
    double output_ivar(0.0);
    double calc_ivar(-1.0);
    double t_fac(2.0/(t[s[0]]-t[0]));

    do
    {
      // get approximation value
      output_ivar = start_ivar+idx*t_inc;

      calc_ivar = (output_ivar-last_start_ivar)*t_fac-1.0;
      nt[idx]   = output_ivar;
      dd[idx]   = t_fac*polynomial.firstDerivative(calc_ivar);
      ++idx;

      if (segment+1 < s.size() && (output_ivar+t_inc) > t[s[segment]])
      {
        last_start_ivar=t[s[segment]];
        t_fac = (2.0/(t[s[segment+1]]-t[s[segment]]));
        polynomial = polys[++segment];      // set new polynomial
      }
    }
    while (idx < nt.size());   // end loop after last point
  } // void CalcPolyApproxFirstDerivative(...)



  //-----------------------------------------------------------------------------
  // second derivatives of approximation at given knot points
  void CalcPolyApproxSecondDerivative(
    const vector<Polynomial>&     polys  ///< vector of approximation polynomials
    ,const std::valarray<double>& t      ///< original time variable array
    ,const std::valarray<uint>&   s      ///< segmentation array
    ,const double                 rf     ///< resolution factor
    ,      std::valarray<double>& nt     ///< new time resolution
    ,      std::valarray<double>& dd     ///< new derivatives array
    )
  {
    if (s.size() != polys.size() || t.size() == 0) return;
    valarray<uint> ss;            // new segmentation
    CompModSegVec(s, ss);         // alter segmentation vector

    dd.resize(int(ceil(rf*(t.size()-1)+1)));
    nt.resize(dd.size());

    const double delta_t(t[t.size()-1]-t[0]); // input span
    const double t_inc(delta_t/(nt.size()-1));
    unsigned int idx(0);                      // index in ivar array
    unsigned int segment(0);                  // segment index
    Polynomial polynomial(polys[segment]);    // polynomial of current segment

    double start_ivar(t[0]);
    double last_start_ivar(start_ivar);
    double output_ivar(0.0);
    double calc_ivar(-1.0);
    double t_fac(2.0/(t[s[0]]-t[0]));
    const double t_fac2(t_fac*t_fac);

    do
    {
      // get approximation value
      output_ivar = start_ivar+idx*t_inc;

      calc_ivar = (output_ivar-last_start_ivar)*t_fac-1.0;
      nt[idx]   = output_ivar;
      dd[idx]   = t_fac2*polynomial.secondDerivative(calc_ivar);
      ++idx;

      if (segment+1 < s.size() && (output_ivar+t_inc) > t[s[segment]])
      {
        last_start_ivar=t[s[segment]];
        t_fac = (2.0/(t[s[segment+1]]-t[s[segment]]));
        polynomial = polys[++segment];      // set new polynomial
      }
    }
    while (idx < nt.size());   // end loop after last point
  } // void CalcPolyApproxSecondDerivative(...)



  //---------------------------------------------------------------------------
  // private template implementation for calculating a spline interpolation
  // see also:
  // (1) Sedgewick Robert, "Algorithmen in C++", Addison-Wesley (Deutschland),1995
  template<class T>
    inline void MakeSpline_local(const T& x, const T& y, T& p)
  {
    if (x.size()==y.size() || x.size()>0) {
      const uint N(x.size());p.resize(N);
      valarray<double> d(N-2), u(N-1), w(N-2);      
      for(uint i=1;i<N-1;++i) d[i-1] = 2.0*(x[i+1]-x[i-1]);
      for(uint i=0;i<N-1;++i) u[i] = x[i+1]-x[i];
      for(uint i=1;i<N-1;++i) w[i-1]=6.0*((y[i+1]-y[i])/u[i]-(y[i]-y[i-1])/u[i-1]);
      p[0] = p[N-1] = 0.0;
      for (uint i=1; i<(N-2);++i)
      {
        w[i] = w[i]-w[i-1]*u[i]/d[i-1];
        d[i] = d[i]-u[i]*u[i]/d[i-1];
      }
      for (uint i=N-2;i>0;--i)
      {
        p[i] = (w[i-1]-u[i]*p[i+1])/d[i-1];
      }
    }
  } // MakeSpline(const T&, const T&, T&)


  //---------------------------------------------------------------------------
  // MakeSpline implementation for std::vector
  void MakeSpline(
    const std::vector<double>& x
    ,const std::vector<double>& y
    ,      std::vector<double>& p)
  {
    MakeSpline_local(x,y,p);
  }


  //---------------------------------------------------------------------------
  // MakeSpline implementation for std::valarrays
  void MakeSpline(
    const std::valarray<double>& x
    ,const std::valarray<double>& y
    ,      std::valarray<double>& p)
  {
    MakeSpline_local(x,y,p);
  }

} // namespace MathLib