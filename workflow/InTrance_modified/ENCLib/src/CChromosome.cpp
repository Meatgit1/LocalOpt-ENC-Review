/*!
\file CChromosome.cpp
\brief A class that provides the functionality for a chromosome of an evolutionary algorithm.

\version

  1.7

\author

  Bernd Dachwald\n
  German Aerospace Center (DLR)\n
  Institute of Space Simulation\n
  Linder Hoehe\n
  51147 Cologne\n
  Germany\n
  bernd.dachwald@dlr.de

\author

  Andreas Ohndorf (AO)\n
  German Aerospace Center (DLR)\n
  Institute of Space Simulation\n
  Linder Hoehe\n
  51147 Cologne\n
  Germany\n
  andreas.ohndorf@dlr.de


\date

  29.08.2005

\b Compiler:

  GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++

\b Copyright:

  (c) 2001-2004 by Bernd Dachwald

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/
// STL
#include <iomanip>

// SpaceLib
#include <Global.h>
#include <MLConst.h>
#include <MLGlobal.h>

// ENCLib
#include <CChromosome.h>


using namespace std;
using namespace Global;
using namespace MathLib;


//-----------------------------------------------------------------------------
void CChromosome::setChromosomeRandom(const double min_value, const double max_value)
{
  int i;
  if ( min_value != max_value )
    for (i = 0; i < n; ++i) v[i] = randInterval(min_value, max_value);
  else for (i=0; i<n; ++i) v[i] = min_value;
};


/* function not used at the moment
double CChromosome::Norm(const VectorNorm norm)
{
  double result = 0.0;
  int i;

  switch (norm)
  {
    case MAXIMUM_SUM_NORM: // (1-norm)

      for (i = 0; i < n; i++)
        result += fabs(v[i]);
      result /= static_cast<double>(n);
      break;

    case EUCLIDEAN_NORM: // (2-norm)

      for (i = 0; i < n; i++)
        result += pow(v[i], 2);
      result = sqrt(result/static_cast<double>(n));
      break;

    case MAXIMUM_ELEMENT_NORM: // (infinity norm)

      for (i = 0; i < n; i++)
        result = (fabs(v[i]) > result) ? fabs(v[i]) : result;
      break;

  } // end switch

  return result;
};*/


//------------------------------------------------------------------------------
void CChromosome::Load(const char* filename)
{
  ifstream in(filename);
  Assure(in, filename);

  int n_loaded = static_cast<int>(ReadNumber(in));

#ifdef development_version
  if (n_loaded != n) FatalError("CChromosome::Load(): loaded chromosome has wrong size");
#endif

  for (int i = 0; i < n; ++i) v[i] = ReadNumber(in);
};


//------------------------------------------------------------------------------
double CChromosome::getAllele(const int locus) const
{
#ifdef development_version
  if (locus < 0 || locus >= n) 
    FatalError("CChromosome::getAllele(): locus out of range");
#endif

  return v[locus];
};


//------------------------------------------------------------------------------
void CChromosome::setAllele(const int locus, const double value)
{
#ifdef development_version
  if (locus < 0 || locus >= n) 
    FatalError("CChromosome::setAllele(): locus out of range");
#endif

  v[locus] = value;
};


//------------------------------------------------------------------------------
double CChromosome::getMappedAllele( const int locus
                                    ,const double min_value
                                    ,const double max_value
                                    ,const double orig_range
                                    ,const bool constrain)
{
#ifdef development_version
  if (locus < 0 || locus >= n) 
    FatalError("CChromosome::getMappedAllele(): locus out of range");
#endif

  double _min = min_value;
  double _max = max_value;

  if (_min > _max) swap(_min, _max);

  const double mid    = (_min + _max) / 2.0;
  const double delta  = _max - _min;

  if (constrain == true) 
  {
    v[locus] = bounded(v[locus], -orig_range, orig_range);
  }

  const double offset = v[locus] / orig_range * delta / 2.0;

  return bounded(mid+offset, _min, _max);
};


//------------------------------------------------------------------------------
double CChromosome::getSphereMappedAllele( const int locus_1
                                          ,const int locus_2
                                          ,const double min_value
                                          ,const double max_value
                                          ,const double orig_range) const
{
#ifdef development_version
  if (locus_1 < 0 || locus_1 >= n)
    FatalError("CChromosome::getSphereMappedAllele(): locus_1 out of range");
  if (locus_2 < 0 || locus_2 >= n)
    FatalError("CChromosome::getSphereMappedAllele(): locus_2 out of range");
#endif

  if (min_value < max_value)
  {
    return bounded(min_value+(max_value-min_value)*
      (1.0-sqrt( pow(v[locus_1]/orig_range,int(2)) + pow(v[locus_2]/orig_range,int(2)) ) ),min_value,max_value);
  } else
  {
    return bounded(max_value+(min_value-max_value)*
      (1.0-sqrt( pow(v[locus_1]/orig_range,int(2)) + pow(v[locus_2]/orig_range,int(2)) ) ),max_value,min_value);
  }
}


//------------------------------------------------------------------------------
double CChromosome::getCircMappedAllele(const int locus, 
                                        const double min_value, 
                                        const double max_value,
                                        const double orig_range) const
{
#ifdef development_version
  if (locus < 0 || locus >= n) FatalError("CChromosome::getCircMappedAllele(): locus out of range");
#endif

  if(min_value <= -PI || min_value >= PI) FatalError("CChromosome::getCircMappedAllele(): min. value out of range");
  if(max_value <= -PI || max_value >= PI) FatalError("CChromosome::getCircMappedAllele(): max. value out of range");

  double _min = min_value;
  double _max = max_value;

  const bool center180 = (_min < _max) ? false : true;

  double mid, delta;

  if (!center180)
  {
    mid    = (_min + _max) / 2.0;
    delta  = _max - _min;
  }
  else
  {
    mid    = ((_min + _max) / 2.0) + PI;
             AssertAngularRange(mid);
    delta  = TWOPI - (_max - _min);
  }
  const double offset = v[locus] / orig_range * delta / 2.0;

  return bounded(mid+offset, _min, _max);
};

//------------------------------------------------------------------------------
double CChromosome::getCircSphereMappedAllele(const int locus_1
                                 ,const int locus_2
                                 ,const double min_value
                                 ,const double max_value
                                 ,const double orig_range) const
{
#ifdef development_version
  if (locus_1 < 0 || locus_1 >= n)
    FatalError("CChromosome::getSphereMappedAllele(): locus_1 out of range");
  if (locus_2 < 0 || locus_2 >= n)
    FatalError("CChromosome::getSphereMappedAllele(): locus_2 out of range");
#endif

  double _min = min_value;
  double _max = max_value;

  const bool center180 = (_min < _max) ? false : true;

  double mid, delta;

  if (!center180)
  {
    mid    = (_min + _max) / 2.0;
    delta  = _max - _min;
  }
  else
  {
    mid    = ((_min + _max) / 2.0) + PI;
             AssertAngularRange(mid);
    delta  = TWOPI - (_max - _min);
  }
  const double offset = (1.0-sqrt(pow(v[locus_1]/orig_range,int(2))+pow(v[locus_2]/orig_range,int(2)))) * delta / 2.0;

  return bounded(mid+offset, _min, _max);

}


//------------------------------------------------------------------------------
void CChromosome::OnePointMutation(const double max_mutation_value)
{
  const int locus = randIntInterval(0, n-1);

  v[locus] = randSymInterval(max_mutation_value);
};


//------------------------------------------------------------------------------
void CChromosome::GeneralMutation(const double max_mutation_value, const double probability)
{
  for (int i=0;i<n;++i)
    if (Flip(probability)) v[i] = randSymInterval(max_mutation_value);
}

/* function not used at the moment
void CChromosome::ShiftMutation(const double mutation_sdev)
{
  for (int i = 0; i < n; i++)
    v[i] += randNormal(0.0, 0.01 * mutation_sdev);
};*/


//------------------------------------------------------------------------------
void CChromosome::Save(const char* filename) const
{
  ofstream out(filename);
  Assure(out, filename);

  out << n << endl;
  out << setprecision(16);

  for (int i = 0; i < n; ++i)
  {
    out << v[i];
    if (i != n-1) out << endl;
  }
};
