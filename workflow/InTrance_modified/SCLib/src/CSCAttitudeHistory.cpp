/*!
\file CSCAttitudeHistory.cpp
\brief A class that provides the functionality for handling the attitude history of a spacecraft.

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

\date

  26.06.2005

\b Compiler:

  Microsoft Visual C++ .NET

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
#include <string>

//SCLib
#include <configSC.h>

// SpaceLib
#include <SpaceLib.h>

// SCLib
#include <CSCAttitudeHistory.h>


using namespace std;
using namespace Global;


// ----------------------------------------------------------------------------
CSCAttitudeHistory::CSCAttitudeHistory(const std::string& filename)
// ----------------------------------------------------------------------------
// public
: m_n(0), m_v(NULL)
{
  const string err_msg1 = "CSCAttitudeHistory::CSCAttitudeHistory(const \
std::string& filename): ";
  const string err_msg2 = "no input file given!";

  if (filename.size() == 0) throw std::logic_error(err_msg1+err_msg2);

  ifstream in(filename.c_str());

  if (in.is_open() == false) throw std::invalid_argument(err_msg1+filename+
    string(" not found!"));

  m_n = static_cast<int>(ReadNumber(in));

  m_v = new CSCAttitude*[m_n];
        Assure(m_v, "CSCAttitudeHistory::CSCAttitudeHistory()");

  double mjd, clock, cone;

  for (int i = 0; i < m_n; i++)
  {
    in >> mjd;
    in >> clock; 
    in >> cone;
    m_v[i] = new CSCAttitude(clock, cone);
             Assure(m_v[i], "CSCAttitudeHistory::CSCAttitudeHistory()");
    m_v[i]->setMJD(mjd);
  }

};

// ----------------------------------------------------------------------------
CSCAttitudeHistory::CSCAttitudeHistory(const CIndividual& solution,
                                       const double       mjd0,
                                       const double       dmjd)
// ----------------------------------------------------------------------------
// public
: m_n(0), m_v(NULL)
{
  m_n = solution.getLength()/2;

  m_v = new CSCAttitude*[m_n];
        Assure(m_v, "CSCAttitudeHistory::CSCAttitudeHistory()");

  double clock, cone;

  for (int i = 0, j = 0; i < m_n; i++, j += 2)
  {
    clock = solution.getAllele(j);
    cone  = solution.getAllele(j+1);
    m_v[i] = new CSCAttitude(clock, cone);
             Assure(m_v[i], "CSCAttitudeHistory::CSCAttitudeHistory()");
    m_v[i]->setMJD(mjd0 + j*dmjd);
  }

}; // end CSCAttitudeHistory::CSCAttitudeHistory


// ----------------------------------------------------------------------------
CSCAttitudeHistory::~CSCAttitudeHistory()
// ----------------------------------------------------------------------------
// public
{
  if (m_v!=0)
  {
    for (int i = 0; i < m_n; ++i) SafeDelete(m_v[i]);
    SafeDeleteArray(m_v);
  }
}; // end CSCAttitudeHistory::~CSCAttitudeHistory


// ----------------------------------------------------------------------------
int CSCAttitudeHistory::getDimension(void) const
// ----------------------------------------------------------------------------
// public
{
  return m_n;

}; // end CSCAttitudeHistory::getDimension


// ----------------------------------------------------------------------------
double CSCAttitudeHistory::getMJD(int element) const
// ----------------------------------------------------------------------------
// public
{
#ifdef development_version
  if(element >= m_n) FatalError("CSCAttitudeHistory::getMJD(): index out of range");
#endif

  return m_v[element]->getMJD();

}; // end CSCAttitudeHistory::getMJD


// ----------------------------------------------------------------------------
double CSCAttitudeHistory::getClockAngle(int element) const
// ----------------------------------------------------------------------------
// public
{
#ifdef development_version
  if(element >= m_n) FatalError("CSCAttitudeHistory::getClockAngle(): index out of range");
#endif

  return m_v[element]->getClockAngle();

}; // end CSCAttitudeHistory::getClockAngle


// ----------------------------------------------------------------------------
double CSCAttitudeHistory::getConeAngle(int element) const
// ----------------------------------------------------------------------------
// public
{
#ifdef development_version
  if(element >= m_n) FatalError("CSCAttitudeHistory::getConeAngle(): index out of range");
#endif

  return m_v[element]->getConeAngle();

}; // end CSCAttitudeHistory::getConeAngle
