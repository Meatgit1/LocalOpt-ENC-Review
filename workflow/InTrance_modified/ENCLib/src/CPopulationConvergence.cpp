/*!
\file CPopulationConvergence.cpp
\brief An auxiliary class to handle convergence issues for a population.

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

  05.05.2005 (AO)

\b Compiler:

  GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++

\b Copyright:

  (c) 2001-2005 by Bernd Dachwald, Andreas Ohndorf

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/
// standard C library
#include <math.h>

// SpaceLib
#include <Global.h>

// ENCLib
#include <CPopulationConvergence.h>


using namespace std;
using namespace Global;


//-----------------------------------------------------------------------------
CPopulationConvergence::CPopulationConvergence(
  const long dim
 ,const double p_limit_fbc_not_met
 ,const double p_limit_fbc_met
 ,const double p_fbc_met )
: m_dim(dim)
 ,m_conv_limit_fbc_not_met(p_limit_fbc_not_met)
 ,m_conv_limit_fbc_met(p_limit_fbc_met)
 ,m_fbc_met(p_fbc_met) 
{

   m_fitness = new double[dim];
              Assure(m_fitness, "CPopulationConvergence::CPopulationConvergence()");
   Reset();
};


//-----------------------------------------------------------------------------
CPopulationConvergence::~CPopulationConvergence() {
  if (m_fitness)
  {
    delete []m_fitness; // delete the array
    m_fitness = NULL;
  } else FatalError("CPopulationConvergence::~CPopulationConvergence()");
};


//-----------------------------------------------------------------------------
void CPopulationConvergence::Reset(void) {
  for (int i = 0; i < m_dim; ++i) m_fitness[i] = -HUGE_VAL;
};


//-----------------------------------------------------------------------------
double CPopulationConvergence::getFitness(const long reproduction) const {

  const long lastindex = (reproduction % m_dim) + 1;
  const long index     = (lastindex == m_dim) ? 0 : lastindex;

  if (index < 0 || index >= m_dim) 
    FatalError("CPopulationConvergence::getFitness: index out of range");

  return m_fitness[index];
};


//-----------------------------------------------------------------------------
double CPopulationConvergence::getImprovement(const long reproduction
                                             ,const double fitness) const {

  const double lastfitness = getFitness(reproduction);
  const double improvement = (lastfitness == -HUGE_VAL) ? 1.0 :
  (fitness-lastfitness)/fabs(lastfitness);
  return improvement;
};


//-----------------------------------------------------------------------------
bool CPopulationConvergence::Convergence(const long reproduction
                                        ,const double fitness) const {
  //double f = getFitness(reproduction);
  //double required_percentage;

  //if (f > 1.0) required_percentage = m_conv_limit_fbc_met;
  //else required_percentage = m_conv_limit_fbc_not_met;

  //if (getImprovement(reproduction, fitness) > required_percentage) return false;
  //else return true;
  if (getFitness(reproduction) > m_fbc_met) // FBC met
  {
    if (getImprovement(reproduction, fitness) > m_conv_limit_fbc_met) 
      return false;
    else return true;
  }
  else  // FBC not met
  {
    if (getImprovement(reproduction, fitness) > m_conv_limit_fbc_not_met) 
      return false;
    else return true;
  }
}


//-----------------------------------------------------------------------------
void CPopulationConvergence::setFitness(const long reproduction
                                       ,const double fitness) {

  const long index = reproduction % m_dim;
  if (index < 0 || index >= m_dim) FatalError("setFitness: index out of range");
  m_fitness[index] = fitness;
};
