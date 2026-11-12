/*!
\file CEAFitness.cpp
\brief An auxiliarary class to measure the fitness of an individual.

\version

  1.7

\author

  Bernd Dachwald (BD)\n
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

  18-08.2005 (AO)

\b Compiler:

  GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++

\b Copyright:

  (c) 2001-2004 by Bernd Dachwald, Andreas Ohndorf

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/

#include <Global.h>

#include <CEAFitness.h>

// ----------------------------------------------------------------------------
double CEAFitness::getFitness(const int index) const {
#ifdef development_version
  if(index < 0 || index >= m_fitness.getDimension()) FatalError("CEAFitness::getFitness(const int index): index out of range");
#endif

  return m_fitness.getComponent(index);

}; // end CEAFitness::getFitness


// ----------------------------------------------------------------------------
void CEAFitness::setFitness(const int index, const double value) {
#ifdef development_version
  if(index < 0 || index >= m_fitness.getDimension()) FatalError("CEAFitness::CheckIndex(): index out of range");
#endif

  m_fitness.setComponent(index, value);

}; // end CEAFitness::setFitness
