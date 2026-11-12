/*!
\file CPopulationConvergence.h
\brief An auxiliary class to handle convergence issues for a population.

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

  05.05.2004 (AO)

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
#ifndef CPOPULATIONCONVERGENCE_H_INCLUDED
#define CPOPULATIONCONVERGENCE_H_INCLUDED


// population convergence class declaration
// ----------------------------------------------------------------------------
class CPopulationConvergence
{
  long    m_dim;
  double  m_conv_limit_fbc_not_met; // convergence limit if final boundary condition is not met
  double  m_conv_limit_fbc_met;     // convergence limit if final boundary condition is met
  double  m_fbc_met;                 // fitness(fbc not met) < m_fbc_met <= fitness(fbc met)
  double* m_fitness;


//=============================================================================
// CPopulationConvergence CLASS INTERFACE
//=============================================================================
public:

  // __constructor__

  CPopulationConvergence(
    long dim
   ,const double p_limit_fbc_not_met
   ,const double p_limit_fbc_met
   ,const double p_fbc_met = 1.0  // this is the former hardwired boundary
  );

  // purpose: create population convergence object and initialize

  // __destructor__

  ~CPopulationConvergence();

  // __functionality__

  void Reset(void);

  double getFitness(const long reproduction) const;

  double getImprovement(const long reproduction,
                        const double fitness) const;

  /*!
   * \todo BD: The decision of convergence in this function depends
   * on the last fitness in convergence queue. If it is higher than
   * a certain boundary (final boundary condition met) an other
   * convergence limit is taken than if it is lower (has not met
   * final boundary conditions). Right now it uses a "hardwired"
   * boundary fitness value of 1.0 to check for final BC. Though
   * it is possible to provide  such a value to CEA class objects,
   * I suggest to use it also here, which means to create a new member
   * field named "m_fbc_met", provide means to set/get it and also
   * include it in constructor parameter lists. Finally this value
   * should be used in the "Convergence"-method.
   */
  bool Convergence(const long reproduction,
                   const double fitness) const;

  void setFitness(const long reproduction,
                  const double fitness);

}; // end class CPopulationConvergence


#endif // CPOPULATIONCONVERGENCE_H_INCLUDED
