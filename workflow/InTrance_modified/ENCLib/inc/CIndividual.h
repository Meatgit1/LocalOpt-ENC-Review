/*!
\file CIndividual.h
\brief A class that provides the functionality for an individual of an
evolutionary algorithm.

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
  andreas.ohndorf@gmx.de


\date

  23.10.2006 (AO)

\b Compiler:

  GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++

\b Copyright:

  (c) 2001-2006 by Bernd Dachwald

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/
#ifndef CINDIVIDUAL_H_INCLUDED
#define CINDIVIDUAL_H_INCLUDED

#ifdef DEBUG
#include <iostream>
#endif

// ENCLib
#include <configENC.h>
#include <CChromosome.h>

class CEA;
class CEAFitness;
class CEAEvalFunc;


/// evolutionary algorithm individual (a population member)
class CIndividual 
{
  friend class CEAOperator;
  friend std::ostream& operator<<(std::ostream& os, const CIndividual& rhs);

  int m_length;                                           ///< chromosome length

  /// size of the initial hypercube (needed for mapping operations)
  double m_initial_hypercube_size;

  /// get initial hypercube size
  double getInitialHypercubeSize(void) const throw();
  
  CEAFitness* m_fitness;                ///< application specific fitness object
  CChromosome m_chromosome;                         ///< individual's chromosome
                
  /// get individual's chromosome
  const CChromosome& getChromosome(void) const throw();

  /// set individual's chromosome
  void setChromosome(const CChromosome& chrom) throw();

  CChromosome m_partial_solution;             ///< individual's partial solution

  /// get individual's partial solution
  const CChromosome& getPartialSolution(void) const throw();

  /// set individual's partial solution
  void setPartialSolution(const CChromosome& ps) throw();
  
  CChromosome m_solution;                            ///< inidividual's solution

  /// set individual's solution
  void setSolution(const CChromosome& sol) throw();

//=============================================================================
// CIndividual CLASS INTERFACE
//=============================================================================
public:
  
  // __constructors__

  /*!
   * \brief standard constructor.
   */
  CIndividual(void);


  /*!
   * \brief Constructor with evolutionary algorithm object.
   *
   *\param ea pointer to CEA class object.
   */
  CIndividual(CEA* ea);


  /*!
   * \brief Constructor with CEA class pointer and init value for all alleles.
   *
   * \param ea pointer to CEA class object
   * \param value Init value to set all the chromosome's alleles with
   */
  CIndividual(CEA* ea,
              // pointer to the evolutionary algorithm object
              const double value);
              // every allele of the individual's chromosome gets this value


  /*!
   * \brief Copy constructor with CEA class object pointer and individual to copy.
   * \param ea pointer to the evolutionary algorithm object
   * \param partial_solution partial solution of the new individual
   */
  CIndividual(CEA* ea, const CIndividual& partial_solution);


  /*!
   * \brief Copy constructor with partial solution and init value for chromosome.
   *
   * \param ea pointer to the evolutionary algorithm object
   * \param allele_value every allele of the individual's chromosome gets this value
   * \param ps partial solution of the new individual
   */
  CIndividual(CEA* ea,
              const double allele_value,
              const CIndividual& ps);


  /*!
   * \brief Constructor to initialize a individual from file.
   *
   * \param ea pointer to the evolutionary algorithm object
   * \param filename file from which individual is to be created
   */
  CIndividual(CEA* ea,
              const char* filename);

  /*!
   * \brief copy constructor
   *
   * \param ind CIndividual class object to initialize the new object with.
   */
  CIndividual(const CIndividual& ind);


  /// assign operator
  const CIndividual& operator=(const CIndividual& ind);


  /// destructor
  ~CIndividual(void);
  
  /*! \brief set allele value at given locus
   *
   *  note: THE CALLING FUNCTION HAS TO UPDATE THE SOLUTION AND TO 
   * RE-EVALUATE THE INDIVIDUAL, IF NECESSARY
   * \param locus chromosome index (of allele to change)
   * \param value new allele value   */
  void setAllele(const int locus, const double value);
  
  /// true, if individual already has a fitness object
  bool hasFitnessObject(void) const throw();

  /// returns inidividual's fitness object pointer
  CEAFitness* getFitnessObject(void) const throw();

  /// set new/initial fitness object
  void setFitnessObject(CEAFitness* fitness_object) throw();

  /*!
   * \brief returns solution embodied by individual as 
   * CChromosome class object pointer
   *
   * The solution is composed of the partial solution and the
   * chromosome. Basically it is nothing more than an addition
   * of every chromosome allele to the allele at the same loci
   * from the partial solution.
   *
   * \return CChromosome class pointer representing the solution. */
  const CChromosome& getSolution(void) const throw();

  CChromosome& getSolution(void) throw();

  /// length of chromosome (=length of solution/partial solution)
  int getLength(void) const throw();

  /*!
   * \brief get certain allele value
   *
   * \param locus the point of interest in chromosome */
  double getAllele(const int locus) const throw();

  /*!
   * \brief save individual's solution to mentioned file
   *
   * \param filename name of the file to save the individuals data to */
  void Save(const char* filename) const throw();

  /*!
   * \brief allele value at given chromosome locus
   *
   * \param locus position of interest in individual's chromosome
   * \param min_value not used
   * \param max_value not used
   */
  double getAllele(const int locus,
                   const double min_value,
                   const double max_value) const;

  /*!
   * \brief Maps current solution allele onto scope defined by a minimial
   *and maximal value
   *
   * Depending on the individuals initial hypercube size (IHCS) the allele value
   * at given position is taken and normed by this initial IHCS. The resulting
   * fractional is multiplied by the scope defined by the two boundary
   * parameters. The result is always in [min_value, max_value].
   *
   * \param locus position of interest in individual's chromosome
   * \param min_value lower boundary of mapping scope
   * \param max_value upper boundaray of mapping scope
   */
  double getMappedAllele(const int locus,
                         const double min_value,
                         const double max_value,
                         const bool constrain = false);

  double getSphereMappedAllele(const int locus_1,
                               const int locus_2,
                               const double min_value,
                               const double max_value) const;

  /*!
   * \brief Maps current solution allele onto angular scope defined by 
   *  minimial and maximal angle values
   *
   * Depending on the individuals initial hypercube size (IHCS) the allele value
   * at given position is taken and normed by this initial IHCS. The resulting
   * fractional is multiplied by the scope defined by the two boundary
   * angle parameters.
   *
   * \param locus position of interest in individual's chromosome
   * \param min_value lower boundary of mapping scope (angle)
   * \param max_value upper boundaray of mapping scope (angle)
   */
  double getCircMappedAllele(const int locus,
                             const double min_value,
                             const double max_value) const;

  double getCircSphereMappedAllele( const int locus_1
                                   ,const int locus_2
                                   ,const double min_value
                                   ,const double max_value) const;

  /*!
   * \brief Combines the CIndividual's chromosome with a possibly stored partial solution.
   *
   *  The result of this method is a combination of the individual's chromosome
   *  and a partial solution whereas the combination usually is a simple addition.
   */
  void UpdateSolution(void);


  /*!
  \brief Sets solution from chromosome and partial solution.

  First the solution of individual will be set by adding every
  allele value of the chromosome to those of the partial solution.
  Then the new solution will also be set as new partial solution
  and finally the chromosome will be cleared.
  */
  void AddChromosomeToPartialSolution(void);


  /*!
   * \brief Modifies one allele at randomly choosen loci
   *
   * A chromosome locus is choosen randomly and its allele value is modified by
   * adding a positive or negative value, that also is choosen randomly
   * but from the scope [-max_muation_value, +max_mutation_value].
   *
   * \param max_mutation_value scope size of which the mutation value will be 
   taken
   */
  void OnePointMutation(const double max_mutation_value);

  void GeneralMutation(const double max_mutation_value,const double probability);

//  void ShiftMutation(const double mutation_sdev);

}; // end class CIndividual

// implementation

//------------------------------------------------------------------------------
inline void CIndividual::setAllele(const int locus, const double value)
{
  m_chromosome.setAllele(locus, value);
} // void setAllele(const int locus, const double value);

//------------------------------------------------------------------------------
inline double CIndividual::getInitialHypercubeSize(void) const throw()
{ 
  return m_initial_hypercube_size; 
} // double getInitialHypercubeSize(void) const throw()

//------------------------------------------------------------------------------
inline const CChromosome& CIndividual::getChromosome(void) const throw()
{ 
  return m_chromosome;
} // const CChromosome& getChromosome(void) const throw()

//------------------------------------------------------------------------------
inline void CIndividual::setChromosome(const CChromosome& chrom) throw()
{ 
  m_chromosome = chrom;
} // void setChromosome(const CChromosome& chrom) throw()

//------------------------------------------------------------------------------
inline const CChromosome& CIndividual::getPartialSolution(void) const throw()
{ 
  return m_partial_solution;
} // const CChromosome& getPartialSolution(void) const throw()

//------------------------------------------------------------------------------
inline void CIndividual::setPartialSolution(const CChromosome& ps) throw()
{ 
  m_partial_solution = ps;
} // void setPartialSolution(const CChromosome& ps) throw()

//------------------------------------------------------------------------------
inline void CIndividual::setSolution(const CChromosome& sol) throw()
{ 
  m_solution = sol;
} // void setSolution(const CChromosome& sol) throw()

//------------------------------------------------------------------------------
inline CIndividual::~CIndividual(void)
{
  SafeDelete(m_fitness);
} // ~CIndividual

//------------------------------------------------------------------------------
inline bool CIndividual::hasFitnessObject(void) const throw()
{
  return m_fitness!=0;
} // bool hasFitnessObject(void) const

//------------------------------------------------------------------------------
inline CEAFitness* CIndividual::getFitnessObject(void) const throw()
{
  return m_fitness;
} // CEAFitness* getFitnessObject(void) const throw()

//------------------------------------------------------------------------------
inline void CIndividual::setFitnessObject(CEAFitness* fitness_object) throw()
{
  if (hasFitnessObject()) SafeDelete(m_fitness);
  m_fitness = fitness_object;
} // void setFitnessObject(CEAFitness* fitness_object) throw()

//------------------------------------------------------------------------------
inline const CChromosome& CIndividual::getSolution(void) const throw()
{ 
  return m_solution; 
} // const CChromosome& getSolution(void) const throw()

//------------------------------------------------------------------------------
inline CChromosome& CIndividual::getSolution(void) throw()
{
  return m_solution;
} // CChromosome& getSolution(void) throw()

//------------------------------------------------------------------------------
inline int CIndividual::getLength(void) const throw()
{ 
  return m_length;
} // int getLength(void) const throw()

//------------------------------------------------------------------------------
inline double CIndividual::getAllele(const int locus) const throw()
{
  return m_chromosome.getAllele(locus);
} // double getAllele(const int locus) const throw()

//------------------------------------------------------------------------------
inline void CIndividual::Save(const char* filename) const throw()
{
  m_solution.Save(filename);
} // void Save(const char* filename) const throw()

std::ostream& operator<<(std::ostream& os, const CIndividual& rhs);

#endif // CINDIVIDUAL_H_INCLUDED
