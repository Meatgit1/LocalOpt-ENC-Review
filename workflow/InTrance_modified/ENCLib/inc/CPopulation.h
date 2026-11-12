/*!
\file CPopulation.h
\brief A class that provides the functionality for a population of an
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

  31.10.2006 (AO)

\b Compiler:

  GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++

\b Copyright:

  (c) 2001-2006 by Bernd Dachwald, Andreas Ohndorf

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/
#ifndef CPOPULATION_H_INCLUDED
#define CPOPULATION_H_INCLUDED

// STL
#include <string>
#include <sstream>
#include <cassert>

// ENCLib
#include <configENC.h>

// classes forward declarations
class CPopulationConvergence;
class CIndividual;
class CEA;


// function declarations

/*!
 *  \brief User provided function for control output purposes.
 *
 *  \param ind
 *  \param improvement
 *  \param nor Number of reproductions.
 *  \param rep Name of report file.
 *  \param data Pointer to user provided data (if needed).
 */
typedef void (*PrintReproductionFunction) (       CIndividual* ind
                                           ,const double improvement
                                           ,const int n_epc_repro
                                           ,const int n_abs_repro
                                           ,void* data);


/// Sets counter of absolute number of reproductions to zero.
void ResetReproCtr(void);

/// population object class (population of individuals)
class CPopulation
{
  friend class CEAOperator;

  /// pointer to the calling evolutionary algorithm
  CEA* m_ea;

  std::string m_rep;  ///< name of report file

  /// external function to call if reproduction was successful
  PrintReproductionFunction PrintReproFunc;

  void* m_data; ///< data to give reproduction function

  // Index into fitness class (which is derived from SpaceLib class Vector).
  // It determines, which fitness value from fitness vector should be used
  // to estimate the winner of a tournament that two individuals fight.
  // Example: m_tfi = 3 => the fourth(3+1; element 1 has index 0) element
  //                       of fitness vector is to be used to determine the
  //                       winner of a tournament
  // If it is not set its default value will be 0.
  int m_tournament_fitness_index;

  // Index into fitness class that should be the time fitness function value.
  // It is used in conjunction with possible negative values of time pressure
  // when determining the winner of a tournament of two individuals.
  // If it is not set its default value will be 0.
  int m_time_fitness_index;


  // general population properties
  int           m_size;
public:
  inline int size(void) {return m_size;};
private:

  double m_fitness_sum;  ///< sum of all population individual fitness values

                ///< number of individuals in the population
  long          m_number_of_reproductions;
                ///< number of reproductions in the current epoch

  CPopulationConvergence *m_convergence;

  CIndividual** m_individual;          ///< array of pointers to the individuals

  // auxiliary class members
  CIndividual*  m_parent1;      ///< first parent individual during reproduction
  CIndividual*  m_parent2;     ///< second parent individual during reproduction
  CIndividual*  m_child1;  ///< auxiliary child individual used for reproduction
  CIndividual*  m_child2;  ///< auxiliary child individual used for reproduction
  bool* m_should_stop;          ///< stop flags; externally set (multithreading)
  bool convergence;

  // __accessors__

  CIndividual* getFirstParent(void) const throw();
  CIndividual* getSecondParent(void) const throw();
  CIndividual* getFirstChild(void) const throw();
  CIndividual* getSecondChild (void) const throw();

public:
  /*! \brief select one individual randomly
   * \return pointer to selected individual
   */  
  CIndividual* getRandomIndividual(void) throw();

private:
  void modifyHypercubeOnGoodEpochReproduction(void);
  void modifyHypercubeOnBadEpochReproduction(void);

public:

  /*!
   *
   */
  void Tournament(CIndividual*  competitor1,
                  CIndividual*  competitor2,
                  CIndividual** winner,
                  CIndividual** loser) throw();

  /*!
   *
   */
  double Reproduction(void) throw();


  /*!
   *
   */
  inline void statusMessage(const std::ostringstream& mess)
  {
    Output(mess.str()); // output through global output object

    // and report file
    if (m_rep.size() != 0)
    {
      std::ofstream report_file(m_rep.c_str(), std::ios::app);
      if (report_file.is_open() == true) report_file << mess.str();
    }
  };


  /*! \brief returns true if flag "m_should_stop" is true, false if otherwise 
   *or not set
   *
   */
  bool abort(void) throw();

  /*!
   * Tbd.
   */
  void checkKeyboard(void);

//=============================================================================
// CPopulation CLASS INTERFACE
//=============================================================================
public:

  inline void abortReproduction(void) { convergence = true;};

  // __constructor__
  /*!
   * \brief Create population and initialize.
   */
  CPopulation(CEA* ea                    // pointer to the calling evolutionary algorithm
             ,const int popsize         // population size
             ,const std::string& rep_file = 0 // name of report file
             ,PrintReproductionFunction PRF = 0 // function to print a reproduction
             ,void* data = 0);           // additionally data for printing of reproduction


  /*!
   * \brief Destructor
   */
  ~CPopulation();

  /*!
   *
   */
  void PrintReproduction( CIndividual* ind) const;


  /*!
   * \brief Sets index of fitness function that is used in tournaments.
   *
   * \param idx Index to be used in tournaments of two individuals.
   */
  void setTournamentFitnessIndex(const int idx) {m_tournament_fitness_index = idx;};


  /*!
   * \brief Returns index of fitness function currently used in tournaments.
   *
   * \return Used index to determine the winner of a tournament of two individuals.
   */
  int  getTournamentFitnessIndex(void) const {return m_tournament_fitness_index;};


  /*!\brief Sets index of time fitness function that is used in tournaments.
   * \param idx Index of time fitness value to be used in tournaments of two
   *individuals. */
  void setTimeFitnessIndex(const int idx) throw();

  /*! \brief index of time fitness function currently used in tournaments.
   * \return Used index of time fitness value. */
  int getTimeFitnessIndex(void) const throw();

  /*! \brief select parent individual
   *  \return pointer to the selected parent individual */
  CIndividual* SelectParent(const CIndividual* const not_this1,
                            const CIndividual* const not_this2,
                            const CIndividual* const not_this3);

  /// reproduction until population convergence has been reached
  double Reproduce(void) throw();

  /// set stop variable pointer
  void setStopBoolPointer(bool*) throw();

  /// get stop variable pointer
  bool* getStopBoolPointer(void) throw();

  /// individual object adressed by given index
  CIndividual* getIndividual(int which) throw();

}; // end class CPopulation

// implementation

//------------------------------------------------------------------------------
inline CIndividual* CPopulation::getIndividual(int which) throw()
{
#ifdef development_version
  assert(which<size());
#endif // development_version
  return m_individual[which];
} // CIndividual* getIndividual(unsigned int which) throw


//------------------------------------------------------------------------------
inline CIndividual* CPopulation::getFirstParent(void) const throw()
{
  return m_parent1;
} // CIndividual* getFirstParent(void) const  throw()

//------------------------------------------------------------------------------
inline CIndividual* CPopulation::getSecondParent(void) const throw()
{
  return m_parent2;
} // CIndividual* CPopulation::getSecondParent(void) const throw()

//------------------------------------------------------------------------------
inline CIndividual* CPopulation::getFirstChild(void) const throw()
{
  return m_child1;
} // CIndividual* CPopulation::getFirstChild(void) const throw()

//------------------------------------------------------------------------------
inline CIndividual* CPopulation::getSecondChild (void) const throw()
{
  return m_child2;
} // inline CIndividual* CPopulation::getSecondChild (void) const throw()

//------------------------------------------------------------------------------
inline bool CPopulation::abort(void) throw()
{
  return ((m_should_stop!=0?(*m_should_stop):false)==true);
} // inline bool abort(void)

//------------------------------------------------------------------------------
inline int CPopulation::getTimeFitnessIndex(void) const throw()
{
  return m_time_fitness_index;
} // int getTimeFitnessIndex(void) const throw()

//------------------------------------------------------------------------------
inline void CPopulation::setTimeFitnessIndex(const int idx) throw()
{
  m_time_fitness_index = idx;
} // void setTimeFitnessIndex(const int idx) throw()

#endif // CPOPULATION_H_INCLUDED
