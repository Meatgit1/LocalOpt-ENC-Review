/*!
\file CInputEA.h
\brief An auxiliarary class that contains the values of the
       evolutionary algorithm parameter file.

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

  05.10.2006 (AO)

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
#ifndef CINPUTEA_H_INCLUDED
#define CINPUTEA_H_INCLUDED

// C++ STL
#include<string>
#include<valarray>

// ENCLib
#include <configENC.h>
#include <CEATypes.h>

// SpaceLib
#include <Global.h>


/*!
 * \brief Enumeration for handling convergence limits.
 */
enum CONV_LIMIT
{
   FBC_MET           ///< final boundary condition (FBC) met
  ,POP_CL_FBC_MET    ///< population convergence limit 1 (FBC met)
  ,POP_CL_FBC_NMET   ///< population convergence limit 2 (FBC not met)
  ,EA_CL_FBC_MET     ///< evol. algo. convergence limit 1 (FBC met)
  ,EA_CL_FBC_NMET    ///< evol. algo. convergence limit 2 (FBC not met)
  ,N_LIMITS          ///< number of limit values (needed data structure size)
};

typedef std::valarray<double> conv_limits;
///< definition of convergence limits data type


///< default convergence limit values
const double DEFAULT_CONV_LIMITS[N_LIMITS] =
{
  5.0,5.0,5.0,5.0,5.0
};

const VectorNorm   DEFAULT_VECTOR_NORM = MAXIMUM_SUM_NORM;
const double       DEFAULT_MIN_FITNESS_VALUE = -1.0e10;
const std::string  DEFAULT_FITNESS_FUNCTION_TYPE = "NO_FFT";
const int          DEFAULT_FITNESS_FUNCTION_TYPE2 = 0;
const int          DEFAULT_FITNESS_FUNCTION_TYPE3 = 0;
const double       DEFAULT_HYPERCUBE_SHRINKING_FACTOR = 0.5;
const double       DEFAULT_HYPERCUBE_UPPER_LIMIT = 0.5;
const double       DEFAULT_HYPERCUBE_INIT_SIZE = 1.0;
const double       DEFAULT_HYPERCUBE_START_SIZE = DEFAULT_HYPERCUBE_INIT_SIZE;
const std::string  EAP_FILE_IDENTIFIER = "INTRANCE1.7EVOLUTIONARYALGORITHMPARAMETERFILE";
const unsigned int DEFAULT_MAX_HYPERCUBE_SSS_EPOCHS = 5;
const int          DEFAULT_MAX_N_BAD_EPOCHS = 2;
const double       DEFAULT_CHROM_MUTATION_PROBABILITY = 0.90;
const double       DEFAULT_GEN_MUTATION_PROBABILITY = 0.0;
const unsigned int DEFAULT_POP_SIZE = 25;
const double       DEFAULT_TIME_SELECTION_PRESSURE = 0.0;
const double       DEFAULT_BAD_REPRODUCTION_HC_FACTOR = 1.0;
const double       DEFAULT_GOOD_REPRODUCTION_HC_FACTOR = 1.0;


/*! \brief evolutionary algorithm input object class
 *
 * This is a auxiliary class containing the values of the evolutionary algorithm 
 * input file for initializing CEA class objects.
 */
class CInputEA
{
//=============================================================================
// CInputEA CLASS INTERFACE
//=============================================================================
public:

  std::string m_identifier;
  ///< identifier string for validity checks of input files

  VectorNorm      m_vectornorm;

  double min_fitness_value;      ///< minimum fitness value
  /// set minimum fitness value
  void setMinFitnessValue(const double val) { min_fitness_value = val; }
  /// get minimum fitness value
  double getMinFitnessValue(void) const { return min_fitness_value; }
  

  double initial_hypercube_size;
  ///< size of initial search space hypercube

  double hypercube_start_size;
  ///< start value of search space hypercube size

  double hypercube_shrinking_factor;
  ///< hypercube shrinking factor

  double hypercube_upper_limit;
  ///< maximum size of hypercube after optimization

  double hypercube_resize_factor_on_good_repro;
  ///< hypercube resize factor after good reproduction
  ///< (new best epoch solution)

  double hypercube_resize_factor_on_bad_repro;
  ///< hypercube resize factor after a bad epoch reproduction
  ///< (no new best epoch reproduction)

  int chromosome_length; ///< # of loci in chromsome
  ///< chromosome length, as defined by the neurocontroller object
  ///< (= number of parameters to be optimized)

  /// set length of chromosome = number of chromosome's loci
  void setChromosomeLength(const int value)
  {
    chromosome_length = value;
  }
  /// chromosome loci count
  int getChromosomeLength(void) const
  {
    return chromosome_length;
  }


  int population_size;
  ///< population size

  int maximum_number_of_sss_epochs;
  ///< maximum number of search space scan epochs

  int max_n_bad_epochs;
  ///< max. # of bad epochs before hypercube size will be changed

  ///  mutation_probability per chromosome
  double mutation_probability_per_chromosome;

  /// probability that chromosome gen is changed
  double mutation_probability_per_gen;

  double time_selection_pressure;
  ///< selection pressure on transfer time

  std::string fitness_function_type;
  ///< primary fitness function type

  int fft2;                      ///< fitness function type 2 (must not be used)
  int fft3;                      ///< fitness function type 2 (must not be used)

  conv_limits convergence_limits;
  ///< array of convergence limit values

private:
  void init(void);
public:


  /*!
   * \brief Standard constructor.
   * \todo AO set standard values for convergence limits
   */
  //----------------------------------------------------------------------------
  CInputEA(void)
  {
    init();
  }


  /*!
   * \brief Parameter list constructor (use with configuration files).
   */
  //----------------------------------------------------------------------------
  CInputEA(ParamList& parameters, const std::string& trailing_chars = " \t") 
    throw(std::invalid_argument);



  /*!
   * \brief Direct constructor.
   */
  //----------------------------------------------------------------------------
  CInputEA( const double                 p_init_hc_sz
           ,const double                 p_hc_shrink_fac
           ,const int                    p_pop_sz
           ,const long                   p_max_n_sss_epoches
           ,const double                 p_mut_prob_per_chrom
	   ,const double                 p_mut_prob_per_gen
           ,const double                 p_t_selec_press
          );

  /*!
   * \brief Helper function to print current setting to output object ( file, screen)
   */
  //----------------------------------------------------------------------------
  void printSettings(void) const;

  /// initial size of search space hypercube
  double getHypercubeInitialSize(void) const
  {
    return initial_hypercube_size;
  }


  /// set new initial hypercube dimension
  /*! \param value new initial hypercube dimension*/
  void setHypercubeInitialSize(const double value)
  {
    initial_hypercube_size = value;
  }

  /*!
   * \brief hypercube size at the beginning of optimisation
   */
  //----------------------------------------------------------------------------
  double getHypercubeStartSize(void) const
  {
    return hypercube_start_size;
  }

  /// number of max. bad epochs after which a hypercube size change is made
  int getMaxBadEpochsCount(void) const { return max_n_bad_epochs; };

  /*! \brief sets max. count value of bad epochs
   * \param count the maximum count value to set*/
  void setMaxBadEpochsCount(const int count) { max_n_bad_epochs = count; };


  /*! \brief file identifier string
   *
   * This method returns the string, by which a configuration file is recognized
   * as a evolutionary algorithm input file.\n
   * So, to make a text file to such a configuration file it must contain the
   * following line:\n
   * FILE_TYPE=identifer string returned by this method
   *
   * \return evolutionary algorithm configuration file identifier string
   */
  //----------------------------------------------------------------------------
  std::string& getIdentifier(void)
  {
    return m_identifier;
  }


  /*! \brief Returns fitness function type string.
   *
   */
  //----------------------------------------------------------------------------
  const std::string& getFitnessFunctionType(void) const
  {
    return fitness_function_type;
  }


  /*! \brief Sets fitness function type string.
   */
  //----------------------------------------------------------------------------
  std::string& setFitnessFunctionType(const std::string& s)
  {
    fitness_function_type=s;
  }

  //----------------------------------------------------------------------------
  /// \brief Return secondary fitness function type
  int getFitnessFunctionType2(void) const { return fft2;}

  //----------------------------------------------------------------------------
  /// \brief Return fitness function type 3
  int getFitnessFunctionType3(void) const { return fft3;}

  //----------------------------------------------------------------------------
  /// \brief Set secondary fitness function type
  int setFitnessFunctionType2(const int p_fft2) { fft2 = p_fft2;}

  //----------------------------------------------------------------------------
  /// \brief Set secondary fitness function type
  int setFitnessFunctionType3(const int p_fft3) { fft3 = p_fft3;}

  /*!\brief Sets limits of improvement for convergence determination.
   *
   */
  //----------------------------------------------------------------------------
  void setConvergenceLimits(const conv_limits& limits)
  {
    convergence_limits = limits;
  }

  /*! \brief Set single convergence limit.
  */
  //----------------------------------------------------------------------------
  void setConvergenceLimit(const CONV_LIMIT i, const double val)
  {
    if (i<N_LIMITS) convergence_limits[i]=val;
  }

  /*!\brief Get limits of improvement for convergence determination.
   * \return desired limits
   */
  //----------------------------------------------------------------------------
  const conv_limits& getConvergenceLimits(void) const
  {
    return convergence_limits;
  }

  /*! \brief Returns single convergence limit.
  */
  //----------------------------------------------------------------------------
  double getConvergenceLimit(const CONV_LIMIT i) const
  {
    return convergence_limits[i];
  }

  /*! \brief Returns whether default convergence limit set or not.
   */
  bool isDefaultConvergenceLimit(const CONV_LIMIT i)
  {
    return (convergence_limits[i] == DEFAULT_CONV_LIMITS[i]);
  }

   //accessor/modifier functions

   /*! \brief return hypercube reproduction modification factor 
   * (after bad reproduction)
   *
   * Current hypercube size is multiplied with this factor after every 
   * bad reproduction.
   *
   * \sa getHypercubeFactorOnBadRepro
   */
   double getHypercubeFactorOnBadEpochRepro(void) const
   {
     return hypercube_resize_factor_on_bad_repro;
   };

   /*! \brief return hypercube reproduction modification factor 
   * (after good reproduction)
   *
   * Current hypercube size is multiplied with this factor after every 
   * good reproduction.
   *
   * \sa getHypercubeFactorOnGoodRepro
   */
   double getHypercubeFactorOnGoodEpochRepro(void) const
   {
     return hypercube_resize_factor_on_good_repro;
   };


   /*! \brief set hypercube reproduction modification factor (bad reproduction)
    *
    * current hypercube size is multiplied with this factor after bad reproduction*/
   void setHypercubeFactorOnBadEpochRepro(const double factor)
   {
     hypercube_resize_factor_on_bad_repro = factor;
   };

   /*! \brief set hypercube reproduction modification factor (good reproduction)
    *
    * current hypercube size is multiplied with this factor after good reproduction*/
   void setHypercubeFactorOnGoodEpochRepro(const double factor)
   {
     hypercube_resize_factor_on_good_repro = factor;
   };


   /*!
    * \brief Set initial size of search space hypercube.
    */
   //----------------------------------------------------------------------------
   void setHypercubeSize(const double p_init_hc_sz) 
   {
     initial_hypercube_size = p_init_hc_sz;
   }


   /*!
    * \brief Set start value of hypercube.
    */
   //----------------------------------------------------------------------------
   void setHypercubeStartSize(const double value)
   {
     hypercube_start_size = value;
   }


  /*!
   * \brief Return size of search space hypercube.
   */
  //----------------------------------------------------------------------------
   double getHypercubeSize(void) const 
   {
     return initial_hypercube_size;
   }


  /*!
   * \brief Set hypercube shrinking factor.
   */
  //----------------------------------------------------------------------------
  void   setHypercubeShrinkingFactor(const double p_hc_shrink_fac) {
    hypercube_shrinking_factor = p_hc_shrink_fac;
  }


  /*!
   * \brief Return hypercube shrinking factor.
   */
  //----------------------------------------------------------------------------
  double getHypercubeShrinkingFactor(void) const {
    return hypercube_shrinking_factor;
  }


  /*!
   * \brief Set hypercube upper limit.
   */
  //----------------------------------------------------------------------------
  void   setHypercubeUpperLimit(const double p_hc_upper_limit) {
    hypercube_upper_limit = p_hc_upper_limit;
  }


  /*!
   * \brief Return hypercube upper limit.
   */
  //----------------------------------------------------------------------------
  double getHypercubeUpperLimit(void) const {
    return hypercube_upper_limit;
  }
  

  /*!
   * \brief Set size of population of individuals.
   */
  //----------------------------------------------------------------------------
  void setPopulationSize(const int p_pop_sz) {
    population_size = p_pop_sz;
  }


  /*!
   * \brief Return size of population (number of individuals in population).
   */
  //----------------------------------------------------------------------------
  int getPopulationSize(void) const { return population_size;};


  /*!
   * \brief Set maximum number of search space scan epoches.
   */
  //----------------------------------------------------------------------------
  void setMaxSSSEpochs(const long p_max_n_sss_epochs) {
    maximum_number_of_sss_epochs = p_max_n_sss_epochs;
  }


  /*!
   * \brief Set maximum number of search space scan epoches.
   */
  //----------------------------------------------------------------------------
  long getMaxSSSEpochs(void) const {
    return maximum_number_of_sss_epochs;
  }


  /*!
   * \brief Set mutation probability per chromosome.
   */
  //----------------------------------------------------------------------------
  void setMutProbPerChrom(const double p_mut_prob_per_chrom) {
    mutation_probability_per_chromosome = p_mut_prob_per_chrom;
  }

  /*!
   * \brief Set mutation probability per gen.
   */
  //----------------------------------------------------------------------------
  void setMutProbPerGen(const double p_mut_prob_per_gen) {
    mutation_probability_per_gen = p_mut_prob_per_gen;
  }


  /*!
   * \brief Return mutation probability per chromosome.
   */
  //----------------------------------------------------------------------------
  double getMutProbPerChrom(void) const {
    return mutation_probability_per_chromosome;
  }

  /// probability of mutation for every allele of chromosome
  //----------------------------------------------------------------------------
  double getMutProbPerGen(void) const
  {
    return mutation_probability_per_gen;
  }


  /*!
   * \brief Set time selection pressure.
   */
  //----------------------------------------------------------------------------
  void setTimeSelectionPressure(const double p_t_selec_press) {
    time_selection_pressure = p_t_selec_press;
  }


  /*!
   * \brief Return time selection pressure.
   */
  //----------------------------------------------------------------------------
  double getTimeSelectionPressure(void) const {
    return time_selection_pressure;
  }

  /// used vector norm
  //----------------------------------------------------------------------------
  VectorNorm getVectorNorm(void) const 
  {
    return m_vectornorm;
  }

  /// set new vector norm
  //----------------------------------------------------------------------------
  void setVectorNorm(const VectorNorm norm)
  {
    m_vectornorm = norm;
  }

  /// destructor
  //----------------------------------------------------------------------------
  ~CInputEA() {};

  /// resets input object
  //----------------------------------------------------------------------------
  inline void reset(void) {init();};

}; // end class CInputEA


#endif // CINPUTEA_H_INCLUDED
