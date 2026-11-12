/*!
\file CInputEA.cpp
\brief An auxiliarary class that contains the values of the evolutionary algorithm parameter file.

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

// ENCLib
#include <CInputEA.h>

// STL
#include <algorithm>
#include <iomanip>
#include <sstream>

// SpaceLib
#include <Global.h>

using namespace std;
using namespace Global;

// parameter indices to access different arrays
enum PARAM_IDX
{
  TYPE       ///< file type
 ,SS_HCS     ///< search space hypercube size
 ,SS_HCSS    ///< search space hypercube start size
 ,SS_HCSF    ///< search space hypercube shrinking factor
 ,SS_HC_MOD_GOOD_REPRO ///< hypercube modifier after good epoch reproduction
 ,SS_HC_MOD_BAD_REPRO  ///< hypercube modifier after bad epoch reproduction
 ,POP_SIZE   ///< population size
 ,SSS_EPOCHS ///< number of search space scan epochs
 ,CHROM_MP   ///< chromosome mutation probability
 ,GEN_MP     ///< gen mutation probability
 ,SPT        ///< selection pressure on time
 ,FFT        ///< fitness function type
 ,FFT2       ///< fitness function type 2
 ,FFT3       ///< fitness function type 3
 ,HUL        ///< hypercube upper limit
 ,FBC_LIMIT  ///< final boundary condition (FBC) limit
 ,PCL1       ///< population convergence limit 1 (FBC met)
 ,PCL2       ///< population convergence limit 2 (FBC not met)
 ,EAL1       ///< evol. alg. limit 1 (FBC met)
 ,EAL2       ///< evol. alg. limit 2 (FBC not met)
 ,N_MSGS
};


// parameter names
const string PARAM_NAME[N_MSGS]=
{
  "FILE_TYPE"
 ,"SEARCH_SPACE_HYPERCUBE_SIZE"
 ,"HYPERCUBE_START_SIZE"
 ,"HYPERCUBE_SHRINKING_FACTOR"
 ,"GOOD_REPRO_HYPERCUBE_MOD_VALUE"
 ,"BAD_REPRO_HYPERCUBE_MOD_VALUE"
 ,"POPULATION_SIZE"
 ,"SEARCH_SCAN_EPOCHS"
 ,"CHROMOSOME_MUTATION_PROBABILITY"
 ,"GEN_MUTATION_PROBABILITY"
 ,"SELECTION_PRESSURE_ON_TIME"
 ,"FITNESS_FUNCTION_TYPE"
 ,"FITNESS_FUNCTION_TYPE2"
 ,"FITNESS_FUNCTION_TYPE3"
 ,"HYPERCUBE_UPPER_LIMIT"
 ,"FBC_MET_FITNESS"
 ,"IL_POP_CONV_FBC_MET"
 ,"IL_POP_CONV_FBC_NOT_MET"
 ,"IL_EA_CONV_FBC_MET"
 ,"IL_EA_CONV_FBC_NOT_MET"
};


// error messages in case a keyword could not be found
const string ERR_MSG[N_MSGS] =
{
  "CInputEA::CInputEA(): not an InTrance 1.7 evolutionary algorithm parameter file\n"
 ,"CInputEA::CInputEA(): missing parameter (SEARCH_SPACE_HYPERCUBE_SIZE)\n"
 ,""  // search space hypercube start size must not be provided
 ,"CInputEA::CInputEA(): missing parameter (HYPERCUBE_SHRINKING_FACTOR)\n"
 ,"" // hypercube modification after good reproduction mustn't be provided
 ,"" // hypercube modification after bad reproduction mustn't be provided
 ,"CInputEA::CInputEA(): missing parameter (POPULATION_SIZE)\n"
 ,"CInputEA::CInputEA(): missing parameter (SEARCH_SCAN_EPOCHS)\n"
 ,"CInputEA::CInputEA(): missing parameter (CHROMOSOME_MUTATION_PROBABILITY)\n"
 ,""  // gen mutation probability must not be provided
 ,"CInputEA::CInputEA(): missing parameter (SELECTION_PRESSURE_ON_TIME)\n"
 ,""  // fitness function type must not be provided
 ,""  // fitness function type 2 must not be provided
 ,""  // fitness function type 3 must not be provided
 ,""  // hypercube upper limit must not be provided
 ,""  // final boundary condition met limit must not be provided
 ,""  // epoch convergence limit 1 (FBC met) must not be provided
 ,""  // epoch convergence limit 2 (FBC not met) must not be provided
 ,""  // ea convergence limit 1 (FBC met) must not be provided
 ,""  // ea convergence limit 2 (FBC not met) must not be provided
};


//------------------------------------------------------------------------------
void CInputEA::init(void)
{
  // convergence limits
  if (convergence_limits.size() != N_LIMITS) 
    convergence_limits.resize(N_LIMITS);
  for (unsigned int i=0;i<convergence_limits.size();++i)
    convergence_limits[i]=DEFAULT_CONV_LIMITS[i];

  m_vectornorm                       = DEFAULT_VECTOR_NORM;
  min_fitness_value                   = DEFAULT_MIN_FITNESS_VALUE;
  fitness_function_type               = DEFAULT_FITNESS_FUNCTION_TYPE;
  fft2                                = DEFAULT_FITNESS_FUNCTION_TYPE2;
  fft3                                = DEFAULT_FITNESS_FUNCTION_TYPE3;
  hypercube_shrinking_factor          = DEFAULT_HYPERCUBE_SHRINKING_FACTOR;
  hypercube_upper_limit               = DEFAULT_HYPERCUBE_UPPER_LIMIT;
  initial_hypercube_size              = DEFAULT_HYPERCUBE_INIT_SIZE;
  hypercube_start_size                = DEFAULT_HYPERCUBE_START_SIZE;

  setHypercubeFactorOnGoodEpochRepro(DEFAULT_GOOD_REPRODUCTION_HC_FACTOR);
  setHypercubeFactorOnBadEpochRepro(DEFAULT_BAD_REPRODUCTION_HC_FACTOR);

  m_identifier                        = EAP_FILE_IDENTIFIER;
  maximum_number_of_sss_epochs        = DEFAULT_MAX_HYPERCUBE_SSS_EPOCHS;
  max_n_bad_epochs                    = DEFAULT_MAX_N_BAD_EPOCHS;
  mutation_probability_per_chromosome = DEFAULT_CHROM_MUTATION_PROBABILITY;
  mutation_probability_per_gen        = DEFAULT_GEN_MUTATION_PROBABILITY;
  population_size                     = DEFAULT_POP_SIZE;
  time_selection_pressure             = DEFAULT_TIME_SELECTION_PRESSURE;
}



/*!
 * \brief Parameter list constructor (use with configuration files).
 */
//------------------------------------------------------------------------------
CInputEA::CInputEA(ParamList& p, const string& trailing_chars)
  throw(std::invalid_argument)
{
  init();    // initialize member fields

  //check for file type and version
  // Usually "InputStringParameter" should be put into a try-catch
  // block because it may throw a std::invalid_argument type
  // exception if there's no parameter named "FILE_TYPE" could be
  // found. The right way to react to this, is to end this constructor and
  // to let the exception be handled by the caller.
  string version;
  InputStringParameter(PARAM_NAME[TYPE], version, p
    , ERR_MSG[TYPE], trailing_chars);
  removeChars(version, " ");
  upper(version);
  if (version!=m_identifier) throw invalid_argument(ERR_MSG[TYPE]);

  // template for including of new parameters
  //// your comment
  //InputParameter(PARAM_NAME[new_parameter_index], your_member field, p
  //  , ERR_MSG[new_parameter_index], trailing_chars, false_or_true);

  // initial search space hypercube size
  InputParameter(PARAM_NAME[SS_HCS], initial_hypercube_size, p, ERR_MSG[SS_HCS], trailing_chars);
  // initialize start size with initial size
  setHypercubeStartSize(getHypercubeInitialSize());

  // and overwrite with search space hypercube start size, if given
  InputParameter(PARAM_NAME[SS_HCSS], hypercube_start_size, p, ERR_MSG[SS_HCSS], trailing_chars, false);

  // hypercube shrinking factor
  InputParameter(PARAM_NAME[SS_HCSF], hypercube_shrinking_factor, p, ERR_MSG[SS_HCSF]
  , trailing_chars);

  // hypercubesize modification factor after good epoch reproduction
  InputParameter(PARAM_NAME[SS_HC_MOD_GOOD_REPRO]
   , hypercube_resize_factor_on_good_repro, p, ERR_MSG[SS_HC_MOD_GOOD_REPRO]
   , trailing_chars, false);

  // hypercubesize modification factor after good epoch reproduction
  InputParameter(PARAM_NAME[SS_HC_MOD_BAD_REPRO]
   , hypercube_resize_factor_on_bad_repro, p, ERR_MSG[SS_HC_MOD_BAD_REPRO]
   , trailing_chars, false);

  // population size
  InputParameter(PARAM_NAME[POP_SIZE], population_size, p, ERR_MSG[POP_SIZE]
  , trailing_chars);

  // number of search space scan epochs
  InputParameter(PARAM_NAME[SSS_EPOCHS], maximum_number_of_sss_epochs, p, ERR_MSG[SSS_EPOCHS]
  , trailing_chars);

  // chromosome mutation probability
  InputParameter(PARAM_NAME[CHROM_MP], mutation_probability_per_chromosome, p, ERR_MSG[CHROM_MP]
  , trailing_chars);

  // gen  mutation probability
  InputParameter(PARAM_NAME[GEN_MP], mutation_probability_per_gen, p, ERR_MSG[GEN_MP]
  , trailing_chars, false);

  // selection pressure on transfer time
  InputParameter(PARAM_NAME[SPT], time_selection_pressure, p, ERR_MSG[SPT]
  , trailing_chars);

  // fitness function type
  try
  {
    InputStringParameter(PARAM_NAME[FFT], fitness_function_type, p, ERR_MSG[FFT]
    , trailing_chars);
    upper(fitness_function_type);
  }
  catch (std::invalid_argument) {} // do nothing, because it mustn't be provided


  // secondary fitness function type
  InputParameter(PARAM_NAME[FFT2], fft2, p, ERR_MSG[FFT2]
  , trailing_chars, false);

  // fitness function type 3
  InputParameter(PARAM_NAME[FFT3], fft3, p, ERR_MSG[FFT3]
  , trailing_chars, false);

  // hypercube upper limit
  InputParameter(PARAM_NAME[HUL], hypercube_upper_limit, p, ERR_MSG[HUL]
  , trailing_chars, false);

  // final boundary conditions met
  InputParameter(PARAM_NAME[FBC_LIMIT], convergence_limits[FBC_MET]
  , p, ERR_MSG[FBC_LIMIT], trailing_chars, false);

  // epoch convergence limit 1 (FBC met)
  InputParameter(PARAM_NAME[PCL1], convergence_limits[POP_CL_FBC_MET]
  , p, ERR_MSG[PCL1], trailing_chars, false);

  // epoch convergence limit 2 (FBC not met)
  InputParameter(PARAM_NAME[PCL2], convergence_limits[POP_CL_FBC_NMET]
  , p, ERR_MSG[PCL2], trailing_chars, false);

  // evolutionary algorithm convergence limit 1 (FBC met)
  InputParameter(PARAM_NAME[EAL1], convergence_limits[EA_CL_FBC_MET]
  , p, ERR_MSG[EAL1], trailing_chars,false);

  // evolutionary algorithm convergence limit 2 (FBC not met)
  InputParameter(PARAM_NAME[EAL2], convergence_limits[EA_CL_FBC_NMET]
  , p, ERR_MSG[EAL2], trailing_chars, false);
} // CInputEA::CInputEA(ParamList&,...)


/*!
 * \brief direct constructor.
 */
//------------------------------------------------------------------------------
CInputEA::CInputEA( const double              p_init_hc_sz
                   ,const double              p_hc_shrink_fac
                   ,const int                 p_pop_sz
                   ,const long                p_max_n_sss_epoches
                   ,const double              p_mut_prob_per_chrom
	           ,const double              p_mut_prob_per_gen
                   ,const double              p_t_selec_press
                   )
{
  init(); // pre - initialize member fields
  setHypercubeSize(p_init_hc_sz);
  setHypercubeStartSize(getHypercubeSize());
  setHypercubeShrinkingFactor(p_hc_shrink_fac);
  setPopulationSize(p_pop_sz);
  setMaxSSSEpochs(p_max_n_sss_epoches);
  setMutProbPerChrom(p_mut_prob_per_chrom);
  setMutProbPerGen(p_mut_prob_per_gen);
  setTimeSelectionPressure(p_t_selec_press);
}


//------------------------------------------------------------------------------
void CInputEA::printSettings(void) const
{
  const unsigned int w(37);

  // initial hypercube size
  ostringstream ss;
  ss << left << setw(w) << "EA: initial hypercube size"<<": " 
     << initial_hypercube_size << endl;Output(ss.str());

  // hypercube start size (for warmstart purposes)
  ss.str("");
  ss << left << setw(w) << "EA: hypercube start size"<<": " 
     << hypercube_start_size <<endl;Output(ss.str());

  // minimum fitness value
  ss.str("");
  ss << left << setw(w) << "EA: minimum fitness value"<<": " 
     << min_fitness_value <<endl;Output(ss.str());

  // hypercube shrinking factor
  ss.str("");
  ss << left << setw(w) << "EA: hypercube shrinking factor"<<": " 
     << hypercube_shrinking_factor <<endl;Output(ss.str());

  // population size
  ss.str("");
  ss << left << setw(w) << "EA: population size" << ": " << population_size <<endl;Output(ss.str());

  // # of search space epochs
  ss.str("");
  ss << left << setw(w) << "EA: search space scan epochs" << ": " << maximum_number_of_sss_epochs 
     << endl;Output(ss.str());

  // mutation propability
  ss.str("");
  ss << left << setw(w) << "EA: chromosome mutation probability" << ": " 
     << mutation_probability_per_chromosome << endl;Output(ss.str());

  // gen mutation probability
  ss.str("");
  ss << left << setw(w) << "EA: gen mutation probability" << ": "
     << mutation_probability_per_gen << endl;Output(ss.str());

  // time selection pressure
  ss.str("");
  ss << left << setw(w) << "EA: selection pressure on time" << ": " 
     << time_selection_pressure << endl;Output(ss.str());

  // primary fitness function type
  ss.str("");
  ss << left << setw(w) << "EA: fitness function type" << ": " << fitness_function_type;
  if (fitness_function_type.compare(DEFAULT_FITNESS_FUNCTION_TYPE)==0 ) 
    ss << " (default)";
  ss << endl;
  Output(ss.str());

  // secondary fitness function type
  ss.str("");
  ss << left << setw(w) << "EA: fitness function type 2" << ": " << fft2;
  if (fft2==DEFAULT_FITNESS_FUNCTION_TYPE2) ss<< " (default)";
  ss << endl;Output(ss.str());

  // third fitness function type
  ss.str("");
  ss << left << setw(w) << "EA: fitness function type 3" << ": " << fft3;
  if (fft3==DEFAULT_FITNESS_FUNCTION_TYPE3) ss<< " (default)";
  ss << endl;Output(ss.str());

  // hypercube upper limit
  ss.str("");
  ss << left << setw(w) << "EA: hypercube upper limit" << ": " << hypercube_upper_limit;
  if (hypercube_upper_limit==DEFAULT_HYPERCUBE_UPPER_LIMIT) ss<< " (default)";
  ss << endl;Output(ss.str());

  // final boundary conditions met
  ss.str("");
  ss << left << setw(w) << "EA: final bound. cond. (FBC) limit" << ": " 
     << convergence_limits[FBC_MET];
  if (convergence_limits[FBC_MET]==0.0) ss<< " (default)";
  ss << endl;Output(ss.str());

  // epoch convergence limit 1 (FBC met)
  ss.str("");
  ss << left << setw(w) << "EA: pop. conv. limit 1 (FBC met)" << ": " 
     << convergence_limits[POP_CL_FBC_MET];
  if (convergence_limits[POP_CL_FBC_MET]==0.0) ss << " (default)";
  ss<<endl;Output(ss.str());

  // epoch convergence limit 2 (FBC not met)
  ss.str("");
  ss << left << setw(w) << "EA: pop. conv. limit 2 (FBC not met)" << ": " 
     << convergence_limits[POP_CL_FBC_NMET];
  if (convergence_limits[POP_CL_FBC_NMET]==0.0) ss << " (default)";
  ss<<endl;Output(ss.str());

  // evolutionary algorithm convergence limit 1 (FBC met)
  ss.str("");
  ss << left << setw(w) << "EA: ea conv. limit 1 (FBC met)" << ": " 
     << convergence_limits[EA_CL_FBC_MET];
  if (convergence_limits[EA_CL_FBC_MET]==0.0) ss << " (default)";
  ss<<endl;Output(ss.str());

  // evolutionary algorithm convergence limit 2 (FBC not met)
  ss.str("");
  ss << left << setw(w) << "EA: ea conv. limit 2(FBC not met)" << ": " 
     << convergence_limits[EA_CL_FBC_NMET];
  if (convergence_limits[EA_CL_FBC_NMET]==0.0) ss<< " (default)";
  ss<<endl;Output(ss.str());

} // CInputEA::printSettings(void) const

