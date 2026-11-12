/*!
\file CEA.cpp
\brief A class that provides the functionality for an evolutionary
       algorithm.

\version 1.7

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


\date 28.10.2005 (AO)

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

// STANDARD C
#include <assert.h>

// STL
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

// SpaceLib
#include <Global.h>
#include <MLGlobal.h>

// ENCLib
#include <CNC.h>
#include <CInputEA.h>
#include <CEAFitness.h>
#include <CEA.h>
#include <CIndividual.h>
#include <CPopulation.h>
#include <CEAOperator.h>
#include <CEAEvalFunc.h>
#include <ENCObserver.h>

// used namespaces
using namespace std;
using namespace Global;
using namespace MathLib;

//------------------------------------------------------------------------------
void CEA::uniformSpreadSolutionAlleles(const unsigned int locus,
                                       const double l_b, const double u_b)
{
  CPopulation* pop = getCurrentPopulationObject();
  if (pop != 0)
  {
    const unsigned int sz = pop->size();
    if (sz>2)
    {
      unsigned int i;
      double ub = u_b;//m_data.getHypercubeInitialSize();
      double lb = l_b;//-ub;
      double span = (ub-lb)/(sz-1);
      double value = (lb-span);
      for (i=0;i<sz;++i)
      {
        pop->getRandomIndividual()->setAllele(locus,value+=span);
        //pop->getIndividual(i)->setAllele(locus,value+=span);
        //value+=span;
      }
    } // if (sz>2)
  } // if (getCurrentPopulationObject() != 0)
} // void uniformSpreadSolutionAlleles(const unsigned int loci)

//------------------------------------------------------------------------------
void CEA::printHypercubeShrinkFactorSetting(const int column_width) const throw()
{
  std::ostringstream ss("");
  ss << left << setw(column_width) << "EA: hypercube shrinking factor" << ": " 
    << getHypercubeShrinkFactor() << endl;Output(ss.str());
} // void printHypercubeShrinkFactorSetting(const int column_width) const


//------------------------------------------------------------------------------
void CEA::printMutationRateSetting(const int column_width) const throw()
{
  std::ostringstream ss("");
  ss << left << setw(column_width)<<"EA: chromosome mutation probability"<<": " 
    << getChromMutationRate() << endl;Output(ss.str());
} // void printMutationRateSetting(const int column_width) const


//------------------------------------------------------------------------------
void CEA::printMaxBadEpochsSetting(const int column_width) const throw()
{
  std::ostringstream ss("");
  ss << left << setw(column_width) << "EA: max. # of bad epochs" << ": " 
    << getMaxBadEpochsCount() << endl;Output(ss.str());
} // void CEA::printMaxBadEpochsSetting(const int column_width) const

//------------------------------------------------------------------------------
void CEA::printBadEpochsStatus(
 const int column_width, const bool with_max) const throw()
{
  std::ostringstream ss(""),helper("");
  helper << "EA: # of bad epochs";
  if (with_max==true) helper << " (max.)";
  ss <<left<<setw(column_width)<<helper.str()<<": "<<getBadEpochsCount();
  if (with_max==true) ss << " (" << getMaxBadEpochsCount() << ")";
  ss << endl;Output(ss.str());
} // void printBadEpochsStatus(const int , const bool ) const


//------------------------------------------------------------------------------
void PrintEaSettings(const CEA* ea)
{
  const unsigned int w(37);

  // initial hypercube size
  ostringstream ss;
  ss << left << setw(w) << "EA: initial hypercube size"<<": " 
    << ea->getHypercubeInitialSize() << endl;Output(ss.str());

  ea->printHypercubeShrinkFactorSetting(w);

  // population size
  ss.str("");
  ss << left << setw(w) << "EA: population size" << ": " 
    << ea->getPopulationSize() <<endl;Output(ss.str());

  // # of search space epochs
  ss.str("");
  ss << left << setw(w) << "EA: search space scan epochs" << ": " 
    << ea->getMaxSSSEpochsCount() << endl;Output(ss.str());

  ea->printBadEpochsStatus(w, true);
  ea->printMutationRateSetting(w);

  // gen mutation probability
  ss.str("");
  ss << left << setw(w) << "EA: gen mutation probability" << ":"
    << ea->getGenMutationRate() << endl;Output(ss.str());

  // time selection pressure
  ss.str("");
  ss << left << setw(w) << "EA: selection pressure on time" << ": " 
    << ea->getTimeSelectionPressure() << endl;Output(ss.str());

  // hypercube upper limit
  ss.str("");
  ss << left << setw(w) << "EA: hypercube upper limit" << ": " 
    << ea->getHypercubeHighLimit();
  ss << endl;Output(ss.str());

  // final boundary conditions met
  ss.str("");
  ss << left << setw(w) << "EA: final bound. cond. (FBC) limit" << ": " 
    << ea->m_final_bound_cond_met_limit;
  ss << endl;Output(ss.str());

  // epoch convergence limit 1 (FBC met)
  ss.str("");
  ss << left << setw(w) << "EA: pop. conv. limit 1 (FBC met)" << ": " 
    << ea->m_pop_conv_fbc_met;
  ss<<endl;Output(ss.str());

  // epoch convergence limit 2 (FBC not met)
  ss.str("");
  ss << left << setw(w) << "EA: pop. conv. limit 2 (FBC not met)" << ": " 
    << ea->m_pop_conv_fbc_not_met;
  ss<<endl;Output(ss.str());

  // evolutionary algorithm convergence limit 1 (FBC met)
  ss.str("");
  ss << left << setw(w) << "EA: ea conv. limit 1 (FBC met)" << ": " 
    << ea->m_epoch_conv_limit_final_bound_cond_met;
  ss<<endl;Output(ss.str());

  // evolutionary algorithm convergence limit 2 (FBC not met)
  ss.str("");
  ss << left << setw(w) << "EA: ea conv. limit 2(FBC not met)" << ": " 
    << ea->m_epoch_conv_limit_final_bound_cond_not_met;
  ss<<endl;Output(ss.str());
} // void PrintEaSettings(const CEA* ea)


//------------------------------------------------------------------------------
void CEA::updateStatusObject(void) throw()
{
  if (m_status!=0)
  {
    m_status->population_size = m_population->size();
    m_status->mutation_probability = float(getChromMutationRate());
    m_status->epoch = m_epoch;
    m_status->bad_epochs = getBadEpochsCount();
    m_status->sss_epoch = m_sss_epoch;
    m_status->n_sss_epochs = getMaxSSSEpochsCount();
    m_status->hc_size = float(m_delta);
    m_status->hc_upper_limit = float(getHypercubeHighLimit());
    m_status->hc_sf = float(getHypercubeShrinkFactor());
    m_status->hc_init_size = float(getHypercubeInitialSize());
    m_status->fbc_met_limit = float(m_final_bound_cond_met_limit);
    m_status->pop_conv_limit1 = float(m_pop_conv_fbc_not_met);
    m_status->pop_conv_limit2 = float(m_pop_conv_fbc_met);
    m_status->ea_conv_limit1= float(m_epoch_conv_limit_final_bound_cond_not_met);
    m_status->ea_conv_limit2 = float(m_epoch_conv_limit_final_bound_cond_met);
    m_status->best_fitness_epoch = float(m_best_individual_epoch_fitness);
    m_status->best_fitness_ever = float(m_best_individual_ever_fitness);
  } // if (m_status!=0)
} // void CEA::updateStatusObject(void)

//------------------------------------------------------------------------------
void CEA::abortCurrentEpoch(void) throw()
{
  if (m_population!=0) m_population->abortReproduction();
} // void CEA::abortCurrentEpoch(void)

// constructor
//------------------------------------------------------------------------------
CEA::CEA(
  const CInputEA& input_eap
 ,CNC* nc
 ,CEAEvalFunc* evalfunc
 ,const string& rep_file
 ,const PrintReproductionFunction PRF
 ,const NewHcSizeHandler new_hc_size_handler
 ,void* p_data
):
   m_data(input_eap)
  ,m_rep(rep_file)
  ,m_evalfunc(evalfunc)
  ,m_population(0)
  ,m_pop_conv_fbc_not_met(input_eap.convergence_limits[POP_CL_FBC_NMET])
  ,m_pop_conv_fbc_met(input_eap.convergence_limits[POP_CL_FBC_MET])
  ,m_final_bound_cond_met_limit(input_eap.convergence_limits[FBC_MET])
  ,m_epoch_conv_limit_final_bound_cond_met(input_eap.convergence_limits[EA_CL_FBC_MET])
  ,m_epoch_conv_limit_final_bound_cond_not_met(input_eap.convergence_limits[EA_CL_FBC_NMET])
  ,m_operator(NULL)
  ,m_reproduction(0)
  ,m_sss_epoch(0)
  ,m_epoch(0)
  ,m_delta(input_eap.getHypercubeStartSize())
  ,m_delta_best(m_delta)
  ,m_should_stop(0)
  ,m_observer(0)
  ,m_preparation_smoke_enabled(false)
  ,m_preparation_stop(false)
  ,m_preparation_search_calls(0)
  ,m_preparation_pairs(0)
  ,m_ea_converged(false)
  ,m_partial_solution()
  ,m_partial_solution_fitness(input_eap.min_fitness_value)
  ,m_best_individual_epoch()
  ,m_best_individual_epoch_fitness(input_eap.min_fitness_value)
  ,m_best_individual_ever()
  ,m_best_individual_ever_fitness(input_eap.min_fitness_value)
  ,m_avg_fitness_convergence(false)
  ,m_trigger_avg_fitness_ratio(0.09)
  ,m_dynamic_mutation_rate(false)
  ,m_nc()
  ,m_status(0)
  ,m_print_reproduction(PRF)
  ,m_print_repro_data(p_data)
  ,m_new_hc_size_handler(new_hc_size_handler)
{
  resetBadEpochsCounter();
  setEventHandler();        // initialize handler pointer

  if (nc != 0) m_nc.push_back(nc);              // set first neurocontroller

  if (m_new_hc_size_handler!=0)
  {
    m_new_hc_size_handler(getHypercubeInitialSize(),m_delta, 
      getHypercubeHighLimit());
  }

  m_partial_solution=CIndividual(this,0.0); // initial partial solution
//  m_evalfunc->evaluate(m_partial_solution); // first evaluation
//  m_partial_solution_fitness=m_partial_solution.getFitnessObject()->getMixedFitness();
  m_partial_solution_fitness = MIN_FITNESS;

  m_best_individual_epoch=CIndividual(this,0.0);
  m_best_individual_ever=CIndividual(this,0.0);

  m_operator=new CEAOperator(this);
  Assure(m_operator, "CEA::CEA()");
} // CEA::CEA(void)


//------------------------------------------------------------------------------
void CEA::setLimits(
   const double p_pop_conv_fbc_not_met
  ,const double p_pop_conv_fbc_met
  ,const double p_hyp_large_limit
  ,const double p_final_bound_cond_met_limit
  ,const double p_epoch_conv_limit_final_bound_cond_met
  ,const double p_epoch_conv_limit_final_bound_cond_not_met)
{
  m_pop_conv_fbc_not_met = p_pop_conv_fbc_not_met;
  m_pop_conv_fbc_met = p_pop_conv_fbc_met;
  setHypercubeHighLimit(p_hyp_large_limit);
  m_final_bound_cond_met_limit = p_final_bound_cond_met_limit;
  m_epoch_conv_limit_final_bound_cond_met = p_epoch_conv_limit_final_bound_cond_met;
  m_epoch_conv_limit_final_bound_cond_not_met = p_epoch_conv_limit_final_bound_cond_not_met;
} // void CEA::setLimits(...)


//------------------------------------------------------------------------------
CEA::~CEA(void)
{
  delete m_observer;
  SafeDelete(m_operator);
} // CEA::~CEA(void)


//------------------------------------------------------------------------------
void CEA::dump(std::string file_name)
{
  std::ofstream dump_file(file_name.c_str());
  if (dump_file.is_open() == true)
  {
    dump_file << "m_chromosome_length:"<< getChromosomeLength() << std::endl;
    dump_file << "m_population_size:"<< getPopulationSize()<< std::endl;
    dump_file << "m_reproduction:"<< m_reproduction<< std::endl;
    dump_file << "m_sss_epoch:"<< m_sss_epoch<< std::endl;
    dump_file <<  "m_epoch:"<< m_epoch << std::endl;
    dump_file <<  "m_initial_hypercube:"<< getHypercubeInitialSize()<< std::endl;
    dump_file <<  "m_delta:"<< m_delta<< std::endl;
    dump_file <<  "m_delta_best:"<< m_delta_best<< std::endl;
    dump_file <<  "m_c_h:"<< getHypercubeShrinkFactor()<< std::endl;
    dump_file <<  "m_vectornorm:"<< getVectorNorm() << std::endl;
    dump_file <<  "m_mutation_probability_per_chromosome:"
      << getChromMutationRate()<< std::endl;
  } // if (dump_file.is_open() == true)
} // void dump(std::string file_name)


//------------------------------------------------------------------------------
// BEGIN FIXED ENC PREPARATION METHODS
void CEA::enablePreparationSmoke(void)
{
  if (m_observer!=0) throw std::runtime_error("ENC_SMOKE_REFUSE_OBSERVER");
  if (m_sss_epoch!=0 || m_epoch!=0)
    throw std::runtime_error("ENC_SMOKE_REFUSE_ALREADY_STARTED");
  m_preparation_smoke_enabled=true;
}

void CEA::checkPreparationPopulation(void) const
{
  if (m_preparation_smoke_enabled && (m_preparation_stop ||
      m_preparation_search_calls!=0 || m_preparation_pairs!=0 || getPopulationSize()!=40))
    throw std::runtime_error("ENC_SMOKE_REFUSE_SECOND_OR_PARTIAL_POPULATION");
}

void CEA::checkPreparationReproduceEntry(void) const
{
  if (m_preparation_smoke_enabled && (m_preparation_stop ||
      m_preparation_search_calls!=40 || m_preparation_pairs!=0))
    throw std::runtime_error("ENC_SMOKE_REFUSE_PAIR_ENTRY");
}

void CEA::recordPreparationSearchCall(void)
{
  if (!m_preparation_smoke_enabled) return;
  if (m_preparation_stop || m_preparation_search_calls>=42)
    throw std::runtime_error("ENC_SMOKE_REFUSE_EXTRA_SEARCH_CALL");
  ++m_preparation_search_calls;
}

void CEA::recordPreparationPair(void)
{
  if (!m_preparation_smoke_enabled) return;
  if (m_preparation_stop || m_preparation_search_calls!=42 || m_preparation_pairs!=0)
    throw std::runtime_error("ENC_SMOKE_REFUSE_INCOMPLETE_PAIR");
  ++m_preparation_pairs;
}

void CEA::requestPreparationStop(void)
{
  if (!m_preparation_smoke_enabled) return;
  if (m_preparation_stop || m_preparation_search_calls!=42 || m_preparation_pairs!=1)
    throw std::runtime_error("ENC_SMOKE_REFUSE_STOP_COUNTS");
  // A genuine user/simulation abort remains authoritative and is never cleared.
  if (m_should_stop!=0 && *m_should_stop) return;
  m_preparation_stop=true;
}
// END FIXED ENC PREPARATION METHODS

void CEA::enableObserver(const std::string& filename,unsigned long seed)
{
  if (m_observer || m_preparation_smoke_enabled || m_sss_epoch || m_epoch ||
      getPopulationSize()!=40 || getChromosomeLength()!=881)
    throw std::runtime_error("ENC_OBSERVER_REFUSE_CONFIGURATION_OR_REENTRY");
  ENCObservation::Config config=ENCObservation::Config::load(filename);
  ENCObservation::require(config.seed==seed,"SEED_MISMATCH");
  m_observer=new ENCObservation::Observer(config);
}

bool CEA::observerStopped(void) const
{ return m_observer && m_observer->stopped(); }

bool CEA::reserveObserverGroup(unsigned long count)
{ return !m_observer || m_observer->reserve(count); }

void CEA::recordObserverEvaluation(void)
{ if(m_observer) m_observer->returnedEvaluation(); }

void CEA::recordObserverBoundary(bool population)
{
  if(!m_observer) return;
  // Copy absolute values, not a newly randomised CIndividual; never call the live saver.
  const CIndividual* best=&m_best_individual_ever;
  double fitness=m_best_individual_ever_fitness;
  if(!best->hasFitnessObject() || m_best_individual_epoch_fitness>fitness) {
    best=&m_best_individual_epoch; fitness=m_best_individual_epoch_fitness;
  }
  ENCObservation::require(best->hasFitnessObject(),"BEST_UNAVAILABLE");
  std::vector<double> values;
  for(int i=0;i<best->getLength();++i) values.push_back(best->getSolution().getAllele(i));
  m_observer->boundary(population,values,fitness,m_sss_epoch,m_epoch,m_delta,
    population ? 1.0 : m_delta_factor);
}

void CEA::finishObserver(void)
{
  if(!m_observer) return;
  const CIndividual* best=&m_best_individual_ever;
  double fitness=m_best_individual_ever_fitness;
  if(!best->hasFitnessObject() || m_best_individual_epoch_fitness>fitness) {
    best=&m_best_individual_epoch; fitness=m_best_individual_epoch_fitness;
  }
  ENCObservation::require(best->hasFitnessObject(),"BEST_UNAVAILABLE");
  std::vector<double> values;
  for(int i=0;i<best->getLength();++i) values.push_back(best->getSolution().getAllele(i));
  m_observer->finish(m_should_stop && *m_should_stop,m_ea_converged,values,fitness,
    m_sss_epoch,m_epoch,m_delta,m_observer->completedPairs() ? m_delta_factor : 1.0);
}

void CEA::Optimize(const string& filename /* = "" */)
{
  if(m_observer && !filename.empty()) throw std::runtime_error("ENC_OBSERVER_REFUSE_WARMSTART");
  m_preparation_search_calls=0;
  m_preparation_pairs=0;
  m_preparation_stop=false;
  if (m_preparation_smoke_enabled && (filename.size()!=0 ||
      getPopulationSize()!=40 || getMaxSSSEpochsCount()!=1 || m_sss_epoch!=0 || m_epoch!=0))
    throw std::runtime_error("ENC_SMOKE_REFUSE_WARM_OR_CHANGED_CONFIGURATION");
  if (filename.size() != 0)
  {
    m_partial_solution=CIndividual(this, filename.c_str());
    m_evalfunc->evaluate(m_partial_solution, false);
    m_partial_solution_fitness = 
      m_partial_solution.getFitnessObject()->getMixedFitness();

    m_evalfunc->evaluate(m_best_individual_epoch, false);
    m_best_individual_epoch_fitness = 
      m_best_individual_epoch.getFitnessObject()->getMixedFitness();

    m_evalfunc->evaluate(m_best_individual_ever, false);
    m_best_individual_ever_fitness = 
      m_best_individual_ever.getFitnessObject()->getMixedFitness();
  } // if (filename.size() != 0)
  else
  {
    SearchSpaceScan();
    if (!m_preparation_stop && !observerStopped()) modifyHypercube();
  } // if (filename.size() != 0) else branch
  if (!m_preparation_stop && !observerStopped()) Run();
  finishObserver();
} // void Optimize(const string& filename /* = "" */)


//------------------------------------------------------------------------------
void CEA::modifyHypercube(void)
{
  m_delta*=getHypercubeShrinkFactor();
  callEventHandler(HYPERCUBE_SIZE_CHANGE);
} // void CEA::shrinkHypercube(void)


//------------------------------------------------------------------------------
void CEA::SearchSpaceScan(void) 
{
  m_delta=getHypercubeInitialSize();

  // for explanation of the following check see CEA::Run()
  while (getMaxSSSEpochsCount() > m_sss_epoch
    && (m_should_stop!=0?(*m_should_stop):false) == false
    && m_preparation_stop == false && !observerStopped())
  {
    if(!reserveObserverGroup(getPopulationSize())) break;
    ++m_sss_epoch;
    // print search space scan epoch information
    ostringstream mess("");mess<<"\n... search space scan epoch ";
    mess<< m_sss_epoch <<" of "<< getMaxSSSEpochsCount() 
      <<" [H = " <<setprecision(9)<<m_delta<<"]: \n";
    statusMessage(mess);

    // create population and start new reproduction epoch
    m_population = new CPopulation(
      this
      ,getPopulationSize()
      ,m_rep
      ,m_print_reproduction
      ,m_print_repro_data
      );
    Assure(m_population, 
     "CEA::SearchSpaceScan() => Error during allocation of CPopulation object!");
    m_population->setStopBoolPointer(m_should_stop);
    callEventHandler(NEW_POPULATION_INITIALIZED);
    updateStatusObject();
    recordObserverBoundary(true);

    // Start reproducing
    m_population->Reproduce();

    // now m_best_individual_epoch holds the best chromosome found in 
    // this epoch, so set new best individual ever
    if (m_best_individual_epoch_fitness > m_best_individual_ever_fitness) 
    {
      callEventHandler(SUCCESSFUL_SSS);
      m_best_individual_ever=m_best_individual_epoch;
      m_best_individual_ever_fitness=m_best_individual_epoch_fitness;
    }
    else
    {
      callEventHandler(UNSUCCESSFUL_SSS);
    }

    // set up status message string and give out
    mess.str("");mess<<setiosflags(ios::fixed)<<setprecision(9);
    mess<<"best fitness within this epoch: " << m_best_individual_epoch_fitness;
    mess<<"\nbest fitness so far: " <<m_best_individual_ever_fitness<<endl;
    mess<<resetiosflags(ios::fixed);statusMessage(mess);

    SafeDelete(m_population); // destroy old population
    // m_operator->PrintCredits(); // optional, used not at the moment
  } // search space scan epoch loop

  m_partial_solution=m_best_individual_ever;
  m_partial_solution.AddChromosomeToPartialSolution();
  m_partial_solution_fitness=m_best_individual_ever_fitness;    

  // now m_partial_solution (and m_best_individual_ever) hold the best 
  // individual found by search space scanning of populations
} // SearchSpaceScan()


//------------------------------------------------------------------------------
void CEA::goodEpochPostprocess(void) throw()
{
  callEventHandler(SUCCESSFUL_EPOCH);
  // save best individual of epoch as new best individual ever
  m_best_individual_ever=m_best_individual_epoch;
  m_best_individual_ever_fitness=m_best_individual_epoch_fitness;

  // use chromosome of best individual of epoch 
  // as partial solution for new population
  m_partial_solution=m_best_individual_epoch;
  m_partial_solution.AddChromosomeToPartialSolution();
  m_partial_solution_fitness=m_best_individual_epoch_fitness;

  m_delta_best=m_delta;   // save new best hypercube size
  modifyHypercube();
  resetBadEpochsCounter(); 
} // void CEA::goodEpochPostprocess(void)


//------------------------------------------------------------------------------
void CEA::badEpochPostprocess(void) throw()
{
  callEventHandler(UNSUCCESSFUL_EPOCH);
  // reinitialize partial solution with best individual ever
  m_partial_solution=m_best_individual_ever;
  m_partial_solution_fitness=m_best_individual_ever_fitness;

  modifyHypercube();
  resetBadEpochsCounter();
} // void CEA::badEpochPostprocess(void)


//------------------------------------------------------------------------------
void CEA::lastBadEpochPostprocess(void) throw()
{
  callEventHandler(LAST_UNSUCCESSFUL_EPOCH);
  // re-use chromosome of best individual of last epoch
  // as partial solution for initializatin of new population
  m_partial_solution=m_best_individual_epoch;
  m_partial_solution.AddChromosomeToPartialSolution();
  m_partial_solution_fitness=m_best_individual_epoch_fitness;
  incBadEpochsCounter();  // count bad epochs
} // void CEA::lastBadEpochPostprocess(void)


//----------------------------------------------------------------------------
void CEA::Run(void) 
{
  CEA::resetBadEpochsCounter();
  // shrinking factor * start size  
  //m_delta = getHypercubeShrinkFactor() * getHypercubeStartSize();

  // The following loop condition seems a little bit complicated but it is not
  // The first check is, whether convergence has already been archieved or not
  // and the second one checks a flag. Because this is flag is usually provided
  // by other applications, it's done with a pointer to a boolean. Thus you've
  // got to check the pointer for validity and then check the value itself.
  // This is done with a C-style trigraph.
  while ( (m_ea_converged == false) &&         // set in CPopulation::Reproduce
    (m_should_stop!=0?(*m_should_stop):false) ==false
    && m_preparation_stop == false && !observerStopped())
  {
    if(!reserveObserverGroup(getPopulationSize())) break;
    if (m_new_hc_size_handler!=0)
    {
      m_new_hc_size_handler(getHypercubeInitialSize(),m_delta, 
        m_data.getHypercubeUpperLimit());
    } // if (m_new_hc_size_handler!=0)

    // aborting of algorithm by setting flag to signal premature convergence
    if ((m_should_stop != 0) && ((*m_should_stop) == true)) m_ea_converged=true;

    ++m_epoch;                                         // increase epoch counter

    // set up epoch information string
    ostringstream mess("");            // message string stream
    mess << "\n... epoch "<<m_epoch<<" [H = " <<m_delta;
    if (getBadEpochsCount()==0) mess << "]:\n"; 
    else mess<<" after "<<getBadEpochsCount() << " unsuccessful epoch(s)]: \n";

    // Output through global output object and report file
    statusMessage(mess);

    // save best fitness from last epoch
    double best_individual_ever_of_last_epoch_fitness(m_best_individual_ever_fitness);

    // create new poulation
    m_population = new CPopulation(
      this
      ,getPopulationSize()
      ,m_rep
      ,m_print_reproduction
      ,m_print_repro_data
      );
    Assure(m_population, "CEA::Run()=> Failure during allocation of CPopulation object!");
    
    m_population->setStopBoolPointer(m_should_stop);// and set stop flag pointer
    callEventHandler(NEW_POPULATION_INITIALIZED);
    updateStatusObject();
    recordObserverBoundary(true);

    // optimize through reproduction and evaluating
    double epoch_improvement = m_population->Reproduce();
    // now m_best_individual_epoch holds the best individual found in this epoch

    if (m_best_individual_epoch_fitness > m_best_individual_ever_fitness)
      goodEpochPostprocess();
    else if (getBadEpochsCount() >= getMaxBadEpochsCount())
      badEpochPostprocess();
    else 
      lastBadEpochPostprocess();

    updateStatusObject();
    // Output some statistics after each epoch. This is done
    // regardless of failure (bad epoch) or success (new best individual).
    mess.seekp(0,ios_base::beg);
    mess<<scientific<<setprecision(9);
    mess<<"best fitness within this epoch: "<<m_best_individual_epoch_fitness;
    mess<<"\nbest fitness so far: "<<m_best_individual_ever_fitness<<endl;
    mess<<"absolute improvement w.r.t. last epoch: "
      << (m_best_individual_epoch_fitness-best_individual_ever_of_last_epoch_fitness);
    mess<<"\nrelative improvement within this epoch: " <<epoch_improvement<<endl;
    mess<<resetiosflags(ios::scientific);statusMessage(mess);

    // delete old population (is easier and safer than re-initializing the old 
    // one)
    SafeDelete(m_population);

    // m_operator->PrintCredits(); // optional, not used at the moment
  } // end of "while (!m_ea_converged) {" - scope

  m_partial_solution=m_best_individual_ever;
  m_partial_solution.AddChromosomeToPartialSolution();
  m_partial_solution_fitness=m_best_individual_ever_fitness;
  // Now m_best_individual_ever holds the best individual found by the
  // evolutionary algorithm. This can be saved to file using the 
  // SaveBestSolution() method.
} // void CEA::Run(void)


//------------------------------------------------------------------------------
void CEA::SaveBestSolution(const char* filename)
{
  if (m_best_individual_ever.getFitnessObject()==0)
  {
    if (m_best_individual_epoch.getFitnessObject()!=0)
    {
      m_best_individual_ever = m_best_individual_epoch;
      m_best_individual_ever_fitness = m_best_individual_epoch_fitness;
    } // if (m_best_individual_epoch.getFitnessObject()!=0)
  } // if (m_best_individual_ever.getFitnessObject()==0)

  // Progress message
  Output(string("\n... saving best solution to file\n"));

  // Combine partial solution and chromosome to solution.
  m_best_individual_ever.AddChromosomeToPartialSolution();
  m_best_individual_ever.Save(filename);

} // void CEA::SaveBestSolution(const char*filename)
