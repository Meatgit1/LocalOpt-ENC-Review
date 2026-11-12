/*!
\file CPopulation.cpp
\brief A class that provides the functionality for a population of an evolutionary algorithm.

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

  NU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++

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

// C standard library
#include <math.h>
#include <conio.h>

//STL
#include <iomanip>
#include <string>
#include <sstream>
#include <cassert>
#include <algorithm>

//SpaceLib
#include <Global.h>
#include <MLGlobal.h>


//ENCLib
#include <configENC.h>
#include <CNC.h>
#include <CPopulationConvergence.h>
#include <CEAFitness.h>
#include <CEA.h>
#include <CIndividual.h>
#include <CPopulation.h>
#include <CEAOperator.h>
#include <CEAEvalFunc.h>


using namespace std;
using namespace Global;
using namespace MathLib;


static unsigned long int repro_ctr(0);

//-----------------------------------------------------------------------------
void ResetReproCtr(void) 
{
  repro_ctr=0;
}

//------------------------------------------------------------------------------
CIndividual *CPopulation::getRandomIndividual(void)
{
  return m_individual[randIntInterval(0, m_size-1)];
} // CIndividual *CPopulation::getRandomIndividual(void)

//-----------------------------------------------------------------------------
void CPopulation::modifyHypercubeOnGoodEpochReproduction(void)
{
  m_ea->m_delta_factor = std::min( 2.0,
    m_ea->m_delta_factor * m_ea->getHypercubeModFactorOnGoodEpochRepro());
} // void CPopulation::modifyHypercubeOnGoodEpochReproduction(void)

//-----------------------------------------------------------------------------
void CPopulation::modifyHypercubeOnBadEpochReproduction(void)
{
  m_ea->m_delta_factor *= m_ea->getHypercubeModFactorOnBadEpochRepro();
} // void CPopulation::modifyHypercubeOnBadEpochReproduction(void)

//-----------------------------------------------------------------------------
void CPopulation::Tournament(CIndividual*  comp1
                            ,CIndividual*  comp2
                            ,CIndividual** winner
                            ,CIndividual** loser
                            ) {
  double f1, f2; // fitness values of competitor one and two

  // tournament according to time and the mixed fitness
  const double r = rand01();  // a (pseudo-) random number out of [0,1]
  const double tp = m_ea->getTimeSelectionPressure();

  if (r < fabs(tp)) 
  { // time fitness
    if (tp < 0.0) 
    {
      f1 = comp1->getFitnessObject()->getFitness(m_time_fitness_index);
      f2 = comp2->getFitnessObject()->getFitness(m_time_fitness_index);
    } // if (tp < 0.0) 
    else 
    {
      f1 = comp1->getFitnessObject()->getFitness(m_tournament_fitness_index);
      f2 = comp2->getFitnessObject()->getFitness(m_tournament_fitness_index);
    } // else branch of if (tp < 0.0)
    if (tp < 0.0) // negative time pressure
      swap(f1, f2);
  } // if (r < fabs(tp)) 

  if (r >= fabs(tp) || f1 == f2) {// mixed fitness, 100% if time/mass fitnesses are equal
    f1 = comp1->getFitnessObject()->getMixedFitness();
    f2 = comp2->getFitnessObject()->getMixedFitness();
  }

  if (f1 > f2) 
  {
    *winner = comp1;
    *loser  = comp2;
    m_fitness_sum-=f2;
  } else if (f2 > f1) 
  {
    *winner = comp2;
    *loser  = comp1;
    m_fitness_sum-=f1;
  } else   // if f1 == f2
  {
    if (FlipCoin()) 
    {
      *winner = comp1;
      *loser  = comp2;
      m_fitness_sum-=f2;
    } 
    else 
    {
      *winner = comp2;
      *loser  = comp1;
      m_fitness_sum-=f1;
    }
  }  
} //void CPopulation::Tournament(CIndividual*  comp1....



//-----------------------------------------------------------------------------
void CPopulation::PrintReproduction( CIndividual* ind) const 
{
  if (PrintReproFunc != 0)
   PrintReproFunc(ind
    ,m_convergence->getImprovement(m_number_of_reproductions
                                  ,m_ea->m_best_individual_epoch_fitness)
    ,m_number_of_reproductions
    ,repro_ctr
    ,m_data);
};


//-----------------------------------------------------------------------------
double CPopulation::Reproduction(void) 
{
  // selection process; first select two individuals
  CIndividual* competitor1 = SelectParent(NULL, NULL, NULL);
  CIndividual* competitor2 = SelectParent(competitor1, NULL, NULL);

  // tournament one to determine parent individual 1 and replaced 1 (will be child 1)
  Tournament(competitor1, competitor2, &m_parent1, &m_child1);

  // select another two individuals
  CIndividual* competitor3 = SelectParent(competitor1, competitor2, NULL);
  CIndividual* competitor4 = SelectParent(competitor1, competitor2, competitor3);

  // tournament two to determine parent individual 2 and replaced 2 (will be child 2)
  Tournament(competitor3, competitor4, &m_parent2, &m_child2);

  // apply evolutionary operator (result in child 1 and 2)
  EvOperator op = m_ea->m_operator->Operator(m_ea->m_nc);


  // mutation of child individual 1
  if (Flip(m_ea->getChromMutationRate()))
    m_child1->OnePointMutation(m_ea->m_delta*m_ea->m_delta_factor);
  else if (m_ea->getGenMutationRate() > 0.0)
    m_child1->GeneralMutation(
    m_ea->m_delta*m_ea->m_delta_factor,m_ea->getGenMutationRate());

  // evaluation of new child individual 1
  m_ea->getEvalFunc()->evaluate(*m_child1);
  m_ea->recordObserverEvaluation();
  m_ea->recordPreparationSearchCall();
  const double child_1_fitness = m_child1->getFitnessObject()->getMixedFitness();
  m_fitness_sum+=child_1_fitness;


  // mutation of child individual 2
  if (Flip(m_ea->getChromMutationRate()))
    m_child2->OnePointMutation(m_ea->m_delta*m_ea->m_delta_factor);
  else if (m_ea->getGenMutationRate() > 0.0)
    m_child2->GeneralMutation(
    m_ea->m_delta*m_ea->m_delta_factor,m_ea->getGenMutationRate());

  // evaluation of new child individual 2
  m_ea->getEvalFunc()->evaluate(*m_child2);
  m_ea->recordObserverEvaluation();
  m_ea->recordPreparationSearchCall();
  const double child_2_fitness = m_child2->getFitnessObject()->getMixedFitness();
  m_fitness_sum+=child_2_fitness;


  // Is any of the two new child individuals better than the best population
  // individual so far (fitness value comparison)?
  if ((child_1_fitness > m_ea->m_best_individual_epoch_fitness) ||
      (child_2_fitness > m_ea->m_best_individual_epoch_fitness))
  {
    if (child_1_fitness < child_2_fitness) // child 2 is best of two children
    {
      m_ea->m_best_individual_epoch = *m_child2;
      m_ea->m_best_individual_epoch_fitness = child_2_fitness;
      PrintReproduction(m_child2);
    }             
    else // child 1 is better than child 2
    {
      m_ea->m_best_individual_epoch = *m_child1;
      m_ea->m_best_individual_epoch_fitness = child_1_fitness;
      PrintReproduction(m_child1);
    }
    modifyHypercubeOnGoodEpochReproduction();
  }
  else modifyHypercubeOnBadEpochReproduction();

  // store best population individual fitness so far for evolutionary
  // algorithm convergence determination
  m_convergence->setFitness(m_number_of_reproductions
    , m_ea->m_best_individual_epoch_fitness);

  ++repro_ctr;  // increase global reproductions counter
  m_ea->recordPreparationPair();

  return m_ea->m_best_individual_epoch_fitness;
};



//----------------------------------------------------------------------------
CPopulation::CPopulation( CEA* ea
                         ,const int popsize
                         ,const string& rep_file
                         ,PrintReproductionFunction PRF
                         ,void* data
                        )
: m_child1(0)
, m_child2(0),
  m_ea(ea),
  m_rep(rep_file),
  PrintReproFunc(PRF),
  m_data(data),
  m_tournament_fitness_index(0),
  m_time_fitness_index(0),
  m_size(popsize),
  m_fitness_sum(0.0),
  m_convergence(0),
  m_number_of_reproductions(0),
  m_parent1(NULL),
  m_parent2(NULL),
  m_should_stop(0)
{

  m_ea->checkPreparationPopulation();

  // calculate number of fitness values to hold for convergence
  // estimation (n_conv=20.0*n_pop^0.8)
  const int conv_size(static_cast<int>(20.0*pow(double(popsize), 0.8)));

  // initialize new convergence object
  m_convergence = new CPopulationConvergence(
    conv_size
   ,m_ea->m_pop_conv_fbc_not_met
   ,m_ea->m_pop_conv_fbc_met
   ,m_ea->m_final_bound_cond_met_limit
  );
  Assure(m_convergence, "CPopulation::CPopulation()");

  // create and evaluate individuals
  m_individual = new CIndividual*[popsize];
                 Assure(m_individual, "CPopulation::CPopulation()");

  // initialize fitness of best individual of epoch to a very bad one
  m_ea->m_best_individual_epoch_fitness = MIN_FITNESS;

  CEAEvalFunc* evaluator = m_ea->getEvalFunc();
  double act_fitness; // fitness value of new created individuals

  for (int i = 0; i < popsize; ++i) 
  {
    m_individual[i] = new CIndividual(m_ea, m_ea->m_partial_solution);
                      Assure(m_individual[i], "CPopulation::CPopulation()");

    // initial evaluation of individual of new population
    evaluator->evaluate(*(m_individual[i]));
    m_ea->recordObserverEvaluation();
    m_ea->recordPreparationSearchCall();
    act_fitness = std::max(MIN_FITNESS, m_individual[i]->getFitnessObject()->getMixedFitness());
    m_fitness_sum+=act_fitness;

    // if one of the initial population's inidviduals is better than
    // the currently saved best individual of epoch make this new individual
    // the new best inidividual of epoch
    if (act_fitness >= m_ea->m_best_individual_epoch_fitness) 
    {
      m_ea->m_best_individual_epoch=*(m_individual[i]);
      assert(m_ea->m_best_individual_epoch.getFitnessObject()!=0);
      m_ea->m_best_individual_epoch_fitness = act_fitness;
    }
    assert(m_ea->m_best_individual_epoch.getFitnessObject()!=0);
  }
  assert(m_ea->m_best_individual_epoch.getFitnessObject()!=0);
  if (m_ea->hasStatusObject()==true)
  {
    assert(m_ea->m_best_individual_epoch.getFitnessObject()!=0);
    m_ea->getStatusObject()->avg_fitness = float(m_fitness_sum/m_size);
  }    
  PrintReproduction(&(m_ea->m_best_individual_epoch));
};



//----------------------------------------------------------------------------
CPopulation::~CPopulation() 
{
  for (int i = 0; i < m_size; ++i) SafeDelete(m_individual[i]);
  SafeDeleteArray(m_individual);
  SafeDelete(m_convergence);
};


//----------------------------------------------------------------------------
CIndividual* CPopulation::SelectParent(const CIndividual* const not_this1,
                                       const CIndividual* const not_this2,
                                       const CIndividual* const not_this3) {
  CIndividual* parent = NULL;
  do {
    parent = getRandomIndividual();
  } while (parent == not_this1 || parent == not_this2 || parent == not_this3);
  return parent;
};

//-----------------------------------------------------------------------------
void CPopulation::checkKeyboard(void)
{
  if (_kbhit()!=0) 
  {
    bool extended = false;
    int first = _getch();
    int second = 0;
    if (first==0)
    {
      extended = true;
      second = _getch();
    }
    cout << endl << first << " " << second << " (" << (extended==true?"true":"false") << ")\n";
    KeyboardHandler kbh = m_ea->getKeyboardHandler();
    if (kbh!=0) kbh(first,second);
  }
}


//-----------------------------------------------------------------------------
double CPopulation::Reproduce(void) 
{
  // A population can itself exhaust the count budget; do not enter a pair then.
  if(m_ea->observerStopped()) return 0.0;
  m_ea->checkPreparationReproduceEntry();
  m_convergence->Reset();

  double act_fitness;
  double first_fitness = std::max(MIN_FITNESS, m_ea->m_best_individual_epoch_fitness);
  double best_fitness(first_fitness);
  double avg_fitness(0.0);
  convergence = false;
  
  if (m_ea->hasStatusObject())
    m_ea->getStatusObject()->gen_mutation_probability =
    float(m_ea->getGenMutationRate());
  m_ea->m_delta_factor = 1.0;

  do
  {
    // Reserve both evaluations before selection/operator/mutation consumes RNG.
    if(!m_ea->reserveObserverGroup(2)) break;
    checkKeyboard();

    act_fitness  = std::max(Reproduction(), MIN_FITNESS);

    avg_fitness = m_fitness_sum/m_size;

    if (act_fitness > best_fitness)
    {
      best_fitness = act_fitness;

      if (m_ea->hasStatusObject() == true)
      {
        m_ea->getStatusObject()->current_best_fitness=float(act_fitness);
        m_ea->getStatusObject()->avg_fitness=float(avg_fitness);
        m_ea->getStatusObject()->mutation_probability = 
          float(m_ea->getChromMutationRate());
        m_ea->getStatusObject()->hc_size = 
          float(m_ea->m_delta*m_ea->m_delta_factor);
      }

    } // if ( act_fitness > best_fitness)

    if (m_ea->getAvgFitnessConvergence()==false && (convergence==false))
      convergence = m_convergence->Convergence(m_number_of_reproductions++
      ,act_fitness);
    else
      if ((best_fitness-avg_fitness)/best_fitness <= 
        m_ea->getTriggerAvgFitnessRatio()) convergence = true;
    if (m_ea->getDynamicMutation()==true) 
      m_ea->setChromMutationRate(
      bounded(1.0-(best_fitness-avg_fitness)/best_fitness,0.0,1.0));

    // Stop only after both children and the native convergence/mutation updates.
    m_ea->requestPreparationStop();
    m_ea->recordObserverBoundary(false);
  } while ( (convergence==false) && (abort()==false)
    && (m_ea->preparationStopRequested()==false) && !m_ea->observerStopped() );

  double improvement = (best_fitness-first_fitness);
  if (first_fitness!=0.0) improvement/=fabs(first_fitness);
  else improvement = 1.0;

  // Budget exhaustion is not scientific convergence. Genuine abort handling below stays intact.
  if (!m_ea->preparationStopRequested() && !m_ea->observerStopped())
  {
  if (first_fitness > m_ea->m_final_bound_cond_met_limit)
  {
    // final boundary conditions are met at beginning of epoch
    if (improvement < m_ea->m_epoch_conv_limit_final_bound_cond_met)
      m_ea->m_ea_converged = true;
  }
  else
  {
    // final boundary conditions are not met at beginning of epoch
    if (improvement < m_ea->m_epoch_conv_limit_final_bound_cond_not_met)
      m_ea->m_ea_converged = true;
  }

  // no convergence, if hypercube is still large
  if (m_ea->m_delta > m_ea->getHypercubeHighLimit() 
    || (m_ea->getBadEpochsCount() != 0) ) m_ea->m_ea_converged = false; 

  } // unchanged convergence logic only when not stopped by the preparation probe

  if (abort() == true) m_ea->m_ea_converged = true;

  return improvement;
};



//-----------------------------------------------------------------------------
void CPopulation::setStopBoolPointer(bool* should_stop_pointer)
{
  m_should_stop = should_stop_pointer;
}


//-----------------------------------------------------------------------------
bool* CPopulation::getStopBoolPointer(void)
{
  return m_should_stop;
}
