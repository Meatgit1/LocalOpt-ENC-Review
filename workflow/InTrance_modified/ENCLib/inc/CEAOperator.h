/*!
\file CEAOperator.h
\brief A class that provides the functionality for evolutionary operators.

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

  09.05.2006 (AO)

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
#ifndef CEAOPERATOR_H_INCLUDED
#define CEAOPERATOR_H_INCLUDED

// SpaceLib
#include <MLGlobal.h>

// ENCLib
#include <configENC.h>


static const int N_EVOPERATORS = 4;


enum EvOperator {
  ONE_POINT_CROSSOVER    = 0,
  UNIFORM_CROSSOVER      = 1,
  ARITHMETICAL_CROSSOVER = 2,
  CROSSOVER_NODES        = 3
};


// evolutionary operator class declaration
// ----------------------------------------------------------------------------
class CEAOperator
{
  //-----------------------------------------------------------------
  CEA*   m_ea;
  int    m_chromosome_length;

  int    m_number_of_evoperators;
  double m_probability[N_EVOPERATORS];

  double m_credit[N_EVOPERATORS];
  //-----------------------------------------------------------------

  // ~~~~~ begin inline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //----------------------------------------------------------------------------
  void OnePointCrossover(std::vector<CNC*>& nc) {
    CIndividual* ch1 = m_ea->m_population->getFirstChild(); // need for speed
    CIndividual* ch2 = m_ea->m_population->getSecondChild();
    CIndividual* pa1 = m_ea->m_population->getFirstParent();
    CIndividual* pa2 = m_ea->m_population->getSecondParent();

    const int crossover_locus = MathLib::randIntInterval(0, m_chromosome_length-1);

    for(int i = 0; i < m_chromosome_length; ++i) 
    {
      if (i <= crossover_locus) 
      {
        (ch1->m_chromosome)[i] = (pa1->m_chromosome)[i];
        (ch2->m_chromosome)[i] = (pa2->m_chromosome)[i];
      }
      else 
      {
        (ch1->m_chromosome)[i] = (pa2->m_chromosome)[i];
        (ch2->m_chromosome)[i] = (pa1->m_chromosome)[i];
      }
    }
    ch1->UpdateSolution();
    ch2->UpdateSolution();
  };


  //----------------------------------------------------------------------------
  void UniformCrossover(std::vector<CNC*>& nc) 
  {
    CIndividual* ch1 = m_ea->m_population->getFirstChild(); // need for speed
    CIndividual* ch2 = m_ea->m_population->getSecondChild();
    CIndividual* pa1 = m_ea->m_population->getFirstParent();
    CIndividual* pa2 = m_ea->m_population->getSecondParent();

    for(int i = 0; i < m_chromosome_length; ++i) 
    {
      if (MathLib::FlipCoin()) 
      {
        (ch1->m_chromosome)[i] = (pa1->m_chromosome)[i];
        (ch2->m_chromosome)[i] = (pa2->m_chromosome)[i];
      }
      else 
      {
        (ch1->m_chromosome)[i] = (pa2->m_chromosome)[i];
        (ch2->m_chromosome)[i] = (pa1->m_chromosome)[i];
      }
    }
    ch1->UpdateSolution();
    ch2->UpdateSolution();
  };


  //----------------------------------------------------------------------------
  void ArithmeticalCrossover(std::vector<CNC*>& nc) {
    CIndividual* ch1 = m_ea->m_population->getFirstChild(); // need for speed
    CIndividual* ch2 = m_ea->m_population->getSecondChild();
    CIndividual* pa1 = m_ea->m_population->getFirstParent();
    CIndividual* pa2 = m_ea->m_population->getSecondParent();

    const double a = MathLib::rand01();

    for (int i = 0; i < m_chromosome_length; ++i)
    {
      (ch1->m_chromosome)[i] = a*(pa1->m_chromosome)[i] + (1.0-a)*(pa2->m_chromosome)[i];
      (ch2->m_chromosome)[i] = a*(pa2->m_chromosome)[i] + (1.0-a)*(pa1->m_chromosome)[i];
    }
    ch1->UpdateSolution();
    ch2->UpdateSolution();
  };


  //----------------------------------------------------------------------------
  // Node cossover operator had to be adopted for the special application
  // of chromosomes containing the net occupation of at least two neuronal
  // nets.
  void CrossoverNodes(std::vector<CNC*>& nc) 
  {
    CIndividual* ch1 = m_ea->m_population->getFirstChild(); // need for speed
    CIndividual* ch2 = m_ea->m_population->getSecondChild();
    CIndividual* pa1 = m_ea->m_population->getFirstParent();
    CIndividual* pa2 = m_ea->m_population->getSecondParent();

    bool coin, next_is_same_neuron = true;

    int i_new_nc(0);// index from where a new nc is to use for node estimation
    int i_last_nc(0);// start index in chromosome of actual nc's net occ.
    int i_nc(-1); // nc index (actual used nc)
    CNC* l_nc(0); // local copy of nc object pointer (from parameter vector)

    for (int i = 0; i < m_chromosome_length; ++i) 
    {

      // this is new  (AO)
      //------------------
      if (i == i_new_nc)  // do we have to use the next controller
      {
        ++i_nc;                             // increase nc index
        if (i_nc < int(nc.size()))    // is there still a nc pointer left?
        {
          l_nc = nc[i_nc];                // set local vaiable
          next_is_same_neuron = true;  // this always holds true for new ncs
        } else
        {
          l_nc = 0;                    // there's no nc left to use
          next_is_same_neuron = false; // the rest must be special parameters and not nodes
        }
        if (l_nc != 0)    // there's is still at least one nc
        {
          i_last_nc = i;  // set start index in chromosome of new nc data
          i_new_nc += l_nc->getNumberOfParameters(); // increase nc data end index in chromosome
        }
      }

      if (i == 0 || !next_is_same_neuron) coin = MathLib::FlipCoin();

      if (coin) 
      {
        (ch1->m_chromosome)[i] = (pa1->m_chromosome)[i];
        (ch2->m_chromosome)[i] = (pa2->m_chromosome)[i];
      } else 
      {
        (ch1->m_chromosome)[i] = (pa2->m_chromosome)[i];
        (ch2->m_chromosome)[i] = (pa1->m_chromosome)[i];
      }

      // the following is only done in case of valid nc pointer
      //if ( (i < m_chromosome_length-1)  && (l_nc!=0) ) 
      //{
      //  if (l_nc->m_parameter_neuron_number[i-i_last_nc  ] 
      //   == l_nc->m_parameter_neuron_number[i-i_last_nc+1])
      //    next_is_same_neuron = true;
      //  else next_is_same_neuron = false;
      //}
      if ( (i < m_chromosome_length-1)  && (l_nc!=0) ) 
      {
        if (l_nc->getParameterNeuronNumber(i-i_last_nc)
         == l_nc->getParameterNeuronNumber(i-i_last_nc+1))
          next_is_same_neuron = true;
        else next_is_same_neuron = false;
      }
    }
    ch1->UpdateSolution();
    ch2->UpdateSolution();
  };

  // ~~~~~ end inline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // pointer array to evolutionary operators

  void (CEAOperator::*evop_ptr[N_EVOPERATORS])(std::vector<CNC*>& nc);

//=============================================================================
// CEAOperator CLASS INTERFACE
//=============================================================================
public:

  // ~~~~~ begin inline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  CEAOperator(CEA* ea)
  : m_ea(ea),
    m_chromosome_length(ea->getChromosomeLength()),
    m_number_of_evoperators(N_EVOPERATORS)
  {
    setOperatorProbability(ONE_POINT_CROSSOVER,    0.25);
    setOperatorProbability(UNIFORM_CROSSOVER,      0.25);
    setOperatorProbability(ARITHMETICAL_CROSSOVER, 0.0);
    setOperatorProbability(CROSSOVER_NODES,        0.5);

    m_credit[ONE_POINT_CROSSOVER]    = 0.0;
    m_credit[UNIFORM_CROSSOVER]      = 0.0;
    m_credit[ARITHMETICAL_CROSSOVER] = 0.0;
    m_credit[CROSSOVER_NODES]        = 0.0;

    evop_ptr[ONE_POINT_CROSSOVER]    = &CEAOperator::OnePointCrossover;
    evop_ptr[UNIFORM_CROSSOVER]      = &CEAOperator::UniformCrossover;
    evop_ptr[ARITHMETICAL_CROSSOVER] = &CEAOperator::ArithmeticalCrossover;
    evop_ptr[CROSSOVER_NODES]        = &CEAOperator::CrossoverNodes;
  };

  //----------------------------------------------------------------------------
  void setOperatorProbability(const EvOperator op,
                              double prob)
  {
    m_probability[op] = prob;
  };


  //----------------------------------------------------------------------------
  void AssignCredit(const EvOperator op
                   ,      double     credit) 
  {
    m_credit[op] += credit;
  };


  //----------------------------------------------------------------------------
  double Credit(const EvOperator op) const 
  {
    return m_credit[op];
  };


  //----------------------------------------------------------------------------
  void PrintCredits(void) const 
  {
    std::cout << std::endl;
    std::cout << m_credit[ONE_POINT_CROSSOVER]    << " / ";
    std::cout << m_credit[UNIFORM_CROSSOVER]      << " / ";
    std::cout << m_credit[ARITHMETICAL_CROSSOVER] << " / ";
    std::cout << m_credit[CROSSOVER_NODES]        << std::endl;
  };

  //----------------------------------------------------------------------------
  EvOperator Operator(std::vector<CNC*>& nc) 
  {
    const double spin = MathLib::rand01ex();
    double partsum = 0.0;
    int op = 0;

    for (int i = 0; partsum <= spin; ++i) 
    {
      partsum += m_probability[i];
      op = i;
    }

#ifdef development_version
    if(op < 0 || op >= m_number_of_evoperators) FatalError ("CEAOperator::Operator(): index out of range");
#endif

    (this->*evop_ptr[op])(nc);

    return static_cast<EvOperator>(op);
  };

  // ~~~~~ end inline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

}; // end class CEAOperator


#endif // CEAOPERATOR_H_INCLUDED
