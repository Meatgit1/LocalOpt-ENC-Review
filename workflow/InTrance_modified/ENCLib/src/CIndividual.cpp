/*!
\file CIndividual.cpp
\brief A class that provides the functionality for an individual of an evolutionary algorithm.

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

  (c) 2001-2004 by Bernd Dachwald, Andreas Ohndorf

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/

// STL
#include <iomanip>

// SpaceLib
#include <Global.h>

// ENCLib
#include <CEAEvalFunc.h>
#include <CEAFitness.h>
#include <CIndividual.h>
#include <CEA.h>


// namespaces
using namespace std;
using namespace Global;
using namespace MathLib;


//------------------------------------------------------------------------------
CIndividual::CIndividual(void)
: m_length(0)
, m_initial_hypercube_size(0.0)
, m_fitness(0)
, m_chromosome()
, m_partial_solution()
, m_solution() {}


//------------------------------------------------------------------------------
CIndividual::CIndividual(CEA* ea)
: 
  m_length(ea->getChromosomeLength())
, m_initial_hypercube_size(ea->getHypercubeInitialSize())
, m_fitness(0)
, m_chromosome(m_length, -ea->m_delta, ea->m_delta)
, m_partial_solution(m_length)
, m_solution(m_length)
{
  UpdateSolution();
};


//------------------------------------------------------------------------------
CIndividual::CIndividual(CEA* ea, const double value)
: m_length(ea->getChromosomeLength())
, m_initial_hypercube_size(ea->getHypercubeInitialSize())
, m_fitness(0)
, m_chromosome(m_length, value, value)
, m_partial_solution(m_length)
, m_solution(m_length)
{
  UpdateSolution();
}; // CIndividual::CIndividual(CEA* ea, const double value)


//------------------------------------------------------------------------------
CIndividual::CIndividual(CEA* ea, const CIndividual& ps)
: m_length(ea->getChromosomeLength())
, m_initial_hypercube_size(ea->getHypercubeInitialSize()) 
, m_fitness(ps.hasFitnessObject()==true?ps.m_fitness->Duplicate():0)
, m_chromosome(m_length, -ea->m_delta, ea->m_delta)
, m_partial_solution(ps.m_solution)
, m_solution(m_length)
{  
  UpdateSolution();
};


//------------------------------------------------------------------------------
CIndividual::CIndividual(CEA* ea, const double value, const CIndividual& ps)
: m_length(ea->getChromosomeLength())
, m_initial_hypercube_size(ea->getHypercubeInitialSize()) 
, m_fitness(ps.hasFitnessObject()==true?ps.m_fitness->Duplicate():0)
, m_chromosome(m_length, value, value)
, m_partial_solution(ps.m_solution)
, m_solution(m_length)
{
  UpdateSolution();
};


//------------------------------------------------------------------------------
CIndividual::CIndividual(CEA* ea, const char* filename)
: m_length(ea->getChromosomeLength())
, m_initial_hypercube_size(ea->getHypercubeInitialSize())
, m_fitness(0)
, m_chromosome(m_length)
, m_partial_solution(m_length, filename)
, m_solution(m_length)
{
  UpdateSolution();
};


//------------------------------------------------------------------------------
CIndividual::CIndividual(const CIndividual& ind)
: m_length(ind.getLength())
, m_initial_hypercube_size(ind.m_initial_hypercube_size)
, m_fitness(ind.hasFitnessObject()==true?ind.m_fitness->Duplicate():0)
, m_chromosome(ind.m_chromosome)
, m_partial_solution(ind.m_partial_solution)
, m_solution(m_length)
{
  UpdateSolution();
};


//------------------------------------------------------------------------------
const CIndividual& CIndividual::operator=(const CIndividual& rhs) 
{
  if (&rhs != this) 
  {   
    if (hasFitnessObject()==false)   // copy fitness object
    { 
      if (rhs.hasFitnessObject()==true) 
        setFitnessObject(rhs.getFitnessObject()->Duplicate());
    }
    else if (rhs.hasFitnessObject()==true) (*m_fitness)=(*rhs.m_fitness);
    else SafeDelete(m_fitness);

    // copy "admin" data
    m_initial_hypercube_size = rhs.m_initial_hypercube_size;
    m_length                 = rhs.m_length;

    // copy chromosome/solution data
    m_chromosome             = rhs.m_chromosome;
    m_partial_solution       = rhs.m_partial_solution;
    m_solution               = rhs.m_solution;
  }
  return (*this);
};


//------------------------------------------------------------------------------
double CIndividual::getAllele(const int locus
                             ,const double min_value
                             ,const double max_value
                             ) const 
{
  return getSolution().getAllele(locus);
};



//------------------------------------------------------------------------------
double CIndividual::getMappedAllele(const int locus
                                   ,const double min_value
                                   ,const double max_value
                                   ,const bool constrain
                                   )  
{
  return getSolution().getMappedAllele(locus, min_value, max_value, getInitialHypercubeSize(), constrain);
};


double CIndividual::getSphereMappedAllele(const int locus_1,
                                          const int locus_2,
                                          const double min_value,
                                          const double max_value) const
{
  return getSolution().getSphereMappedAllele(locus_1, locus_2, min_value, max_value,
    getInitialHypercubeSize());
}



//------------------------------------------------------------------------------
double CIndividual::getCircMappedAllele(const int locus
                                       ,const double min_value
                                       ,const double max_value
                                       ) const 
{
  return getSolution().getCircMappedAllele(locus, min_value, max_value, getInitialHypercubeSize());
};


//------------------------------------------------------------------------------
double CIndividual::getCircSphereMappedAllele( const int locus_1
                                 ,const int locus_2
                                 ,const double min_value
                                 ,const double max_value) const
{
  return getSolution().getCircSphereMappedAllele(locus_1, locus_2, min_value,
    max_value, getInitialHypercubeSize());
}


//------------------------------------------------------------------------------
void CIndividual::UpdateSolution(void)
{
  m_solution = m_chromosome;
  m_solution += m_partial_solution;
};


//------------------------------------------------------------------------------
void CIndividual::AddChromosomeToPartialSolution(void) 
{
  m_solution = (m_partial_solution+=m_chromosome);
  m_chromosome.Set(0.0);
};


//------------------------------------------------------------------------------
void CIndividual::OnePointMutation(const double max_mutation_value) 
{
  m_chromosome.OnePointMutation(max_mutation_value);
  UpdateSolution();
};

//------------------------------------------------------------------------------
void CIndividual::GeneralMutation(const double max_mutation_value, const double probability)
{
  m_chromosome.GeneralMutation(max_mutation_value, probability);
  UpdateSolution();
}

//------------------------------------------------------------------------------
ostream& operator<<(ostream& os, const CIndividual& rhs)
{
  return os << rhs.m_length 
         << ", " << rhs.m_fitness << ", "
     << rhs.getInitialHypercubeSize() << ", "<< rhs.m_solution;
}


/* function not used at the moment
void CIndividual::ShiftMutation(const double mutation_sdev)
{
  getChromosome()->ShiftMutation(mutation_sdev);
  UpdateSolution();
};*/
