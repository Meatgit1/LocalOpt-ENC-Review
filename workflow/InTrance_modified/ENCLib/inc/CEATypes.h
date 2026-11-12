/*!
\file CEATypes.h
\brief A header that provides type declarations for the evolutionary algorithm.

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

  28.12.2005 (AO)

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
#ifndef CEATYPES_H_INCLUDED
#define CEATYPES_H_INCLUDED


// available vector norms
enum VectorNorm {
  MAXIMUM_SUM_NORM     = 0, // (1-norm)
  EUCLIDEAN_NORM       = 1, // (2-norm)
  MAXIMUM_ELEMENT_NORM = 2  // (infinity norm)
};


/*!
\brief Evolutionary algorithm status information record structure type
*/
struct EAStatusRecStruct
{
  unsigned int population_size;
  float mutation_probability;
  float gen_mutation_probability;
  unsigned int epoch;
  unsigned int bad_epochs;
  unsigned int sss_epoch;
  unsigned int n_sss_epochs;
  float hc_size;
  float hc_upper_limit;
  float hc_sf;
  float hc_init_size;
  float fbc_met_limit;
  float pop_conv_limit1;
  float pop_conv_limit2;
  float ea_conv_limit1;
  float ea_conv_limit2;
  float best_fitness_epoch;
  float best_fitness_ever;
  float current_best_fitness;
  float improvement;
  float avg_fitness;
  float reproduction_rate;
  unsigned long int reproduction;
  unsigned long int epoch_reproduction;

  /// standard constructor
  EAStatusRecStruct(void)
    : population_size(0)
    , mutation_probability(0.0)
    , gen_mutation_probability(0.0)
    , epoch(0)
    , bad_epochs(0)
    , sss_epoch(0)
    , n_sss_epochs(0)
    , hc_size(0.0f)
    , hc_upper_limit(0.0f)
    , hc_sf(0.0f)
    , hc_init_size(0.0f)
    , fbc_met_limit(0.0f)
    , pop_conv_limit1(0.0f)
    , pop_conv_limit2(0.0f)
    , ea_conv_limit1(0.0f)
    , ea_conv_limit2(0.0f)
    , best_fitness_epoch(0.0f)
    , best_fitness_ever(0.0f)
    , current_best_fitness(0.0f)
    , improvement(0.0f)
    , avg_fitness(0.0f)
    , reproduction_rate(0.0)
    , reproduction(0)
    , epoch_reproduction(0)
  {}
};
typedef struct EAStatusRecStruct EAStatusRec;

#endif // CEATYPES_H_INCLUDED
