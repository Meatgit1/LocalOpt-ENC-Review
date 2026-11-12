/*!
\file CInTranceFitness.h
\version 1.7
\brief CInTranceFitness class interface.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 24.06.2005 (AO)

\b Purpose:
   This file provides the interface for the class CInTranceFitness.
   CInTranceFitness implements the multi-objective fitness function of InTrance.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.

\b Copyright:
   (c) 2001-2005 Bernd Dachwald, Andreas Ohndorf
*/

#ifndef CINTRANCEFITNESS_H_INCLUDED
#define CINTRANCEFITNESS_H_INCLUDED

// InTrance
#include <InTranceConst.h>

// ENCLib
#include <CEAFitness.h>
#include <string>
#include <vector>
#include <utility>

/// returns fitness value of actual optimization case where boundary
/// conditions are matched
double GetCriticalFitnessValue(void);

/// Class that implements the multi-objective fitness function of InTrance.
class CInTranceFitness: virtual public CEAFitness
{

public:

  /*!
   * \brief Constructor.
   */
  CInTranceFitness(const int dim = N_SUB_FITNESS) : CEAFitness(dim) {};

  /*!
   * \brief Destructor.
   */
  virtual ~CInTranceFitness() {};

  /*!
   * \brief Returns overall fitness.
   */
  virtual double getMixedFitness(void);

  /*!
   * \brief Returns pointer to a copy of itself.
   */
  virtual CEAFitness* Duplicate(void) const;

}; // end class CEAFitness


/*!
\brief 
*/
struct FitnessRecStruct
{
  std::vector<std::pair< std::string, float> > data;

  /// standard constructor
  FitnessRecStruct(void) : data() {}
};
typedef struct FitnessRecStruct EAFitnessRec;


#endif // CEVFITNESS_H_INCLUDED
