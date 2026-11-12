/*!
\file CChromosome.h
\brief A class that provides the functionality for individual chromosomes 
of evolutionary algorithm.

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

  29.08.2006 (AO)

\b Compiler:

  GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++

\b Copyright:

  (c) 2001-2006 by Bernd Dachwald, Andreas Ohndorf (AO)

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/
#ifndef CCHROMOSOME_H_INCLUDED
#define CCHROMOSOME_H_INCLUDED

// SpaceLib
#include <MLVector.h>
#include <Global.h>

// ENCLib
#include <configENC.h>


/// chromosome class (used in CIndividual)
class CChromosome : public MathLib::Vector
{

  /*!\brief set random allele values within given scope to all chromosome loci
   * Assigns random values within the interval [<min_value>, <max_value>]
   * to all alleles of the chromosome.
   * \param min_value minimum allele value to set
   * \param max_value maximum allele value to set
   */
  void setChromosomeRandom(const double min_value, const double max_value);

  /* function not used at the moment
  double Norm(const VectorNorm norm);
              // type of vector norm

  // purpose: calculate the norm of the chromosome
  // return : norm
  // note   : the available vector norms are declared in CEATypes.h*/

  void Load(const char* filename);

  // purpose: load chromosome from file <filename>
  // note   : the file must have the following structure:
  //            1st line: number of values
  //            2nd line: value[0]
  //            nth line: value[n-2]

//=============================================================================
// CChromosome CLASS INTERFACE
//=============================================================================
public:

  /// standard constructor
  CChromosome(void) : MathLib::Vector() {};

  /// copy constructor (reference)
  CChromosome(const CChromosome& chrom) : MathLib::Vector(chrom) {};

  /*! \brief length, lower and upper allele value boundary constructor
   *
   *  The created chromosome's allele values are randomly in the intervall
   *  [<min_value>,<max_value>]. If no minimum and maximum values are provided
   *  a default value of zero will be used.
   *
   *  \param length Length of chromosom; number of real values.
   *  \param min_value Lower boundary for random numbers (default: 0).
   *  \param min_value Upper boundary for random numbers (default: 0).
   */
  CChromosome(const int length,
              const double min_value = 0.0,
              const double max_value = 0.0)
  : MathLib::Vector(length)
  {
    setChromosomeRandom(min_value, max_value);
  };

  
  /*! \brief data file constructor (const char*)
   *
   * Creates a chromosome from file with given filename. The structure
   * of given file should obey the following notation:\n
   * 1st line: number of values\n
   * 2nd line: value[0]\n
   * nth line: value[n-2]\n
   * \param length number of chromosome values (chromosome length)
   * \param filename name of chromosome data file
   * 
   */
  CChromosome(const int length, const char* filename) : MathLib::Vector(length)
  { 
    Load(filename);
  };


  ///data file constructor (const std::string)
  CChromosome(const int length, const std::string& filename)
    : MathLib::Vector(length)
  {
    Load(filename.c_str());
  }


  /// add/assign operator
  CChromosome& operator+=(const CChromosome& chrom)  
  {
    Vector::operator +=(chrom);
    return (*this);
  };


  /// destructor
  ~CChromosome(void) {};


  /*! \brief allele value at locus <locus>
   *
   * \note range of <locus> is checked [#ifdef development_version]
   * \return value of allele at locus <locus>
   */
  double getAllele(const int locus) const;


  /*! \brief set value <value> of allele at locus <locus>
   *
   * \note range of <locus> is checked [#ifdef development_version]
   */
  void setAllele(const int locus, const double value);

  /*! \brief maps allele value onto given interval
   *
   * Allele value at desired locus will be set into relation
   * to original range (orig_range) projected onto interval
   * determined by given lower and upper boundary
   *
   * \param locus allele to access
   * \param min_value lower boundary of mapping scope
   * \param max_value upper boundary of mapping scope
   * \param orig_range former max allele value
   * \return mapped allele value
   */
  double getMappedAllele(const int locus,
                         const double min_value,
                         const double max_value,
                         const double orig_range,
                         const bool constrain = false);

  double getSphereMappedAllele(const int locus_1,
                               const int locus_2,
                               const double min_value,
                               const double max_value,
                               const double orig_range) const;

  /*! \brief maps allele value onto given angle interval
   *
   * Allele value at desired locus will be set into relation
   * to original range (orig_range) projected onto interval
   * determined by given lower and upper boundary angles.
   * The mapping interval may be thought of a unit circle
   * and the allele value is mapped on the circumference.
   *
   * \param locus allele to access
   * \param min_value lower angle boundary of mapping scope
   * \param max_value upper angle boundary of mapping scope
   * \param orig_range former max allele value
   * \return mapped allele value
   */
  double getCircMappedAllele(const int locus,
                             const double min_value,
                             const double max_value,
                             const double orig_range) const;

  double getCircSphereMappedAllele(const int locus_1
                                  ,const int locus_2
                                  ,const double min_value
                                  ,const double max_value
                                  ,const double orig_range) const;



  /*! \brief mutate a single allele value
   *
   * Select random chromosome locus and assign a value in the
   * interval [-<max_mutation_value>,<max_mutation_value>] to that allele.
   * \param max_mutation_value mutation value intervall boundary
   */
  void OnePointMutation(const double max_mutation_value);

  /*! \brief mutate all chromosome alleles with given probability
   *
   * Check for all alleles of a chromosome whether it should be mutated or
   * not (parameter probability) and set the new allele value from the
   * interval [-<max_mutation_value>,<max_mutation_value>].
   * \param max_mutation_value mutation value interval boundary
   * \param probability a value in (0.0,1.0) where 0.0 means no mutation and 1.0
   * always mutate
   */
  void GeneralMutation(const double max_mutation_value, const double probability);


  /* function not used at the moment
  void ShiftMutation(const double mutation_sdev);


  /*! \brief save chromosome content to file (const char*)
   *
   * note   : the file has the following structure:\n
   *            1st line: number of values\n
   *            2nd line: value[0]\n
   *            nth line: value[n-2]   
   * \param filename data file name to save chromosome data into
   */
  void Save(const char* filename) const;


  /// save chromosome content to file (const std::string&)
  void Save(const std::string& filename) const;

}; // end of CChromosome class declaration

#endif // CCHROMOSOME_H_INCLUDED
