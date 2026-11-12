/*!
\file CEAFitness.h
\brief An auxiliarary class to measure the fitness of an individual.

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

  22.05.2005 (AO)

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

#ifndef CEAFITNESS_H_INCLUDED
#define CEAFITNESS_H_INCLUDED

// STL
#include <sstream>

// SpaceLib
#include <Global.h>
#include <MLVector.h>

// ENCLib
#include <configENC.h>


/*! \brief evolutionary algorithm multi-objective fitness object
 * base class
 *
 */
class CEAFitness 
{

protected:

  MathLib::Vector m_fitness; // vector of fitness function values

public:

   ///brief standard constructor.
  CEAFitness(void) : m_fitness() {};


  /// copy constructor
  CEAFitness(const CEAFitness& to_copy) : m_fitness(to_copy.m_fitness) {};


  /*!
   * \brief size constructor
   * \param dim Number of desired fitness function values.
   */
  CEAFitness(const int dim) : m_fitness(dim)
  {
    m_fitness.Set(MIN_FITNESS);
  };


  /*!
   * \brief initialization constructor
   * \param dim number of desired fitness function values.
   * \param init_val value to set all sub fitnesses with
   */
  CEAFitness(const int dim, const double init_val) : m_fitness(dim) 
  {
    if (init_val!=0) m_fitness.Set(init_val);
  }


  /// destructor
  virtual ~CEAFitness() {};


  /*!
   * \brief number of fitness function values.
   */
  inline int getDimension(void) const 
  {
    return m_fitness.getDimension();
  };


  /*!
   * \brief Returns a specific subfitness function value.
   */
  double getFitness(const int component) const;

  /// reset fitness object
  void reset(void)
  {
    m_fitness.Set(MIN_FITNESS);
  } 

  /// reset fitness object with given value
  void reset(const double init_val)
  {
    m_fitness.Set(init_val);
  }


  /*!
   * \brief Sets a specific subfitness function value.
   */
  void   setFitness(const int component
                   ,const double value);


  /*!
   *  \brief Returns overall fitness value.
   */
  virtual double getMixedFitness(void) = 0;


  /*!
   * \brief Returns pointer to a copy of itself.
   */
  virtual CEAFitness* Duplicate(void) const = 0;

  /*!
   * \brief Assign operator.
   */
  virtual const CEAFitness& operator=(const CEAFitness& rhs)
  {
    if (rhs.getDimension() != getDimension())
    {
      std::ostringstream s("CEAFitness::operator=(const CEAFitness&): unequal \
 dimensions ");
      s << "(lhs: " << getDimension() <<"; rhs: "<< rhs.getDimension()<<")\n";
      Warning(s.str());
    }
    m_fitness.Set(rhs.m_fitness);
    return *this;
  };


  /*!
   * \brief Compare operator.
   */
  virtual inline bool operator==(const CEAFitness& rhs) const
  {
    if (m_fitness == rhs.m_fitness) return true;
    else return false;
  }

}; // end class CEAFitness

#endif // CEVFITNESS_H_INCLUDED
