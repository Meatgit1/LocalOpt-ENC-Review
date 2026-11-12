/*!
\file CSCAttitudeHistory.h
\brief A class that provides the functionality for handling the attitude history of a spacecraft.

\version

  1.7

\author

  Bernd Dachwald\n
  German Aerospace Center (DLR)\n
  Institute of Space Simulation\n
  Linder Hoehe\n
  51147 Cologne\n
  Germany\n
  bernd.dachwald@dlr.de

\date

  26.06.2005(AO)

\b Compiler:

  Microsoft Visual C++ .NET

\b Copyright:

  (c) 2001-2005 by Bernd Dachwald

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/
#ifndef CSPACECRAFTATTITUDEHISTORY_H_INCLUDED
#define CSPACECRAFTATTITUDEHISTORY_H_INCLUDED

#include <string>

#include <CIndividual.h>
#include <CSCAttitude.h>


/// A class that provides the functionality for handling the attitude history of a spacecraft.
class CSCAttitudeHistory
{
  /// number of CSCAttitude-elements
  int           m_n;
  /// CSCAttitude-elements
  CSCAttitude** m_v;

// === Class Interface ===

public:
  /// Construct CSCAttitudeHistory-object from file
  /**
    \param filename pointer to the filename
  */
  CSCAttitudeHistory(const std::string& filename);
  /// Construct CSCAttitudeHistory-object from solution
  /**
    \param solution tbd
    \param mjd0 tbd
    \param dmjd tbd
  */
  CSCAttitudeHistory(const CIndividual& solution, const double mjd0, const double dmjd);

  /// Destruct CSCAttitudeHistory-object
  ~CSCAttitudeHistory();

  /// Get the number of attitude history elements (m_n)
  int getDimension(void) const;

  /// wrapper call for STL look 'n feel
  inline int size(void) {return getDimension();};


  /// Get the modified Julian date (MJD) for an element
  /**
    \param element atttitude history element
  */
  double getMJD(int element) const;
  /// Get the clock angle for an element
  /**
    \param element atttitude history element
  */
  double getClockAngle(int element) const;
  /// Get the cone angle for an element
  /**
    \param element atttitude history element
  */
  double getConeAngle(int element) const;

}; // end class CSCAttitudeHistory

#endif // CSPACECRAFTATTITUDEHISTORY_H_INCLUDED
