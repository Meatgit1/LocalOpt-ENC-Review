/*!
\file CSCAttitude.h
\brief A class that provides the functionality for handling the attitude of a spacecraft.

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

  27.08.2006

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
#ifndef CSPACECRAFTATTITUDE_H_INCLUDED
#define CSPACECRAFTATTITUDE_H_INCLUDED


#include <MLVector.h>
#include <SCConst.h>
#include <assert.h>


using namespace MathLib;


class CSC;


/// A class that provides the functionality for handling the attitude of a spacecraft.
class CSCAttitude : public Vector
// ----------------------------------------------------------------------------
// first  vector component: spacecraft time (modified julian date)
// second vector component: spacecraft clock angle (alpha)
// third  vector component: spacecraft cone  angle (beta)
// fourth vector component: x-component of sun distance vector [AU]
// fifth  vector component: y-component of sun distance vector [AU]
// sixth  vector component: z-component of sun distance vector [AU]
// ----------------------------------------------------------------------------
{
  friend class CSC;

  enum COMPONENTS
  {
     TIME
    ,CLOCK
    ,CONE
    ,SUN_R
    ,SUN_T
    ,SUN_H
    ,N_COMPONENTS
  };

  CSC *m_sc; // associated spacecraft (or NULL, if it is a standalone object)


  /*! \brief set steering clock angle
   * \param alpha clock angle*/
  void setClockAngle(const double alpha);

  /*! \brief set steering cone angle
   * \param beta cone angle*/
  void setConeAngle(const double beta);

  /*! \brief sets distance vector from sun
   *  \param l the distance vector from sn to s/c
   *  \param pos index of first value to use from given vector*/
  void setSunVector(const MathLib::Vector& l, const int pos=0);

  /*! \brief set clock and cone angle
   * \param alpha clock angle
   * \param beta cone angle*/
  void setAttitude(const double alpha, const double beta);

  /*! \brief set clock and cone angle according to the optimal thrust vector
   * \param xi_r radial component of the optimal thrust vector (orbit frame)
   * \param xi_t transversal component of the optimal thrust vector 
   (orbit frame)
   * \param xi_h orbit normal component of the optimal thrust vector*/
  void setAttitude(const double xi_r, const double xi_t, const double xi_h);



public:

  /*! \brief constructor to create nullvector attitude
   * \param sc associated spacecraft*/
  CSCAttitude(CSC* sc = NULL);

  /*! \brief constructor to create from clock and cone angle
   *  \param alpha clock angle
   *  \param beta cone angle
   *  \param sc associated spacecraft*/
  CSCAttitude(const double alpha, const double beta, CSC* sc = NULL);

  // set time (modified Julian date)
  void setMJD(const double mjd);

  /// get type of propulsion system for associated spacecraft
  PropulsionSystemType getPropulsionSystemType(void) const;

  /// get time (modified julian date)
  double getMJD(void) const;

  /// get clock angle
  double getClockAngle(void) const;

  /// get cone angle
  double getConeAngle(void) const;

  /*! \brief stored distance vector from sun
   *  \param l vector data structure to store sun distance vector in
   *  \param pos index of first value to use for storage*/
  void getSunVector(MathLib::Vector& l, const int pos=0);

}; // end class CSCAttitude

//------------------------------------------------------------------------------
inline void CSCAttitude::setSunVector(const MathLib::Vector& l, const int pos)
{
#ifdef development_version
  assert((l.size()-pos)>=3);
#endif
  v[SUN_R] = l[pos];
  v[SUN_T] = l[pos+1];
  v[SUN_H] = l[pos+2];
} //  void setSunVector(const MathLib::Vector& l)

//------------------------------------------------------------------------------
inline void CSCAttitude::getSunVector(MathLib::Vector& l, const int pos)
{
#ifdef development_version
  assert((l.size()-pos)>=3);
#endif
  l[pos]   = v[SUN_R];
  l[pos+1] = v[SUN_T];
  l[pos+2] = v[SUN_H];
} // void getSunVector(MathLib::Vector& l)

// ----------------------------------------------------------------------------
inline double CSCAttitude::getMJD(void) const
{
  return v[TIME];
}; // double CSCAttitude::getMJD


// ----------------------------------------------------------------------------
inline double CSCAttitude::getClockAngle(void) const
{
  return v[CLOCK];
}; // double CSCAttitude::getClockAngle


// ----------------------------------------------------------------------------
inline double CSCAttitude::getConeAngle(void) const
{
  return v[CONE];
}; // double CSCAttitude::getConeAngle

#endif // CSPACECRAFTATTITUDE_H_INCLUDED
