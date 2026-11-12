/*!
\file CLagrange.h
\version 1.7
\brief A class derived from class CBody for Lagrange equilibrium points (CLagrange).
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 09.08.2004 (AO)

\b Purpose:
   The class CLagrange makes it possible to treat a Lagrangian equlibrium point like
   any other body in solar system (planets, asteroids, comets). It uses some
   functionality of the CSolarSystem class to calculate the position and velocity
   data of the Lagrange point and sets its own data members accordingly. So you
   can get the same information as from any other object stored within a CSolarSystem
   object. The name of a Lagrange point is rather simple and straightforward. It
   starts with the capital letter "L" followed by the number of the Lagrange point,
   then an underscore, the name of Body 1 (usually the sun), followed by another
   underscore, and finally the name of Body 2 (usually a planet). As an example,
   the name of L4 in the Sun-Jupiter-system would be: "L4_sun_jupiter". The name may
   change, if the getName member function of the CBody class changes.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.
   
\b Copyright:
   (c) 2004 Andreas Ohndorf
*/

#ifndef CLAGRANGE_H
#define CLAGRANGE_H

#include <CBody.h>
#include <CSolarSystem.h>

namespace SpaceLib
{
  /** \class CLagrange
      The class CLagrange makes it possible to treat a Lagrangian equlibrium point like
      any other body in solar system (planets, asteroids, comets). It uses some
      functionality of the CSolarSystem class to calculate the position and velocity
      data of the Lagrange point and sets its own data members accordingly. So you
      can get the same information as from any other object stored within a CSolarSystem
      object. The name of a Lagrange point is rather simple and straightforward. It
      starts with the capital letter "L" followed by the number of the Lagrange point,
      then an underscore, the name of Body 1 (usually the sun), followed by another
      underscore, and finally the name of Body 2 (usually a planet). As an example,
      the name of L4 in the Sun-Jupiter-system would be: "L4_sun_jupiter". The name may
      change, if the getName member function of the CBody class changes.*/
  class CLagrange : virtual public CBody
  {
  private:

    // The following two BodyIDs define the equilibrium point system because
    // the five Lagrange points are defined in this system. In general Body 1 is
    // the Sun or Earth and Body 2 is the planet or moon. The computation of
    // Lagrange equilibrium points makes only sense if Body 1 is the central
    // body for Body 2.

    BodyID L_Body_1; ///< ID of body 1 (this is the central body, the one with the greater mass)
    BodyID L_Body_2; ///< ID of body 2 (the body with the smaller mass)
    int    m_LPoint; ///< number of the Lagrange point (from 1 to 5)

    /** Tbd. */
    virtual void SetState(const double mjd);

   public:

   /** \brief Constructor for Lagrange equlibrium point.
       Initializes the CLagrange object. In order to create a meaningful
       object, a few checks are performed. The arguments b1 and b2 may be
       obtained by calling the getBody member function of the CSolarSystem
       class with the appropriate BodyIDs.
       \param b1 pointer to CBody object 1 (center body, e.g. sun)
       \param b2 pointer to CBody object 2 (e.g. jupiter)
       \param which number of the Lagrange point (from 1 to 5) */
    explicit CLagrange(CBody * b1, CBody * b2, const int which);

    /** Destructor. */
    virtual ~CLagrange(void);
  
  }; // end class CLagrange

} // end of namespace SpaceLib

#endif // CLAGRANGE_H
