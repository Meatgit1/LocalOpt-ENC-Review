/*!
\file CNEP.h
\brief A class that provides the functionality for a nuclear electric propulsion (NEP) system.

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

\author

  Andreas Ohndorf\n
  German Aerospace Center (DLR)\n
  Institute of Space Simulation\n
  Linder Hoehe\n
  51147 Cologne\n
  Germany\n
  andreas.ohndorf@gmx.de

\date

  12.01.2006 (AO)

\b Compiler:

  Microsoft Visual C++ .NET

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
#ifndef CNEP_H_INCLUDED
#define CNEP_H_INCLUDED

#include <CEP.h>

class CSC;

/// A class that provides the functionality for a nuclear electric propulsion (NEP) system.
class CNEP : public CEP
{
  friend class CSC;

               /// maximum thrust of the NEP sytem [N]
  double       m_thrust_max;
      public:  double getMaxThrust(void) { return m_thrust_max; };
      public:  void setMaxThrust(const double thrust_max) { m_thrust_max = thrust_max; };
      private:

               /// specific impulse [s]
  double       m_isp;
      public:  double getIsp(void) { return m_isp; };
      public:  void setIsp(const double isp) { m_isp = isp; };
      private:

public:

  /// Constructor with given host state, i.e. a CSC object 
  /// (which is also a state object)
  CNEP(CState* state = 0) : CEP(NEP, state)
  {
    m_thrust_max = 0.0;
    m_isp        = 0.0;
  };

  /// Copy constructor
  CNEP(const CNEP& nep, CState* state = 0) : CEP(nep, state)
  {
    m_thrust_max = nep.m_thrust_max;
    m_isp        = nep.m_isp;
  };


  /// assign operator
  CNEP& operator=(const CNEP& nep)
  {
    CEP::operator =(nep);
    m_thrust_max = nep.m_thrust_max;
    m_isp        = nep.m_isp;
    return *this;
  };

public:  
  //double getSolarRadiationDose(void) const { return 0.0; };
  //void   setSolarRadiationDose(const double srd) {};

  /*!
   \brief make copy of object and set given host object pointer

    Because polymorphism is used in implementation of different propulsion
    systems a clone() - method is needed. Due to the fact that only the object
    itself knows its real type (except you try a dynamic_cast<..> with all
    possible object types and handle possibly thrown exceptions) it makes sense
    to have such a clone function which returns a base class pointer to the
    new object which is a copy of itself.

   \param state pointer to "host" object
   \return CPropSys Pointer to new CNEP object
  */
  virtual CPropSys* clone(CState* state)
  {
    try 
    {
      return (new CNEP(*this, state));
    } 
    catch (std::bad_alloc ) 
    {
      return 0;
    }
  };


  ~CNEP(void) {};

  void CheckRemainingPropellant(const double dt);
  double getThrust(void);
  double getMassFlow(void);

public:

  //------------------------------------------------------------------------------
  //
  //   equations of motion for nuclear electric propulsion system in RK-form
  //
  // note:
  //
  //   pAux is expected to point to the CSC object
  //
  //
  //       ddr = + r * dtheta^2
  //
  //             + r * dphi^2 * cos^2(theta)
  //
  //               mu    thrust
  //             - --- + ------ * cos(beta)
  //               r^2    mass
  //
  //
  //                 dr * dphi
  //     ddphi = - 2 ---------
  //                     r
  //
  //             + 2 * dphi * dtheta * tan(theta)
  //
  //               thrust   cos(alpha+zeta) * sin(beta)
  //             + ------ * ---------------------------
  //                mass          r * cos(theta)
  //
  //
  //                 dr * dtheta
  //   ddtheta = - 2 -----------
  //                      r
  //
  //             - dphi^2 * sin(theta) * cos(theta)
  //
  //               thrust   sin(alpha+zeta) * sin(beta)
  //             + ------ * ---------------------------
  //                mass                r
  //
  //
  //       dmp = - mflow
  //
  //
  // where:
  //
  //      mass = msc + mp
  //
  //    thrust = throttle * thrust_max
  //
  //              thrust
  //     mflow = --------
  //             g0 * Isp
  //
  //
  //         r = y(0)
  //       phi = y(1)
  //     theta = y(2)
  //        dr = y(3)
  //      dphi = y(4)
  //    dtheta = y(5)
  //        mp = y(6)
  //       ddr = s.a.
  //     ddphi = s.a.
  //   ddtheta = s.a.
  //
  //------------------------------------------------------------------------------
  void EOM_RK_POLAR(double t, const Vector& y, Vector& yp, void* pAux);

  //----------------------------------------------------------------------------
  void EOM_RK_CARTESIAN(double t, const Vector& y, Vector& yp, void* pAux);

}; // end class CNEP



#endif // CNEP_H_INCLUDED
