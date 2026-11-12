/*!
\file CSolarSail.h
\brief A class that provides the functionality for a solar sail.

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

  29.10.2006 (AO)

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
#ifndef CSOLARSAIL_H_INCLUDED
#define CSOLARSAIL_H_INCLUDED

#include <CPropSys.h>
#include <new>

class CSC;


/// A class that provides the functionality for a solar sail.
class CSolarSail : public CPropSys
{
  friend class CSC;

               /// solar sail(craft) lightness number (lambda)
  double       m_lightness;
               double getLightness(void) const { return m_lightness; };
               void setLightness(const double lambda) { m_lightness = lambda; };

               /// is it an ideal or a real sail
  bool         m_ideal_or_real_sail;
               // true, if ideal sail
               // false, if real sail
               bool isIdealSail(void) const { return m_ideal_or_real_sail; };
               bool isRealSail (void) const { return !m_ideal_or_real_sail; };

               /// solar radiation dose [solar constants * years]
  double       m_srd;
      public:  double getSolarRadiationDose(void) const { return m_srd; };
               void setSolarRadiationDose(const double srd) { m_srd = srd; };
      private:
               
               /// maximum tolerable temerature of the sail film (only used for real sails)
  double       m_max_temperature;
      public:  double getMaxTemperature(void) const { return m_max_temperature; };
      private:

               /// sail's reflection coefficient, when sail is new
  double       m_rho_0;
               /// sail's reflection coefficient, when sail is aged
  double       m_rho_infty;
               /// half life for sail's reflection coefficient
  double       m_rho_half;
               /// degradation constant for sail's reflection coefficient
  double       m_rho_lambda;
      public:  double getrho(void) const { return m_rho_infty + (m_rho_0 - m_rho_infty) * exp(-m_rho_lambda*m_srd); };
      private:
               /// sail's specular reflection coefficient, when sail is new
  double       m_s_0;
               /// sail's specular reflection coefficient, when sail is aged
  double       m_s_infty;
               /// half life for sail's specular reflection coefficient
  double       m_s_half;
               /// degradation constant for sails's specular reflection coefficient
  double       m_s_lambda;
      public:  double gets(void) const { return m_s_infty + (m_s_0 - m_s_infty) * exp(-m_s_lambda*m_srd); };
      private:

               /// sail's front side emission coefficient, when sail is new
  double       m_epsilonf_0;
               /// sail's front side emission coefficient, when sail is aged
  double       m_epsilonf_infty;
               /// half life for sail's front side emission coefficient
  double       m_epsilonf_half;
               /// degradation constant for sail's front side emission coefficient
  double       m_epsilonf_lambda;
      public:  double getepsilonf(void) const { return m_epsilonf_infty + (m_epsilonf_0 - m_epsilonf_infty) * exp(-m_epsilonf_lambda*m_srd); };
      private:

               /// sail's back side emission coefficient, when sail is new
  double       m_epsilonb_0;
               /// sail's back side emission coefficient, when sail is aged
  double       m_epsilonb_infty;
               /// half life for sail's back side emission coefficient
  double       m_epsilonb_half;
               /// degradation constant for sail's back side emission coefficient
  double       m_epsilonb_lambda;
      public:  double getepsilonb(void) const { return m_epsilonb_infty + (m_epsilonb_0 - m_epsilonb_infty) * exp(-m_epsilonb_lambda*m_srd); };
      private:

               /// sail's front side non-Lambertian coefficient
  double       m_Bf;
      public:  double getBf(void) const { return m_Bf; };
      private:

               /// sail's back side non-Lambertian coefficient
  double       m_Bb;
      public:  double getBb(void) const { return m_Bb; };
      private:

  bool CSolarSail::operator==(const CSolarSail& rhs)
  {
    return (
      (m_lightness          == rhs.m_lightness) &&
      (m_ideal_or_real_sail == rhs.m_ideal_or_real_sail) &&
      (m_srd                == rhs.m_srd) &&
      (m_max_temperature    == rhs.m_max_temperature) &&
      (m_rho_0              == rhs.m_rho_0) &&
      (m_s_0                == rhs.m_s_0) &&
      (m_epsilonf_0         == rhs.m_epsilonf_0) &&
      (m_epsilonb_0         == rhs.m_epsilonb_0) &&
      (m_rho_infty          == rhs.m_rho_infty) &&
      (m_s_infty            == rhs.m_s_infty) &&
      (m_epsilonf_infty     == rhs.m_epsilonf_infty) &&
      (m_epsilonb_infty     == rhs.m_epsilonb_infty) &&
      (m_rho_half           == rhs.m_rho_half) &&
      (m_s_half             == rhs.m_s_half) &&
      (m_epsilonf_half      == rhs.m_epsilonf_half) &&
      (m_epsilonb_half      == rhs.m_epsilonb_half) &&
      (m_rho_lambda         == rhs.m_rho_lambda) &&
      (m_s_lambda           == rhs.m_s_lambda) &&
      (m_epsilonf_lambda    == rhs.m_epsilonf_lambda) &&
      (m_epsilonb_lambda    == rhs.m_epsilonb_lambda) &&
      (m_Bf                 == rhs.m_Bf) &&
      (m_Bb                 == rhs.m_Bb) &&
      (CPropSys::operator==(rhs)));
  };


public:

  /// Constructor with given host state, i.e. a CSC object 
  /// (which is also a state object)
  CSolarSail(CState* state = 0) : CPropSys(SOLARSAIL, state)
  {
    m_lightness          = 0.0;
    m_ideal_or_real_sail = true;
    m_srd                = 0.0;
    m_max_temperature    = HUGE_VAL;
    m_rho_0              = 0.0;
    m_s_0                = 0.0;
    m_epsilonf_0         = 0.0;
    m_epsilonb_0         = 0.0;
    m_rho_infty          = 0.0;
    m_s_infty            = 0.0;
    m_epsilonf_infty     = 0.0;
    m_epsilonb_infty     = 0.0;
    m_rho_half           = 0.0;
    m_s_half             = 0.0;
    m_epsilonf_half      = 0.0;
    m_epsilonb_half      = 0.0;
    m_rho_lambda         = 0.0;
    m_s_lambda           = 0.0;
    m_epsilonf_lambda    = 0.0;
    m_epsilonb_lambda    = 0.0;
    m_Bf                 = 0.0;
    m_Bb                 = 0.0;
  }

  /// copy constructor
  CSolarSail(const CSolarSail& sail, CState* state = 0) : CPropSys(sail, state) 
  {
    m_lightness          = sail.m_lightness;
    m_ideal_or_real_sail = sail.m_ideal_or_real_sail;
    m_srd                = sail.m_srd;
    m_max_temperature    = sail.m_max_temperature;
    m_rho_0              = sail.m_rho_0;
    m_s_0                = sail.m_s_0;
    m_epsilonf_0         = sail.m_epsilonf_0;
    m_epsilonb_0         = sail.m_epsilonb_0;
    m_rho_infty          = sail.m_rho_infty;
    m_s_infty            = sail.m_s_infty;
    m_epsilonf_infty     = sail.m_epsilonf_infty;
    m_epsilonb_infty     = sail.m_epsilonb_infty;
    m_rho_half           = sail.m_rho_half;
    m_s_half             = sail.m_s_half;
    m_epsilonf_half      = sail.m_epsilonf_half;
    m_epsilonb_half      = sail.m_epsilonb_half;
    m_rho_lambda         = sail.m_rho_lambda;
    m_s_lambda           = sail.m_s_lambda;
    m_epsilonf_lambda    = sail.m_epsilonf_lambda;
    m_epsilonb_lambda    = sail.m_epsilonb_lambda;
    m_Bf                 = sail.m_Bf;
    m_Bb                 = sail.m_Bb;
  };

  /// copy operator.
  CSolarSail& operator=(const CSolarSail& sail)
  {
    CPropSys::operator=(sail);
    m_lightness          = sail.m_lightness;
    m_ideal_or_real_sail = sail.m_ideal_or_real_sail;
    m_srd                = sail.m_srd;
    m_max_temperature    = sail.m_max_temperature;
    m_rho_0              = sail.m_rho_0;
    m_s_0                = sail.m_s_0;
    m_epsilonf_0         = sail.m_epsilonf_0;
    m_epsilonb_0         = sail.m_epsilonb_0;
    m_rho_infty          = sail.m_rho_infty;
    m_s_infty            = sail.m_s_infty;
    m_epsilonf_infty     = sail.m_epsilonf_infty;
    m_epsilonb_infty     = sail.m_epsilonb_infty;
    m_rho_half           = sail.m_rho_half;
    m_s_half             = sail.m_s_half;
    m_epsilonf_half      = sail.m_epsilonf_half;
    m_epsilonb_half      = sail.m_epsilonb_half;
    m_rho_lambda         = sail.m_rho_lambda;
    m_s_lambda           = sail.m_s_lambda;
    m_epsilonf_lambda    = sail.m_epsilonf_lambda;
    m_epsilonb_lambda    = sail.m_epsilonb_lambda;
    m_Bf                 = sail.m_Bf;
    m_Bb                 = sail.m_Bb;
    return *this;
  };


  /*!
   \brief make copy of object and set given host object pointer

    Because polymorphism is used in implementation of different propulsion
    systems a clone() - method is needed. Due to the fact that only the object
    itself knows its real type (except you try a dynamic_cast<..> with all
    possible object types and handle possibly thrown exceptions) it makes sense
    to have such a clone function which returns a base class pointer to the
    new object which is a copy of itself.

   \param state pointer to "host" object
   \return CPropSys Pointer to new CSolarSail object
  */
  virtual CPropSys* clone(CState* state)
  {
    try 
    {
      return (new CSolarSail(*this, state));
    } 
    catch (std::bad_alloc ) 
    {
      return 0;
    }
  };


  /// destructor
  virtual ~CSolarSail(void)  {}

  // __mutators__

  void setOpticalModel(const bool ideal_or_real,
                       const double srd = 0.0,
                       const double max_temperature = HUGE_VAL,
                       const double rho_0 = 0.0,
                       const double s_0 = 0.0,
                       const double epsilonf_0 = 0.0,
                       const double epsilonb_0 = 0.0,
                       const double rho_infty = 0.0,
                       const double s_infty = 0.0,
                       const double epsilonf_infty = 0.0,
                       const double epsilonb_infty = 0.0,
                       const double rho_half = 0.0,
                       const double s_half = 0.0,
                       const double epsilonf_half = 0.0,
                       const double epsilonb_half = 0.0,
                       const double Bf = 0.0,
                       const double Bb = 0.0);

  // __functionality__

  /// Returns the current sail film temperature.
  double getTemperature(const double r, const double coneangle) const;

  /// Returns the current solar energy flux on the sail film (in times of the solar constant).
  double getSolarEnergyFlux(const double r, const double coneangle) const;

  void Init(const double srd) { m_srd = srd; };
  
  double getG(void) const
  {
    return 1.0 + gets() * getrho();
  };

  double getH(void) const
  {
    return 1.0 - gets() * getrho();
  };

  double getK(void) const
  {
    return getBf()*(1.0-gets())*getrho()+(1.0-getrho())*(getepsilonf()*getBf()-getepsilonb()*getBb())/(getepsilonf()+getepsilonb());
  };

  double getG0(void) const
  {
    return 1.0 + m_s_0 * m_rho_0;
  };

  double getH0(void) const
  {
    return 1.0 - m_s_0 * m_rho_0;
  };

  double getK0(void) const
  {
    return m_Bf*(1.0-m_s_0)*m_rho_0+(1.0-m_rho_0)*(m_epsilonf_0*m_Bf-m_epsilonb_0*m_Bb)/(m_epsilonf_0+m_epsilonb_0);
  };

public:

  //------------------------------------------------------------------------------
  //
  //   equations of motion for solar sailcraft in RK-form
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
  //               mu             mu * cos^3(beta)
  //             - --- + lambda * ----------------
  //               r^2                  r^2
  //
  //
  //                 dr * dphi
  //     ddphi = - 2 ---------
  //                     r
  //
  //             + 2 * dphi * dtheta * tan(theta)
  //
  //                        mu * cos(alpha+zeta) * sin(beta) * cos^2(beta)
  //             + lambda * ----------------------------------------------
  //                                       r^3 * cos(theta)
  //
  //
  //                 dr * dtheta
  //   ddtheta = - 2 -----------
  //                      r
  //
  //             - dphi^2 * sin(theta) * cos(theta)
  //
  //                        mu * sin(alpha+zeta) * sin(beta) * cos^2(beta)
  //             + lambda * ----------------------------------------------
  //                                             r^3
  //
  //
  //         r = y(0)
  //       phi = y(1)
  //     theta = y(2)
  //        dr = y(3)
  //      dphi = y(4)
  //    dtheta = y(5)
  //       ddr = s.a.
  //     ddphi = s.a.
  //   ddtheta = s.a.
  //
  //------------------------------------------------------------------------------
  void EOM_RK_POLAR(double t, const MathLib::Vector& y, MathLib::Vector & yp, void * pAux);

  void EOM_RK_CARTESIAN(double t, const MathLib::Vector& y, MathLib::Vector & yp, void * pAux);

  void CheckRemainingPropellant(const double dt) {};

  /*! \brief Returns zero, because a solar sail doesn't have an engine.*/
  double getThrust(void) {return 0.0;};

  virtual void printVrmlAppearance(std::ostream& os, const std::string& name)
    const throw ();

}; // end class CSolarSail


#endif // CSOLARSAIL_H_INCLUDED