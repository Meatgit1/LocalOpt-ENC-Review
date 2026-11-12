/*!
\file CSCAttitude.cpp
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

  27.08.2006 (AO)

\b Compiler:

  Microsoft Visual C++ .NET

\b Copyright:

  (c) 2001-2006 by Bernd Dachwald

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/

// STL
#include <sstream>

//SCLib
#include <configSC.h>
#include <SpaceConst.h>
#include <CSC.h>

using namespace std;

// ----------------------------------------------------------------------------
void CSCAttitude::setClockAngle(const double alpha)
// ----------------------------------------------------------------------------
// private
{
  //bool warning = false;

  double corrected_alpha = alpha;

  // if necessary, make corrections and produce warning

  if (alpha < -PI)
  {
    corrected_alpha = -PI;
    //warning = true;
  }
  else if (alpha >  PI)
  {
    corrected_alpha = PI;
    //warning = true;
  }
  //if (warning)
  //{
  //  ostringstream s;
  //  s << "CSCAttitude::setClockAngle(): clock angle had to be corrected (";
  //  s << alpha << ")";
  //  Warning(s.str());
  //}

  v[CLOCK] = corrected_alpha;

}; // end CSCAttitude::setClockAngle


// ----------------------------------------------------------------------------
void CSCAttitude::setConeAngle(const double beta)
// ----------------------------------------------------------------------------
// private
{
  bool warning = false;
  double corrected_beta = beta;

  // if necessary, make corrections and produce warning

  if (beta < 0.0)
  {
    corrected_beta = 0.0;
    //warning = true;
  }
  else if (beta > PI)
  {
    corrected_beta = PI;
    warning = true;
  }
  if (m_sc) // not a standalone object
  {
    if (m_sc->isSolarSail()) // special restrictions for solar sail
    {
      if (beta > HALFPI)
      {
        corrected_beta = HALFPI;
        warning = true;
      }
    }
  }
  if (warning)
  {
    ostringstream s;
    s << "CSCAttitude::setConeAngle(): cone angle had to be corrected (";
    s << beta << ">" <<HALFPI<<")";
    Warning(s.str());
  }

  v[CONE] = corrected_beta;

}; // end CSCAttitude::setConeAngle


// ----------------------------------------------------------------------------
void CSCAttitude::setAttitude(const double alpha,
                              const double beta)
// ----------------------------------------------------------------------------
// private
{
  setClockAngle(alpha);

  double minbeta = 0.0;
  if (m_sc->isSolarSail() && m_sc->isRealSail())
  {
    const double term1 = (1.0-m_sc->getrho())/(m_sc->getepsilonf()+m_sc->getepsilonb());
    const double r = m_sc->CState::r();
    const double term2 = S0_AUD/(SIGMA_AUD*r*r);
    const double cosminbeta = std::min(pow(m_sc->getMaxTemperature(), 4.0) / (term1*term2) , 1.0);
    minbeta = acos(cosminbeta);
  }
  setConeAngle(std::max(beta,minbeta));
  setConeAngle(beta);
};

// ----------------------------------------------------------------------------
void CSCAttitude::setAttitude(const double xi_r,
                              const double xi_t,
                              const double xi_h)
// ----------------------------------------------------------------------------
// private
{
  const std::string method_name = "CSCAttitude::setAttitude(const double\
                                  , const double, const double):";
  if (m_sc == 0) 
    FatalError(method_name+ std::string("execution not possible for \
                                        standalone objects"));

  // calculate clock angle
  const double clock = ((xi_h==0.0) && (xi_t==0.0) ? 0.0 :atan2(xi_h,xi_t));

  // calculate cone angle
  double cone = 0.0;

  if (m_sc->isSolarSail())
  {
    // calculate required cone angle
    double tanconereq = HUGE_VAL;
    if(xi_r != 0.0) tanconereq = sqrt(xi_t*xi_t+xi_h*xi_h)/xi_r;

    // calculate real cone angle required for maximization in this direction
    // [1] B.Dachwald, "Interplanetary Mission Analysis for Non-Perfectly Reflecting Solar
    // Sailcraft Using Evolutionary Neurocontrol", AAS 03-579, 2003
    if (tanconereq != 0.0)
    {
      if (xi_r > 0.0)
        cone = atan((-3.0+sqrt(9.0+8.0*tanconereq*tanconereq))/(4.0*tanconereq));
      else if (xi_r < 0.0)
        cone = atan((-3.0-sqrt(9.0+8.0*tanconereq*tanconereq))/(4.0*tanconereq));
      else // if (xi_r == 0.0)
        cone = atan(1.0/sqrt(2.0));
    }
    else //(tanconereq==0.0)
    {
      if (xi_r > 0.0) cone = 0.0;
      else if (xi_r < 0.0) cone = HALFPI;
      else cone = atan(1.0/sqrt(2.0));  // if (xi_r == 0.0)
    }

    double minbeta = 0.0;
    if (m_sc->isRealSail())
    {
      // for a real solar sail there is a cone angle's temperature limitation 
      // due to sail exposure to sun's radiation at lower distances from sun
      const double term1 = (1.0-m_sc->getrho())/(m_sc->getepsilonf()+m_sc->getepsilonb());
      const double r = m_sc->CState::r();
      const double term2 = S0_AUD/(SIGMA_AUD*r*r);
      const double cosminbeta = 
          std::min(pow(m_sc->getMaxTemperature(), 4.0) / (term1*term2) , 1.0);
      minbeta = acos(cosminbeta);
      cone = std::max(cone,minbeta);
    }
  }
  else // not a solar sail
  {
    cone = atan2(sqrt(xi_t*xi_t+xi_h*xi_h),xi_r);
  }

  setClockAngle(clock);
  setConeAngle(cone);

}; // end CSCAttitude::setAttitude


// ----------------------------------------------------------------------------
CSCAttitude::CSCAttitude(CSC* sc /* = NULL */)
// ----------------------------------------------------------------------------
// public
: Vector(N_COMPONENTS), m_sc(sc)
{
};

// ----------------------------------------------------------------------------
CSCAttitude::CSCAttitude(const double alpha,
                         const double beta,
                         CSC* sc /* = NULL */)
// ----------------------------------------------------------------------------
// public
: Vector(N_COMPONENTS), m_sc(sc)
{
  setClockAngle(alpha); // use function to assure the correctness of the value
  setConeAngle(beta);   // - " -
}; // end CSCAttitude::CSCAttitude


// ----------------------------------------------------------------------------
void CSCAttitude::setMJD(const double mjd)
// ----------------------------------------------------------------------------
// public
{
  v[TIME] = mjd;
}; // end CSCAttitude::setMJD


// ----------------------------------------------------------------------------
PropulsionSystemType CSCAttitude::getPropulsionSystemType(void) const
// ----------------------------------------------------------------------------
// public
{
  if (!m_sc) FatalError("CSCAttitude::getPropulsionSystemType(): m_sc == NULL");
  return m_sc->getPropulsionSystemType();

}; // end CSCAttitude::getPropulsionSystemType

