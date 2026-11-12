/*!
\file CPropSys.cpp
\brief A class that provides the functionality for a propulsion system.

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

  26.08.2005 (AO)

\b Compiler:

  Microsoft Visual C++ .NET

\b Copyright:

  (c) 2001-2004 by Bernd Dachwald

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/

// SpaceLib
#include <SpaceLib.h>

//SCLib
#include <configSC.h>
#include <CPropSys.h>
#include <CSC.h>

//------------------------------------------------------------------------------
void TEOM_RK_POLAR(double t, const Vector& y, Vector& yp, void* pAux)
{
  CSC* pSpacecraft = static_cast<CSC*>(pAux);
  pSpacecraft->getPropSys()->EOM_RK_POLAR(t, y, yp, pAux);
};

//------------------------------------------------------------------------------
void TEOM_RK_CARTESIAN(double t, const Vector& y, Vector& yp, void* pAux)
{
  CSC* pSpacecraft = static_cast<CSC*>(pAux);
  pSpacecraft->getPropSys()->EOM_RK_CARTESIAN(t, y, yp, pAux);
};


//------------------------------------------------------------------------------
void REOM_RK_Quaternion_AngVel(double t, const Vector&s, Vector& sd, void* pAux)
{
  CDataReom* data = (CDataReom*) pAux; // data structure
  const double I1(data->moments_of_inertia[0]);
  const double I2(data->moments_of_inertia[1]);
  const double I3(data->moments_of_inertia[2]);
  sd[0] =  0.5*(s[3]*s[4]+s[1]*s[6]-s[2]*s[5]);
  sd[1] =  0.5*(s[3]*s[5]+s[2]*s[4]-s[0]*s[6]);
  sd[2] =  0.5*(s[3]*s[6]+s[0]*s[5]-s[1]*s[4]);
  sd[3] = -0.5*(s[0]*s[4]+s[1]*s[5]+s[2]*s[6]);
  sd[4] =  1.0/I1*(s[5]*s[6]*(I2-I3)+data->steering_moments[0]);
  sd[5] =  1.0/I2*(s[6]*s[4]*(I3-I1)+data->steering_moments[1]);
  sd[6] =  1.0/I3*(s[4]*s[5]*(I1-I2)+data->steering_moments[2]);
}

//------------------------------------------------------------------------------
void REOM_RK_Quaternion(double t, const Vector& s, Vector& sd, void* pAux)
{
  CDataReom* data = (CDataReom*) pAux; // data structure
if (s.size()!=8) flush(cout << "REOM_RK_Quaternion");
  Vector q(s.Slice(0,3));
  Vector q_d(s.Slice(4,7));

  // Correct attitude quaternion
  const double q_corr_fac(NormalizeQuaternion(q));

  // Correct derivative of attitude quaternion
  const double q_d_corr_fac(CorrectQuaternionDerivative(q, q_d));

  // angular velocity vector of body in body frame
  Vector om_b(QuaternionDerivative2AngularVelocityVector(q, q_d));

  // derivative of angular velocity vector in body fame
  Vector om_d_b(3);

  // Euler equations (obtain omega dot from omega, moments of inertia, steering
  // moments)
  EulerEquations(om_b, om_d_b,data->moments_of_inertia,data->steering_moments);

  // set attitude quaternions for DES in quaternion form
  CAttitudeQuaternion q1(q);
  CAttitudeQuaternion q2(om_d_b(0), om_d_b(1), om_d_b(2), -2.0*Dot(q_d,q_d));

  // derivative of quaternion
  Vector q_dd(0.5* ((q1 * q2).getDataVector()));

  sd = Stack(q_d, q_dd);
}



//------------------------------------------------------------------------------
void CPropSys::setInitialPropellantMassInterval(const double min_mass, const double max_mass)
{
  m_min_initialpropellantmass = min_mass;
  m_max_initialpropellantmass = max_mass;
};


//------------------------------------------------------------------------------
void CPropSys::setInitialPropellantMass(const double fraction)
{
  m_initialpropellantmass = (m_min_initialpropellantmass + fraction * (m_max_initialpropellantmass - m_min_initialpropellantmass));
};


//------------------------------------------------------------------------------
void CPropSys::Refuel(void)
{
  m_propellantmass = m_initialpropellantmass;
  m_empty = false;
};
