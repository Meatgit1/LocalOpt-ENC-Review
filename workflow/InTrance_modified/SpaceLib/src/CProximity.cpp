/*!
\file CProximity.cpp
\version 1.7
\brief An auxiliary class for proximity (distance and relative velocity
       between two celestial objects) calculations (CProximity).
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 22.10.2006 (AO)

\b Purpose:
   The class CProximity is intended to simplify the need
   for comparing distance and velocity differences between
   two bodies in space. For the scope of SpaceLib state
   information of bodies in space are almost often expressed
   through CState classes. So it is straightforward to use
   CProximity in connection with them. That's why one can
   initialize a CProximity class with two CState objects but
   there are also initializer functions for using double
   values.
   In general the unit system used is AUD.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.

\b Copyright:
   (c) 2001-2006 Bernd Dachwald, Andreas Ohndorf
*/

#include <math.h>

#include <MLGlobal.h>
#include <CState.h>
#include <CProximity.h>

using namespace std;
using namespace MathLib;
using namespace SpaceLib;

// ------------------------------------------------------------------------------
void CProximity::SetProximity(CState & left, CState & right)
{

  double sum;
  sum = pow((left.x()-right.x()), int(2));
  sum += pow((left.y()-right.y()), int(2));
  sum += pow((left.z()-right.z()), int(2));
  r_proximity = sqrt(sum);

  sum = pow((left.dx()-right.dx()), int(2));
  sum += pow((left.dy()-right.dy()), int(2));
  sum += pow((left.dz()-right.dz()), int(2));
  v_proximity = sqrt(sum);
} // void SetProximity(CState & left, CState & right)

// -----------------------------------------------------------------------------
void CProximity::setSemimajorAxisProximity(const double a, const double a_ref)
{
  a_proximity = fabs(a_ref-a);
} // void setSemimajorAxisProximity(const double a, const double a_ref)

// -----------------------------------------------------------------------------
void CProximity::setEccentricityProximity(const double e, const double e_ref)
{
  e_proximity = fabs(e_ref-e);
} // void setEccentricityProximity(const double e, const double e_ref)

// -----------------------------------------------------------------------------
void CProximity::setInclinationProximity(const double i, const double i_ref)
{
  i_proximity = fabs(i_ref-i);
} // void setInclinationProximity(const double i, const double i_ref)

// -----------------------------------------------------------------------------
void CProximity::setNodeProximity(const double O, const double O_ref)
{
  O_proximity = fabs(O_ref-O);
} // void setNodeProximity(const double O, const double O_ref)

// -----------------------------------------------------------------------------
void CProximity::setPericenterProximity(const double o, const double o_ref)
{
  o_proximity = fabs(o_ref-o);
} // void setPericenterProximity(const double o, const double o_ref)

// -----------------------------------------------------------------------------
void CProximity::setMeanAnomalyProximity(const double M, const double M_ref)
{
  M_proximity = fabs(M_ref-M);
} // void setMeanAnomalyProximity(const double M, const double M_ref)

// -----------------------------------------------------------------------------
double CProximity::getSemimajorAxisProximity(const double a, const double a_ref,
                                             SpaceLib::UnitType type)
{
  setSemimajorAxisProximity(a, a_ref);
  return ( (type==AUD?1.0:AU2KM) * getSemimajorAxisProximity() );
} // double getSemimajorAxisProximity(const double,const double,const UnitType)


// -----------------------------------------------------------------------------
double CProximity::getSemimajorAxisProximity(const CState& left, 
                                             const CState& right,
                                             const SpaceLib::UnitType type)
{
  setSemimajorAxisProximity(left, right);
  return getSemimajorAxisProximity(type);
} // double getSemimajorAxisProximity(const CState&,const CState&,const UnitType)

// -----------------------------------------------------------------------------
double CProximity::getEccentricityProximity(const CState& left,
                                            const CState& right)
{
  setEccentricityProximity(left, right);
  return getEccentricityProximity();
} //double getEccentricityProximity(const CState& left, const CState& right)


// -----------------------------------------------------------------------------
double CProximity::getInclinationProximity(const CState& left,
                                           const CState& right)
{
  setInclinationProximity(left, right);
  return getInclinationProximity();
} // double getInclinationProximity(const CState& left,const CState& right)

// -----------------------------------------------------------------------------
double CProximity::getNodeProximity(const CState& left,
                                 const CState& right)
{
  setNodeProximity(left, right);
  return getNodeProximity();
} // double CProximity::getNodeProximity(const CState& left, const CState& right)

// -----------------------------------------------------------------------------
double CProximity::getPericenterProximity(const CState& left,
                                          const CState& right)
{
  setPericenterProximity(left, right);
  return getPericenterProximity();
} // double getPericenterProximity(const CState& left,const CState& right)

// -----------------------------------------------------------------------------
double CProximity::getMeanAnomalyProximity(const CState& left,
                                           const CState& right)
{
  setMeanAnomalyProximity(left, right);
  return getMeanAnomalyProximity();
} // double getMeanAnomalyProximity(const CState& left, const CState& right)

// -----------------------------------------------------------------------------
double CProximity::GetDistanceProximity(CState & left, CState & right)
{
  SetProximity(left, right);
  return GetDistanceProximity();
} // double GetDistanceProximity(CState & left, CState & right)


// ----------------------------------------------------------------------------
double CProximity::GetVelocityProximity(CState & left, CState & right)
{
  SetProximity(left, right);
  return GetVelocityProximity();
} // double GetVelocityProximity(CState & left, CState & right)


// -----------------------------------------------------------------------------
double CProximity::GetDistanceProximity_SI(CState & left, CState & right)
{
  SetProximity(left, right);
  return GetDistanceProximity_SI();
} // double GetDistanceProximity_SI(CState & left, CState & right)


// -----------------------------------------------------------------------------
double CProximity::GetVelocityProximity_SI(CState & left, CState & right)
{
  SetProximity(left, right);
  return GetVelocityProximity_SI();
} // double GetVelocityProximity_SI(CState & left, CState & right)


// -----------------------------------------------------------------------------
double CProximity::GetNormedDistanceProximity(
 CState & left,
 CState & right,
 const double dr_max
)
{
  SetProximity(left, right);
  return GetNormedDistanceProximity(dr_max);
} // double GetNormedDistanceProximity(CState & ,CState & ,const double)

// -----------------------------------------------------------------------------
double CProximity::GetNormedVelocityProximity(
 CState & left,
 CState & right,
 const double dv_max
)
{
  SetProximity(left, right);
  return GetNormedVelocityProximity(dv_max);
} // double GetNormedVelocityProximity(CState & left, CState & right, const double dv_max)

//// -----------------------------------------------------------------------------
//double CProximity::GetNormedProximity(
// const double dr_max,
// const double dv_max
//) const
//{
//  const double dr_norm = r_proximity / dr_max;
//  const double dv_norm = v_proximity / dv_max;
//  return sqrt(dr_norm * dr_norm + dv_norm * dv_norm);
//} // double CProximity::GetNormedProximity(const double dr_max, const double  v_max ) const 

// -----------------------------------------------------------------------------
double CProximity::GetNormedProximity(
 CState & left, CState & right,
 const double dr_max,
 const double dv_max
) {

  SetProximity(left, right);
  return GetNormedProximity(dr_max, dv_max);

} // double GetNormedProximity(...)

// -----------------------------------------------------------------------------
bool CProximity::CheckProximity(
 const double dr_max,
 const double dv_max
) const 
{
  if (r_proximity <= dr_max && v_proximity <= dv_max) return true;
  else return false;
} // bool CheckProximity(const double dr_max, const double dv_max) const 

// -----------------------------------------------------------------------------
bool CProximity::CheckProximity(
 CState & left,
 CState & right,
 const double dr_max,
 const double dv_max
) {
  SetProximity(left, right);
  return CheckProximity(dr_max, dv_max);
} // bool CProximity::CheckProximity(...)

// ----------------------------------------------------------------------------
double CProximity::GetImpactVelocity(CState & left, CState & right)
{
  /*const double dx = left.dx() - right.dx();
  const double dy = left.dy() - right.dy();
  const double dz = left.dz() - right.dz();

  return sqrt(dx*dx+dy*dy+dz*dz);*/

  const Vector v_left(left.dx(),left.dy(),left.dz());
  const Vector v_right(-right.dx(),-right.dy(),-right.dz());
  return Dot(v_left,v_right);
} // double GetImpactVelocity(CState & left, CState & right)

// -----------------------------------------------------------------------------
double CProximity::GetActVelocity(CState & sc, CState & target)
{
  const Vector v_sc(sc.dx(),sc.dy(),sc.dz());
  const Vector v_target(target.dx(),target.dy(),target.dz());

  return fabs(Dot((v_sc-v_target),v_target));
} // double GetActVelocity(CState & sc, CState & target)

// end of implementation

