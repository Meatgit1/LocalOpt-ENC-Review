/*!
\file CProximity.h
\version 1.7
\brief An auxiliary class for proximity (distance and relative velocity
       between two celestial objects) calculations (CProximity).
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 12.02.2006 (AO)

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

#ifndef CPROXIMITY_H
#define CPROXIMITY_H

#include <cmath>
#include <algorithm>
#include <SpaceConst.h>
using namespace std;
namespace SpaceLib
{
  class CState;

  // auxiliary class for handling proximity values between two celestial bodys
  // ----------------------------------------------------------------------------
  /*!
  \brief A class to handle proximity checks between CState objects.

  The class holds information about absolute differences in distance and
  velocity. It may be set/initialized with certain values of the desired
  distance and velocity (direct initialization) or with two CState objects.
  In the second case CProximity calculates the difference in position and
  velocity between the two given states and stores them. In order to get
  useful calculated data the user of CProximity - class has to assure that
  the given CState objects are defined with respect to the same coordinate
  system.\n
  Either values, whether given or calculated, are converted to positive
  values so it makes no sense to provide the class with negative distance
  or velocity data.
  It is also possible to get stored information with respect to given reference
  distances/velocities. The stored or calculated proximity data will be set into
  relation to the reference values given as parameters and returned.
  In this case it is possible to obtain negative return values when negative
  references are used.
  \sa CState
  */
  class CProximity {

  private:

    double r_proximity;         ///< proximity in distance
    double v_proximity;         ///< proximity in velocity
    double a_proximity;         ///< proximity in semimajor axis
    double e_proximity;         ///< proximity in orbit eccentricity
    double i_proximity;         ///< proximity in inclination angle
    double O_proximity;         ///< proximity in longitude of ascending node
    double o_proximity;         ///< proximity in argument of pericenter
    double M_proximity;         ///< proximity in mean anomaly

  public:


    /*!
    \brief Initialize with the difference of two CState objects.

    This function assumes that both states given as arguments are
    expressed with reference to the same coordinate system. The caller
    must assure that this is really the case. Otherwise the result of
    further calls to Proximity member functions will return useless
    results. The absolute values of distance and velocity differences
    between the two states initializes the class' member fields.

    \param left State of body one (should hold data in AUD- system).
    \param right State of body two (should hold data in AUD - system).
    */
    void SetProximity(CState & left, CState & right);



    /*!
    \brief Initialize with given values for distance and velocity proximity.

    SetProximity expects two arguments as the values for distance
    and velocity proximity values for further use. Negative values will be
    converted to positive ones.

    \param dr Desired distance proximity [AU].
    \param dv Desired velocity proximity [AU/d].
    */
    void SetProximity(const double dr, const double dv) {

      r_proximity = fabs(dr);
      v_proximity = fabs(dv);

    };




    /*!
    \brief Sets distance proximity.

    Takes the given argument as new value of distance proximity
    Negative values will be converted to positive ones.

    \param dr Desired distance proximity value [AU].
    */
    void SetDistanceProximity(const double dr) {

      r_proximity = fabs(dr);

    };





    /*!
    \brief Sets velocity proximity.

    Takes the given argument as new value of velocity proximity.
    Possible negative values will be converted to positive ones.

    \param dv Desired velocity difference value [AU/d].
    */
    void SetVelocityProximity(const double dv) {

      v_proximity = fabs(dv);

    };

    void setSemimajorAxisProximity(const double da)
    {
      a_proximity = fabs(da);
    }

    void setEccentricityProximity(const double de)
    {
      e_proximity = fabs(de);
    }

    void setInclinationProximity(const double di)
    {
      i_proximity = di;
    }

    void setNodeProximity(const double dO)
    {
      O_proximity = dO;
    }

    void setPericenterProximity(const double d_o)
    {
      o_proximity = d_o;
    }

    void setMeanAnomalyProximity(const double dM)
    {
      M_proximity = dM;
    }

    void setSemimajorAxisProximity(const double a, const double a_ref);

    void setEccentricityProximity(const double e, const double e_ref);

    void setInclinationProximity(const double i, const double i_ref);

    void setNodeProximity(const double O, const double O_ref);

    void setPericenterProximity(const double o, const double o_ref);

    void setMeanAnomalyProximity(const double M, const double M_ref);

    void setSemimajorAxisProximity(const CState& left, const CState& right)
    {
      setSemimajorAxisProximity(left.a(), right.a());
    }

    void setEccentricityProximity(const CState& left, const CState& right)
    {
      setEccentricityProximity(left.e(), right.e());
    }

    void setInclinationProximity(const CState& left, const CState& right)
    {
      setInclinationProximity(left.i(), right.i());
    }

    void setNodeProximity(const CState& left, const CState& right)
    {
      setNodeProximity(left.Omega(), right.Omega());
    }

    void setPericenterProximity(const CState& left, const CState& right)
    {
      setPericenterProximity(left.omega(), right.omega());
    }

    void setMeanAnomalyProximity(const CState& left, const CState& right)
    {
      setMeanAnomalyProximity(left.M(), right.M());
    }


    /*!
    \brief Returns stored distance proximity value [AU].

    The returned value is always a positive one.

    \return Stored value of distance difference [AU].
    */
    double GetDistanceProximity(void) const {

      return r_proximity;

    };




    /*!
    \brief Returns difference in distance between two states and stores it.

    The absolute distance between two CState objects given as arguments will
    be calculated and stored in a local variable. So in every case the
    returned distance value is positive.
    \em Attention:
    In order to obtain useful data, the caller has to assure that both
    CState - object refer to the same coordinate system and hold data in
    AUD unit system.

    \param left CState object one.
    \param right CState object two.
    \return Difference in distance [AU] (positive value).
    */
    double GetDistanceProximity(CState & left, CState & right);





     /*!
     \brief Returns the stored velocity proximity value [AU/d].

     The returned value is the stored difference in velocity.

     \return Stored value of velocity difference [AU/d].
     */
    double GetVelocityProximity(void) const {
      return v_proximity;
    };





    /*!
    \brief Returns the absolute difference in velocity between two states
    and stores it in local data fields [AU/d].

    The relative velocity between the given argument CState objects will be
    calculated and stored in a local variable. In every case the returned
    value is positive.

    \param left CState object one.
    \param right CState object two.
    \return Difference in velocity [AU/d] (positive value).
    */
    double GetVelocityProximity(CState & left, CState & right);

    double getSemimajorAxisProximity(SpaceLib::UnitType type = AUD) const
    {
      return (type==AUD?1.0:AU2KM)*a_proximity;
    };

    double getEccentricityProximity(void) const { return e_proximity;}

    double getInclinationProximity(void) const { return i_proximity;};

    double getNodeProximity(void) const { return O_proximity;};

    double getPericenterProximity(void) const { return o_proximity;};

    double getMeanAnomalyProximity(void) const { return M_proximity;};


     /*!
     \brief Returns the stored distance difference value [km].

     The function returns the stored value of distance proximity.
     Although the expression .._SI leads to the assumtion that the
     return value is expressed in meter, it's kilometer in reality.
     So the caller has to convert it to meter if he wants to.

     \return Stored value of distance proximity [km].
     \sa GetDistanceProximity, KM2M
     */
    double GetDistanceProximity_SI(void) const {
      return r_proximity * AU2KM;
    };

    double getSemimajorAxisProximitySI(void) const
    {
      return a_proximity * AU2KM;
    }



    /*!
    \brief Returns the difference in distance between two states and
    stores it in local data fields [km].

    The distance between the given argument CState objects will be calculated
    and stored in a local variable. In every case the returned distance
    value is positive.\n
    Although the expression .._SI leads to the assumtion that the
    return value is expressed in meter, it's kilometer in reality.
    So, eventually, the caller has to convert it to desired units.

    \param left CState object one.
    \param right CState object two.
    \return Difference in distance [km] (positive value).
    */
    double GetDistanceProximity_SI(CState & left, CState & right);


    double getSemimajorAxisProximity(const double a, const double a_ref,
                                     SpaceLib::UnitType type = AUD);

    double getEccentricityProximity(const double e, const double e_ref);

    double getInclinationProximity( const double i, const double i_ref);

    double getNodeProximity( const double O, const double O_ref);

    double getPericenterProximity( const double o, const double o_ref);

    double getMeanAnomalyProximity( const double M, const double M_ref);

    double getSemimajorAxisProximity(const CState& left, const CState& right,
                                     const SpaceLib::UnitType type);

    double getEccentricityProximity(const CState& left, const CState& right);

    double getInclinationProximity(const CState& left, const CState& right);

    double getNodeProximity(const CState& left, const CState& right);

    double getPericenterProximity(const CState& left, const CState& right);

    double getMeanAnomalyProximity(const CState& left, const CState& right);

    double getNormedSemimajorAxisProximity(const double max_da)
    {
      return (a_proximity / max(1.0e-40, max_da));
    }

    double getNormedEccentricityProximity(const double max_de)
    {
      return e_proximity / max(1.0e-40, max_de);
    }

    double getNormedInclinationProximity(const double max_di)
    {
      return i_proximity / max(1.0e-40, max_di);
    }

    double getNormedNodeProximity(const double max_dO)
    {
      return O_proximity / max(1.0e-40, max_dO);
    }

    double getNormedPericenterProximity(const double max_do)
    {
      return o_proximity / max(1.0e-40, max_do);
    }

    double getNormedMeanAnomalyProximity(const double max_dM)
    {
      return M_proximity / max(1.0e-40, max_dM);
    }

    double getNormedSemimajorAxisProximity(
      const double a, const double a_ref, const double max_da)
    {
      setSemimajorAxisProximity(a, a_ref);
      return getNormedSemimajorAxisProximity(max_da);
    }

    double getNormedEccentricityProximity(
      const double e, const double e_ref, const double max_de)
    {
      setEccentricityProximity(e, e_ref);
      return getNormedEccentricityProximity(max_de);
    }

    double getNormedInclinationProximity(
      const double i, const double i_ref, const double max_di)
    {
      setInclinationProximity(i, i_ref);
      return getNormedInclinationProximity(max_di);
    }

    double getNormedNodeProximity(
      const double O, const double O_ref, const double max_dO)
    {
      setNodeProximity(O, O_ref);
      return getNormedNodeProximity(max_dO);
    }

    double getNormedPericenterProximity(
      const double o, const double o_ref, const double max_do)
    {
      setPericenterProximity(o, o_ref);
      return getNormedPericenterProximity(max_do);
    }

    double getNormedMeanAnomalyProximity(
      const double M, const double M_ref, const double max_dM)
    {
      setMeanAnomalyProximity(M, M_ref);
      return getNormedMeanAnomalyProximity(max_dM);
    }

    /*!
    \brief The stored velocity proximity value [km/s].

    The returned value is the stored velocity proximity [km/s].\n
    Although the expression .._SI leads to the assumtion that the
    return value is expressed in meter per second, it's kilometer
    per second in reality.
    So, eventually, the caller has to convert it to desired units.

    \return Stored value of velocity difference [km/s].
    \sa GetVelocityProximity, KMS2MS
    */
    double GetVelocityProximity_SI(void) const {
      return v_proximity * AUD2KMS;
    };


    /*!
    \brief Returns the difference in velocity between two states and stores
    it in local data fields [km/s].

    The relative absolute velocity between the given argument CState objects will be
    calculated and stored in a local variable. In every case the returned
    value is positive.\n
    Although the expression .._SI leads to the assumtion that the
    return value is expressed in meter per second, it's kilometer
    per second in reality.
    So, eventually, the caller has to convert it to desired units.

    \param left CState object one.
    \param right CState object two.
    \return Difference in velocity [km/s] (positive value).
    */
    double GetVelocityProximity_SI(CState & left, CState & right);






    /*!
    \brief Returns the stored distance difference with respect to the given
    value.

    Though data held within CProximity is expressed in AU or AU/d yous should
    also use values in those units to get useful data.\n
    \em Attention: Depending on the argument, the result may be negative.

    \param dr_max The distance value to which the stored one should be normalized to.
    \return Stored distance proximity with respect to given value.
    */
    double GetNormedDistanceProximity(const double dr_max) const {
      return r_proximity / dr_max;
    };





    /*!
    \brief Distance proximity between the given CState objects with respect
    to the given value.

    \em Attention: Depending on the argument, the result may be negative.\n
    It should be assured, that the norm dr_max is expressed in [AU].
    The calculated absolute difference will be stored in a local data field.

    \param left  Reference to State one.
    \param right Reference to State two.
    \param dr_max The distance value to which the calculated one should be
    normalized to [AU].
    \return Distance between the two argument CState object with respect
    to given value.
    */
    double GetNormedDistanceProximity(
     CState & left,
     CState & right,
     const double dr_max);





     /*!
     \brief Returns the stored velocity proximity with respect to the
     given value.

     \em Attention: Depending on the argument, the result may be negative.\n
     It should be assured, that the norm dv_max is expressed in [AU/d].
     The calculated absolute difference will be stored in a local data field.

     \param dv_max The velocity to which the stored one should be
     normalized with.
     \return Stored velocity value with respect to given value.
     */
    double GetNormedVelocityProximity(const double dv_max) const {

      return v_proximity / dv_max;

    };


    /*!
    \brief Relative velocity between the given CState objects with respect
    to the given value.

    \em Attention: Depending on the argument, the result may be negative.\n
    It should be assured, that the norm dv_max is expressed in [AU/d].
    The calculated absolute velocity will be stored in a local data field.

    \param left CState object one.
    \param right CState object two.
    \param dv_max The velocity value to which the calculated one should be normalized to [AU/d].

    \return Velocity difference between the two argument CState object with respect to given value.
    */
    double GetNormedVelocityProximity(
      CState & left,
      CState & right,
      const double dv_max
    );



    /*!
    \brief Square of vector norm of normalized distance and velocity proximities

    It should be assured, that the norm dr_max is expressed in [AU] and
    dv_max in [AU/d].

    \param dv_max The velocity value to normalize the stored one with [AU/d].
    \param dr_max The distance value to normalize the stored one with [AU].
    \return Vector norm square of stored velocity and distance proximities with
    respect to given reference values.
    */
    double GetNormedProximitySquare(const double dr_max, const double dv_max) const throw();

    /*!
    \brief Vector norm of normalized distance and velocity proximities.

    It should be assured, that the norm dr_max is expressed in [AU] and
    dv_max in [AU/d].

    \param dv_max The velocity value to normalize the stored one with [AU/d].
    \param dr_max The distance value to normalize the stored one with [AU].
    \return Vector norm of stored velocity and distance proximities with
    respect to given reference values.
    */
    double GetNormedProximity(const double dr_max, const double dv_max) const;





    /*!
    \brief Vector norm of normalized distance and velocity differences between
    given CState object.

    It should be assured, that the norm dr_max is expressed in [AU] and
    dv_max in [AU/d].\n
    \em Attention: CState data should be expressed in AUD - system.

    \param left CState object one.
    \param right CState object two.
    \param dv_max The velocity value to normalize the stored one with [AU/d].
    \param dr_max The distance value to to normalize the stored one with [AU].
    \return Vector norm of calculated velocity and distance differences
    with respect to given values.
    */
    double GetNormedProximity(
     CState & left,
     CState & right,
     const double dr_max,
     const double dv_max
    );





    /*!
    \brief Check function to see if the stored values are lower or equal than
    given maximal values.

    \param dr_max Upper boundary for distance proximity [AU].
    \param dv_max Upper boundary for velocity proximity [AU/d].
    \return True if the stored values are \em both lower/equal than the
    arguments.
    */
    bool CheckProximity(const double dr_max, const double dv_max) const;




    /*!
    \brief Check function to see if the calculated distance and velocity
    differences are lower or equal than the given maximal values.

    It should be assured, that the norm values dr_max and dv_max are expressed
    in [AU] and [AU/d]. \n
    The calculated differences in distance and velocity will be stored in
    local data fields as new proximity values.

    \param left CState object one.
    \param right CState object two.
    \param dr_max Upper boundary for distance proximity [AU].
    \param dv_max Upper boundary for velocity proximity [AU/d].
    \return True if the calculated differences are both lower/equal than the arguments.
    */
    bool CheckProximity(
     CState & left,
     CState & right,
     const double dr_max,
     const double dv_max
    );



    double GetImpactVelocity(
     CState & left,
     CState & right
     );

	// Mixed velocity calculation for the ACT(ESA) Global Optimization Competition
	double GetActVelocity(
     CState & left,
     CState & right
     );


  };                            // class CProximity

  // implementation

  //----------------------------------------------------------------------------
  inline double CProximity::GetNormedProximitySquare(
    const double dr_max, const double dv_max) const throw()
  {
    return ( pow(r_proximity/dr_max,int(2)) + pow(v_proximity/dv_max, int(2)) );
  } // double CProximity::GetNormedProximitySquare(const double, const double)

  //----------------------------------------------------------------------------
  inline double CProximity::GetNormedProximity(
    const double dr_max, const double dv_max) const throw()
  {
    return sqrt(GetNormedProximitySquare(dr_max, dv_max));
  } // double CProximity::GetNormedProximity(const double, const double) const

} // namespace SpaceLib

#endif // CPROXIMITY_H

