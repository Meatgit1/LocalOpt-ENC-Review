/*!
\file CPhysUnit.cpp
\version 0.1
\brief Implementation of CPhysUnit class and global variables.
between them.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 28.09.2006 (AO)

\b Purpose:
   
\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.
   This code is based on the code from:
   Montenbruck, Gill: Satellite Orbits. Springer, 2000 and
   Montenbruck, Pfleger: Astronomie mit dem Personal Computer (3. Aufl.). Springer, 1999

\b Copyright:
   (c) 2001-2006 Bernd Dachwald, Andreas Ohndorf
*/

// STL
#include <vector>

#include <CPhysUnit.h>
#include <SLTime.h>
#include <MLConst.h>
#include <SpaceConst.h>

namespace SpaceLib
{
  // definition of global unit variables

  CPhysUnit u_nothing(CPhysUnit::UNDEFINED_UNIT, "", 0.0, 1.0, 0.0, 0);
  CPhysUnit* U_NO_UNIT = &u_nothing;

  // length variables
  // meter [m]
  CPhysUnit u_m(CPhysUnit::LENGTH_UNIT , "M" , 0.0, 1.0e+0 , 0.0, 0);
  CPhysUnit* U_METER = &u_m;

  // decimeter [dm]
  CPhysUnit u_dm(CPhysUnit::LENGTH_UNIT, "DM", 0.0, 1.0e+1 , 0.0, U_METER);
  CPhysUnit* U_DECIMETER = &u_dm;

  // centimeter [cm]
  CPhysUnit u_cm(CPhysUnit::LENGTH_UNIT, "CM", 0.0, 1.0e+2 , 0.0, U_METER);
  CPhysUnit* U_CENTIMETER = &u_cm;

  // millimeter [mm]
  CPhysUnit u_mm(CPhysUnit::LENGTH_UNIT, "MM", 0.0, 1.0e+1 , 0.0, U_CENTIMETER);
  CPhysUnit* U_MILLIMETER = &u_mm;

  // kilometer [km]
  CPhysUnit u_km(CPhysUnit::LENGTH_UNIT, "KM", 0.0, 1.0e-3 , 0.0, U_METER); 
  CPhysUnit* U_KILOMETER = &u_km;

  // astronomical unit
  CPhysUnit u_au(CPhysUnit::LENGTH_UNIT, "AU", 0.0, KM2AU  , 0.0, U_KILOMETER);
  CPhysUnit* U_AU = &u_au;

  // inch [in]
  CPhysUnit u_in(CPhysUnit::LENGTH_UNIT, "IN", 0.0, 1/2.54 , 0.0, U_CENTIMETER);
  CPhysUnit* U_IN = &u_in;

  // foot [ft]
  CPhysUnit u_ft(CPhysUnit::LENGTH_UNIT, "FT", 0.0, 1/12.0   , 0.0, U_IN);
  CPhysUnit* U_FT = &u_ft;

  // yard [yd]
  CPhysUnit u_yd(CPhysUnit::LENGTH_UNIT, "YD", 0.0, 1/3.0    , 0.0, U_FT);
  CPhysUnit* U_YD = &u_yd;


  // time units
  // second [sec]
  CPhysUnit u_s(CPhysUnit::TIME_UNIT, "SEC", 0.0, 1.0, 0.0, 0);
  CPhysUnit* U_SECOND = &u_s;

  // minute [min]
  CPhysUnit u_min(CPhysUnit::TIME_UNIT, "MIN", 0.0, 1.0/60.0, 0.0, U_SECOND);
  CPhysUnit* U_MINUTE = &u_min;

  // hour [h]
  CPhysUnit u_hour(CPhysUnit::TIME_UNIT, "HOUR", 0.0, 1.0/60.0, 0.0, U_MINUTE);
  CPhysUnit* U_HOUR = &u_hour;

  // day [d]
  CPhysUnit u_day(CPhysUnit::TIME_UNIT, "DAY", 0.0, 1.0/24.0, 0.0, U_HOUR);
  CPhysUnit* U_DAY = &u_day;

  CPhysUnit u_jyear(CPhysUnit::TIME_UNIT, "JYR", 0.0, SpaceLib::DAY2JYR
    , 0.0, U_DAY);
  CPhysUnit* U_JULIAN_YEAR = &u_jyear;


  // linear velocity
  // meter per second [m/s]
  CPhysUnit u_m_per_sec(CPhysUnit::LINEAR_VELOCITY_UNIT, "M/S", 0.0, 1.0, 0.0
    , 0);
  CPhysUnit* U_M_PER_SEC = &u_m_per_sec;

  // kilometer per second [km/s]
  CPhysUnit u_km_per_sec(CPhysUnit::LINEAR_VELOCITY_UNIT, "KM/S", 0.0
    , 1.0e-3, 0.0, U_M_PER_SEC);
  CPhysUnit* U_KM_PER_SEC = &u_km_per_sec;

  // kilometer per hour [km/h]
  CPhysUnit u_km_per_h(CPhysUnit::LINEAR_VELOCITY_UNIT, "KM/H", 0.0
    , 3.6, 0.0, U_M_PER_SEC);
  CPhysUnit* U_KM_PER_H = &u_km_per_h;

  // au per day [AU/D]
  CPhysUnit u_au_per_day(CPhysUnit::LINEAR_VELOCITY_UNIT, "AU/D", 0.0
    , M2AU/SEC2DAY, 0.0, U_M_PER_SEC);
  CPhysUnit* U_AU_PER_DAY = &u_au_per_day;


  // angle
  // radian
  CPhysUnit u_rad(CPhysUnit::ANGLE_UNIT, "RAD", 0.0, 1.0, 0.0, 0);
  CPhysUnit* U_RADIAN = &u_rad;

  // degree
  CPhysUnit u_deg(CPhysUnit::ANGLE_UNIT, "DEG", 0.0, MathLib::RAD2DEG, 0.0
    , U_RADIAN);
  CPhysUnit* U_DEGREE = &u_deg;


  // angular velocity
  // radian/day [rad/d]
  CPhysUnit u_rad_per_day(CPhysUnit::ANGULAR_VELOCITY_UNIT, "RAD/D", 0.0, 1.0
    , 0.0, 0);
  CPhysUnit* U_RAD_PER_DAY = &u_rad_per_day;

  // degree per day
  CPhysUnit u_deg_per_day(CPhysUnit::ANGULAR_VELOCITY_UNIT, "DEG/D", 0.0
    , MathLib::RAD2DEG , 0.0, U_RAD_PER_DAY);
  CPhysUnit* U_DEG_PER_DAY = &u_deg_per_day;


  //---------------------------------------------------------------------------
  CPhysUnit* Units[] =
  {
    U_NO_UNIT

    ,U_METER, U_DECIMETER, U_CENTIMETER, U_MILLIMETER, U_KILOMETER
    , U_AU, U_IN, U_FT, U_YD

    ,U_SECOND, U_MINUTE, U_HOUR, U_DAY, U_JULIAN_YEAR

    ,U_M_PER_SEC, U_KM_PER_SEC, U_KM_PER_H, U_AU_PER_DAY

    ,U_RADIAN, U_DEGREE

    ,U_RAD_PER_DAY, U_DEG_PER_DAY
  };
  const int N_UNITS = sizeof(Units)/sizeof(Units[0]);

  

  //---------------------------------------------------------------------------
  CPhysUnit* GetUnit(const std::string& s_unit)
  {
    int i = 0;
    while (i<N_UNITS)
    {
      if ((*Units[i])==s_unit) return Units[i];
      ++i;
    }
    return U_NO_UNIT;
  } // CPhysUnit* GetUnit(const std::string& s_unit)


  //---------------------------------------------------------------------------
  double CPhysUnit::convertTo( 
    const double value 
    ,const CPhysUnit* target_unit) const
  { 
    if (type == target_unit->type) 
    { 
      double t1 = 0.0;
      if (base_unit == target_unit->base_unit)
      { // same base unit
        t1 = (((value-post_add)/multiplier)-pre_add);
        return ((t1+target_unit->pre_add)*target_unit->multiplier+
          target_unit->post_add);
      }
      else if (base_unit != 0)
      { // not already at base unit so do the descent deeper towards base
        t1 = (((value-post_add)/multiplier)-pre_add);
        return base_unit->convertTo(t1, target_unit);
      }
      else // base_unit == 0 => base
      {
        t1 = convertTo(value, target_unit->base_unit);
        return ((t1+target_unit->pre_add)*target_unit->multiplier+
          target_unit->post_add);
      }
    } 
    else return value;
  } // double convertTo(const double value, const CPhysUnit* target_unit) const

} // namespace SpaceLib