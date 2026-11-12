/*!
\file CPhysUnit.h
\version 0.1
\brief Class interface for incorporating physical units and conversion
between them.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 03.10.2006 (AO)

\b Purpose:
   The class CPhysUnit provides easy means for storing physical unit
   types and conversion operations between them.
   
\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.

\b Copyright:
   (c) 2001-2006 Bernd Dachwald, Andreas Ohndorf
*/

#ifndef CPHYSUNIT_H
#define CPHYSUNIT_H

/// STL
#include <string>

namespace SpaceLib
{

  /// physical unit representation and conversion class
  class CPhysUnit : public std::string
  { 
  public:

    /// supported unit types
    enum PhysUnitType 
    { 
      UNDEFINED_UNIT          ///< standard unit type
      ,TIME_UNIT 
      ,MASS_UNIT 
      ,LENGTH_UNIT 
      ,ANGLE_UNIT
      ,LINEAR_VELOCITY_UNIT
      ,ANGULAR_VELOCITY_UNIT
      ,FORCE_UNIT 
      ,TEMPERATURE_UNIT
      ,PRESSURE_UNIT
    };

  private: 
    PhysUnitType type;    ///< current type of unit
    double pre_add;       ///< value to add before multiplying
    double post_add;      ///< value to add after multiplying
    double multiplier;    ///< conversion factor
    CPhysUnit* base_unit; ///< base unit on that the conversion is done

  public: 

    /// standard constructor 
    CPhysUnit(void) : type(UNDEFINED_UNIT), base_unit(0) {};

    // copy constructor 
    CPhysUnit(const CPhysUnit& unit); 

    // initializer constructor 
    CPhysUnit( 
      const PhysUnitType unit_type ///< the desired unit type
      , const std::string s_unit   ///< the literal unit description (e.g. "m")
      , const double pre_add
      , const double multiplier 
      , const double post_add 
      , CPhysUnit* base 
      ); 

    /// converts value from class unit to target unit 
    virtual double convertTo(const double value
      , const CPhysUnit* target_unit) const; 

    /// destructor 
    ~CPhysUnit(void) {}; 
  }; // CPhysUnit



  //--------------------------------------------------------------------------- 
  inline CPhysUnit::CPhysUnit(const CPhysUnit& unit) 
    : std::string(unit) 
    , type(unit.type) 
    , pre_add(unit.pre_add) 
    , post_add(unit.post_add) 
    , multiplier(unit.multiplier) 
    , base_unit(unit.base_unit) 
  {} 


  //--------------------------------------------------------------------------- 
  inline CPhysUnit::CPhysUnit( 
    const PhysUnitType unit_type 
    , const std::string s_unit 
    , const double _pre_add 
    , const double _multiplier 
    , const double _post_add 
    , CPhysUnit* base 
    )
    : std::string(s_unit)
    , type(unit_type)
    , pre_add(_pre_add)
    , multiplier(_multiplier)
    , post_add(_post_add)
    , base_unit(base)
  {}

  // physical unit (object)
  typedef CPhysUnit* PhysUnit;

  extern PhysUnit U_NO_UNIT;
  extern PhysUnit U_METER;
  extern PhysUnit U_DECIMETER;
  extern PhysUnit U_CENTIMETER;
  extern PhysUnit U_MILLIMETER;
  extern PhysUnit U_KILOMETER;
  extern PhysUnit U_AU;
  extern PhysUnit U_INCH;
  extern PhysUnit U_FOOT;
  extern PhysUnit U_YARD;

  extern PhysUnit U_SECOND;
  extern PhysUnit U_MINUTE;
  extern PhysUnit U_HOUR;
  extern PhysUnit U_DAY;
  extern PhysUnit U_JULIAN_YEAR;

  extern PhysUnit U_M_PER_SEC;
  extern PhysUnit U_KM_PER_SEC;
  extern PhysUnit U_KM_PER_H;
  extern PhysUnit U_AU_PER_DAY;

  extern PhysUnit U_DEGREE;
  extern PhysUnit U_RADIAN;

  extern PhysUnit U_RAD_PER_DAY;
  extern PhysUnit U_DEG_PER_DAY;

  /*! \brief get unit object of given unit string 

    \params_unit string description of desired unit
    \return CPhysUnit object
  */
  CPhysUnit* GetUnit(const std::string& s_unit);

  /// output operator for physical units
  std::ostream& operator<<(std::ostream& os, CPhysUnit* unit);


  inline std::ostream& operator<<(std::ostream& os, CPhysUnit* unit)
  {
    return (os << (*unit));
  }

} // namespace SpaceLib

#endif // CPHYSUNIT_H