/*!
\file CPropSys.h
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

\author

  Andreas Ohndorf\n
  German Aerospace Center (DLR)\n
  Institute of Space Simulation\n
  Linder Hoehe\n
  51147 Cologne\n
  Germany\n
  andreas.ohndorf@gmx.de

\date

  31.10.2006 (AO)

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
#ifndef CPROPSYS_H_INCLUDED
#define CPROPSYS_H_INCLUDED

// standard C/C++ library
#include <new>
#include <math.h>

// MathLib
#include <MLMatrix.h>

// SpaceLib
#include <CBody.h>
#include <CSolarSystem.h>

// SCLib
#include <SCConst.h>

const SpaceLib::CoordType DFLT_EOM_COORD_TYPE = POLAR;

enum RDES_Type
{
  DES_Q_QD,
  DES_Q_OM
};

// typedef for equations of motion function
typedef void (*EOM_RK_func) (double x, const Vector & y, Vector & yp, void *pAux);

/// The equations of motion (EOM) for spacecraft in Runge-Kutta-form (RK-form).
/** \param t actual time
    \param y actual state vector
    \param yp first derivation of the actual state vector
    \param pAux pointer to auxiliary data (must point to the CSpacecraft-object) */
void TEOM_RK_POLAR(double t,const MathLib::Vector& y,MathLib::Vector& yp
                   , void* pAux);

void TEOM_RK_CARTESIAN(double t, const MathLib::Vector& y, MathLib::Vector& yp
                       , void* pAux);

/// Rotational equations of motion in Runga-Kutta-form and quaternion
/// formulation.
/** \param t actual time
    \param s actual rotational state vector
    \param sp first derivation of the actual rotational state vector
    \param pAux pointer to auxiliary data (must point to the CSpacecraft-object) */
void REOM_RK_Quaternion(double t, const Vector& s, Vector& sp, void* pAux);

void REOM_RK_Quaternion_AngVel(double t, const Vector&s, Vector& sd,void* pAux);


/*! \brief Eulers equations of motion (derivative of angular momentum in
principalaxes frame).

  \param om Angular velocity vector of body (principal axes body frame)
  \param om_dot Derivative with respect to time of angular velocity vector
  omega.
  \param moi Moments in inertia in body frame about principle axes.
  \param sm Steering moments about principle axes.
  \return Derivative of angular velocity vector in modified parameter om_d*/
inline void EulerEquations(const MathLib::Vector& om
                          ,      MathLib::Vector& om_dot
                          ,const MathLib::Vector& moi
                          ,const MathLib::Vector& sm);


/// data structure used within integration of rotational equations of motion
typedef struct CDataStructureReom
{
  Vector steering_moments;       ///< steering moments about principal body axis
  Vector max_steering_moments;                     ///< maximum steering moments
  Vector moments_of_inertia;                   ///< principal moments of inertia

  CDataStructureReom(void);                             ///< default constructor

  // copy constructor
  CDataStructureReom(const CDataStructureReom& init)
    : steering_moments(init.steering_moments)
    , max_steering_moments(init.max_steering_moments)
    , moments_of_inertia(init.moments_of_inertia) 
    {};

  // destructor
  ~CDataStructureReom(void) {};

  // copy operator
  const CDataStructureReom& operator=(const CDataStructureReom& rhs)
  {
    steering_moments = rhs.steering_moments;
    max_steering_moments = rhs.max_steering_moments;
    moments_of_inertia = rhs.moments_of_inertia;
    return (*this);
  }

} CDataReom;

//------------------------------------------------------------------------------
inline CDataStructureReom::CDataStructureReom(void)
: steering_moments(3), max_steering_moments(3), moments_of_inertia(3)
{}


/// Propulsion system base class.
class CPropSys
{
  friend class CSC;
             
protected: 
  CState* m_state;                        ///< pointer to the spacecraft (state)

private:   
    PropulsionSystemType m_type;          ///< spacecraft propulsion system type
    SpaceLib::CoordType m_eom_type;  ///< coordinate type of equations of motion

public:    
  
  /// type of incorporated propulsion system type
  PropulsionSystemType getPropulsionSystemType(void) const;

  bool isSolarSail(void) const;
  bool isNEP(void) const;
  bool isNSTAR(void) const;
  bool isT6(void) const;
  bool isRIT(void) const;
  bool isSEP(void) const;
  bool isEP(void) const;

  /*! \brief type of coordinates used in state vector of equations of motion
   * \return used EOM coordinate type (e.g. POLAR, CARTESIAN)*/
  SpaceLib::CoordType getEomCoordinateType(void) const;

  /*! \brief set new type of coordinates to use in equations of motion
   * \param coord_type the new coordinate type
   * \sa getEomCoordinateType*/
  void setEomCoordinateType(const SpaceLib::CoordType coord_type);
      
             /// Minimum initial propellant mass
  private:   double m_min_initialpropellantmass;
             /// Get minimum initial propellant mass
  public:   double getMinInitialPropellantMass(void) const { return m_min_initialpropellantmass; };

             /// Maximum initial propellant mass
  private:   double m_max_initialpropellantmass;
             /// get maximum initial propellant mass 
  public:   double getMaxInitialPropellantMass(void) const { return m_max_initialpropellantmass; };

             /// Actual initial propellant mass
  private:   double m_initialpropellantmass;
             /// Get actual initial propellant mass           
  public:    double getInitialPropellantMass(void) const { return m_initialpropellantmass; };

             /// Actual propellant mass
  private:   double m_propellantmass;
             /// Get actual propellant mass
  public:    double getPropellantMass(void) const { return m_propellantmass; };
             /// Get actual propellant mass
  private:   void setPropellantMass(const double mass) { m_propellantmass = mass; };

             /// No more propellant in tank
  private:   bool m_empty;
  public:    void setEmpty(const bool yes_or_no) { m_empty = yes_or_no; };
  public:    bool isEmpty(void) const { return m_empty; };

  virtual void printVrmlAppearance(std::ostream& os,const std::string& sc_name) 
    const throw () {};
      
protected:

  bool operator==(const CPropSys& rhs)
  {
    return ( 
      (m_state                     == rhs.m_state) &&
      (m_type                      == rhs.m_type) &&
      (getEomCoordinateType() == rhs.getEomCoordinateType()) &&
      (m_min_initialpropellantmass == rhs.m_min_initialpropellantmass) &&
      (m_max_initialpropellantmass == rhs.m_max_initialpropellantmass) &&
      (m_initialpropellantmass     == rhs.m_initialpropellantmass) &&
      (m_propellantmass            == rhs.m_propellantmass) &&
      (m_empty                     == rhs.m_empty) 
     );
  };

public:

  /// Additional copy constructor
  CPropSys(const CPropSys& ps, CState* state = 0)
  {
    operator=(ps);
    m_state = state;
  }

  /// PropulsionSystemType constructor with given host object
  /// , i.e. a CSC object which is also a CState object
  CPropSys(const PropulsionSystemType psystype, CState* state = 0)
  {
    m_state                     = state;
    m_type                      = psystype;
    setEomCoordinateType(DFLT_EOM_COORD_TYPE);
    m_min_initialpropellantmass = 0.0;
    m_max_initialpropellantmass = 0.0;
    m_initialpropellantmass     = 0.0;
    m_propellantmass            = 0.0;
    m_empty                     = true;
  }

  /// Assign operator
  const CPropSys & operator=(const CPropSys & ps)
  {
    m_state                     = ps.m_state;
    m_type                      = ps.m_type;
    setEomCoordinateType(ps.getEomCoordinateType());
    m_min_initialpropellantmass = ps.m_min_initialpropellantmass;
    m_max_initialpropellantmass = ps.m_max_initialpropellantmass;
    m_initialpropellantmass     = ps.m_initialpropellantmass;
    m_propellantmass            = ps.m_propellantmass;
    m_empty                     = ps.m_empty;
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
   \return CPropSys Pointer to new CPropSys object
  */
  virtual CPropSys* clone(CState* state)
  {
    try 
    {
      return (new CPropSys(*this, state));
    } 
    catch (std::bad_alloc ) 
    {
      return 0;
    }
  };

  /// destructor
  virtual ~CPropSys(void) {};

  /// set minimum and maximum propellant masses
  void setInitialPropellantMassInterval(
     const double min_mass
    ,const double max_mass
  );


  /// tbd
  void setInitialPropellantMass(const double fraction);


  /// tbd
  void Refuel(void);

  // ----- virtual functions -----

  /*!\brief pointer to translational equations of motion (TEOM)
   *
   * The type of EOM depends on the type of coodinates being set internally.
   * You may alter this coordinate type with setEomCoordinateType.
   * \return function pointer to EOM function specified by EOM_RK_func.
   * \sa getEomCoordinateType, setEomCoordinateType*/
  virtual EOM_RK_func getTEOM(void);


  /*!\brief pointer to rotation equations of motion (REOM)
   * \param function pointer
   */
  virtual EOM_RK_func getREOM(const RDES_Type type) 
  {
    switch (type)
    {
    case DES_Q_QD: return REOM_RK_Quaternion;break;
    case DES_Q_OM: return REOM_RK_Quaternion_AngVel;break;
    default: FatalError("virtual CPropSys::EOM_RK_func getREOM(const RDES_Type type): don't know which REOM type!\n");
      return REOM_RK_Quaternion;
      break;
    }
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual void EOM_RK_POLAR(double t, const MathLib::Vector& y
    , MathLib::Vector& yp, void* pAux)
  {
    FatalError("class CPropSys: Function EOM_RK() not defined for the selected propulsion system");
  };

  /// Virtual function that gives error if not implemented in derived class
  virtual void EOM_RK_CARTESIAN(double t, const Vector& y, Vector& yp, void* pAux)
  {
    FatalError("class CPropSys: Function EOM_RK_CARTESIAN() not defined for the selected propulsion system");
  };

  /// Virtual function that gives error if not implemented in derived class
  virtual void setLightness(const double lambda)
  {
    FatalError("class CPropSys: Function setLightness() not defined for the selected propulsion system");
  };

  /// Virtual function that gives error if not implemented in derived class
  virtual void setOpticalModel(const bool ideal_or_real,
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
                               const double Bb = 0.0)
  {
    FatalError("class CPropSys: Function setOpticalModel() not defined for the selected propulsion system");
  };

  /// Virtual function that gives error if not implemented in derived class
  virtual void Init(const double srd_or_mp)
  {
    FatalError("class CPropSys: Function Init() not defined for the selected propulsion system");
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getLightness(void) const
  {
    FatalError("class CPropSys: Function getLightness() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual bool isIdealSail(void) const
  {
    FatalError("class CPropSys: Function isIdealSail() not defined for the selected propulsion system");
    return false; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual bool isRealSail(void) const
  {
    FatalError("class CPropSys: Function isRealSail() not defined for the selected propulsion system");
    return false; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getMaxTemperature(void) const
  {
    FatalError("class CPropSys: Function getMaxTemperature() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getrho(void) const
  {
    FatalError("class CPropSys: Function getrho() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double gets(void) const
  {
    FatalError("class CPropSys: Function gets() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getepsilonf(void) const
  {
    FatalError("class CPropSys: Function getepsilonf() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getepsilonb(void) const
  {
    FatalError("class CPropSys: Function getepsilonb() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getBf(void) const
  {
    FatalError("class CPropSys: Function getBf() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getBb(void) const
  {
    FatalError("class CPropSys: Function getBb() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getG(void) const
  {
    FatalError("class CPropSys: Function getG() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getH(void) const
  {
    FatalError("class CPropSys: Function getH() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getK(void) const
  {
    FatalError("class CPropSys: Function getK() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getG0(void) const
  {
    FatalError("class CPropSys: Function getG0() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getH0(void) const
  {
    FatalError("class CPropSys: Function getH0() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getK0(void) const
  {
    FatalError("class CPropSys: Function getK0() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getTemperature(const double r, const double coneangle) const
  {
    FatalError("class CPropSys: Function getTemperature() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getSolarRadiationDose(void) const
  {
    FatalError("class CPropSys: Function getSolarRadiationDose() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual void setSolarRadiationDose(const double srd)
  {
    FatalError("class CPropSys: Function setSolarRadiationDose() not defined for the selected propulsion system");
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getSolarEnergyFlux(const double r, const double coneangle) const
  {
    FatalError("class CPropSys: Function getSolarEnergyFlux() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual void setThrottle(const double throttle)
  {
    FatalError("class CPropSys: Function setThrottle() not defined for the selected propulsion system");
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual void setThrottleType(const ThrottleType throttletype)
  {
    FatalError("class CPropSys: Function setThrottleType() not defined for the selected propulsion system");
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getThrottle(void) const
  {
    FatalError("class CPropSys: Function getThrottle() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual ThrottleType getThrottleType(void) const
  {
    FatalError("class CPropSys: Function getThrottleType() not defined for the selected propulsion system");
    return UNDEFINED_TT; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual void CheckRemainingPropellant(const double dt)
  {
    FatalError("class CPropSys: Function CheckRemainingPropellant() not defined for the selected propulsion system");
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getIsp(void)
  {
    FatalError("class CPropSys: Function getIsp() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual void setMaxThrust(const double thrust_max)
  {
    FatalError("class CPropSys: Function setMaxThrust() not defined for the selected propulsion system");
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual void setIsp(const double isp)
  {
    FatalError("class CPropSys: Function setIsp() not defined for the selected propulsion system");
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual void setCharacteristicPower(const double characteristic_power)
  {
    FatalError("class CPropSys: Function setCharacteristicPower()not defined for the selected propulsion system");
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual void setPowerExponent(const double p)
  {
    FatalError("class CPropSys: Function setPowerExponent() not defined for the selected propulsion system");
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual void setNumberOfThrusters(const double thrusters)
  {
    FatalError("class CPropSys: Function setNumberOfThrusters() not defined for the selected propulsion system");
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getNumberOfThrusters(void) const
  {
    FatalError("class CPropSys: Function getNumberOfThrusters() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getNumberOfUsedThrusters(void) const
  {
    FatalError("class CPropSys: Function getNumberOfUsedThrusters() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getPowerOfUsedThruster(void) const
  {
    FatalError("class CPropSys: Function getPowerOfUsedThruster() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getPower(void) const
  {
    FatalError("class CPropSys: Function getPower() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getMaxPower(void) const
  {
    FatalError("class CPropSys: Function getMaxPower() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getThrust(void)
  {
    FatalError("class CPropSys: Function getThrust() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getMaxThrust(void)
  {
    FatalError("class CPropSys: Function getMaxThrust() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

  /// Virtual function that gives error if not implemented in derived class
  virtual double getMassFlow(void)
  {
    FatalError("class CPropSys: Function getMassFlow() not defined for the selected propulsion system");
    return 0.0; // dummy
  }

}; // end class CPropSys


//------------------------------------------------------------------------------
inline void EulerEquations(const MathLib::Vector& om,MathLib::Vector& om_dot
 ,const MathLib::Vector& moi, const MathLib::Vector& sm)
{
  om_dot(0) = 1/moi(0) * ( om(1)*om(2) * ( moi(1) - moi(2) ) + sm(0) );
  om_dot(1) = 1/moi(1) * ( om(2)*om(0) * ( moi(2) - moi(0) ) + sm(1) );
  om_dot(2) = 1/moi(2) * ( om(0)*om(1) * ( moi(0) - moi(1) ) + sm(2) );
}


//------------------------------------------------------------------------------
inline PropulsionSystemType CPropSys::getPropulsionSystemType(void) const 
{ 
  return m_type;
} // PropulsionSystemType getPropulsionSystemType(void) const

//------------------------------------------------------------------------------
inline bool CPropSys::isSolarSail(void) const
{ 
  return (getPropulsionSystemType() == SOLARSAIL);
} // bool isSolarSail(void) const

//------------------------------------------------------------------------------
inline bool CPropSys::isNEP(void) const 
{ 
  return (getPropulsionSystemType() == NEP);
} // bool isNEP(void) const 

//------------------------------------------------------------------------------
inline bool CPropSys::isNSTAR(void) const
{ 
  return (getPropulsionSystemType() == NSTAR);
} // bool isNSTAR(void) const

//------------------------------------------------------------------------------
inline bool CPropSys::isT6(void) const
{ 
  return (getPropulsionSystemType() == T6);
} // bool isT6(void) const

//------------------------------------------------------------------------------
inline bool CPropSys::isRIT(void) const
{ 
  return (getPropulsionSystemType() == RIT);
} // bool isRIT(void) const

//------------------------------------------------------------------------------
inline bool CPropSys::isSEP(void) const
{ 
  return (isNSTAR() || isT6() || isRIT());
} // bool isSEP(void) const

//------------------------------------------------------------------------------
inline bool CPropSys::isEP(void) const
{ 
  return (isNEP() || isSEP());
} // bool isEP(void) const

//------------------------------------------------------------------------------
inline SpaceLib::CoordType CPropSys::getEomCoordinateType(void) const
{
  return m_eom_type;
} // SpaceLib::CoordType getEomCoordinateType(void) const;

//------------------------------------------------------------------------------
inline void CPropSys::setEomCoordinateType(const SpaceLib::CoordType coord_type)
{
  m_eom_type = coord_type;
} // inline void setEomCoordinateType(const SpaceLib::CoordType coor_type)

//------------------------------------------------------------------------------
inline EOM_RK_func CPropSys::getTEOM(void)
{
  const std::string method_name = "virtual EOM_RK_func CPropSys::getTEOM(void)";
  const std::string err_msg_01 = ": unsupported type of EOM coordinates!\n";

  switch(getEomCoordinateType())
  {
  case SpaceLib::POLAR:
    return TEOM_RK_POLAR;
    break;
  case SpaceLib::CARTESIAN:
    return TEOM_RK_CARTESIAN;
    break;
  default:
    FatalError(method_name+err_msg_01);
    return 0;
    break;
  }
}

#endif // CPROPSYS_H_INCLUDED
