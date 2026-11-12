/*!
\file CSCSteering.h
\brief Two classes that provide the functionality for steering a spacecraft.

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

  09.05.2006

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
#ifndef CSPACECRAFTSTEERING_H_INCLUDED
#define CSPACECRAFTSTEERING_H_INCLUDED


// C/C++ math library
#include <math.h>

// MathLib
#include <MLGlobal.h>

// SpaceLib
#include <CState.h>
#include <CSC.h>


class CNC; // forward declaration

static const int N_LSL  = 15; // number of local steering laws
                              // 5 x minimize elements
                              // 5 x maximize elements
                              // 5 x optimize elements

enum LSLIndex { amax     =  0,
                amin     =  1,
                emax     =  2,
                emin     =  3,
                imax     =  4,
                imin     =  5,
                Omegamax =  6,
                Omegamin =  7,
                omegamax =  8,
                omegamin =  9,
                aopt     = 10,
                eopt     = 11,
                iopt     = 12,
                Omegaopt = 13,
                omegaopt = 14 };

static const int N_STST = 4;   // number of steering strategies

enum StStIndex { DEOPT     = 0,
                 DEMINMAX  = 1,
                 DIRECTION = 2,
                 ALPHABETA = 3 };


/// A class that provides the functionality for a local steering law (LSL) vector.
class CSCLSLVec : public Vector
// ----------------------------------------------------------------------------
// i-th vector component: factor for i-th local steering law
// ----------------------------------------------------------------------------
{
  friend class CSCSteering;

  // __constructors__

  // null vector
  CSCLSLVec()                : Vector(N_LSL) {}

  // Vector copy
  CSCLSLVec(const Vector& V) : Vector(V) {

    if(V.getDimension() != N_LSL) FatalError("CSCLSLVec(): source vector has wrong dimension");

  }

}; // end class CSCLSLVec


/// A class that provides the functionality for steering a spacecraft.
// ----------------------------------------------------------------------------
class CSCSteering
// ----------------------------------------------------------------------------
{
  friend class CSC;

//!!!
  friend class CSimulation;
//!!!

  //-----------------------------------------------------------------

  bool                  m_ga_steering;
                        // gravity assist optimization is included

  // local steering laws (LSL)
  int                   m_number_of_lsl;
                        // number of available local steering laws
  CSCLSLVec             m_lslvec;
                        // actual local steering law vector

  // steering strategies (StSt)
  int                   m_number_of_stst;
                        // number of available steering strategies
  StStIndex             m_stst_index;
                        // index of actual steering strategy

  // optimal thrust vector
  double                m_xi_r;
                        // radial component of the optimal thrust vector
  double                m_xi_t;
                        // transversal component of the optimal thrust vector
  double                m_xi_h;
                        // orbit normal component of the optimal thrust vector
  //-----------------------------------------------------------------

  //
  // local steering laws
  //

  void lsl_amax(const CState& ps, const CState& ts)
  // ----------------------------------------------
  {
    m_xi_r = ps.e()*sin(ps.f());
    m_xi_t = ps.SemilatusRectum()/ps.r();
    m_xi_h = 0.0;
  };

  void lsl_amin(const CState& ps, const CState& ts)
  // ----------------------------------------------
  {
    m_xi_r = -ps.e()*sin(ps.f());
    m_xi_t = -ps.SemilatusRectum()/ps.r();
    m_xi_h = 0.0;
  };

  void lsl_emax(const CState& ps, const CState& ts)
  // ----------------------------------------------
  {
    const double p = ps.SemilatusRectum();
    const double r = ps.r();
    const double f = ps.f();

    m_xi_r = p*sin(f);
    m_xi_t = (p+r)*cos(f)+r*ps.e();
    m_xi_h = 0.0;
  };

  void lsl_emin(const CState& ps, const CState& ts)
  // ----------------------------------------------
  {
    const double p = ps.SemilatusRectum();
    const double r = ps.r();
    const double f = ps.f();

    m_xi_r = -p*sin(f);
    m_xi_t = -((p+r)*cos(f)+r*ps.e());
    m_xi_h = 0.0;
  };

  void lsl_imax(const CState& ps, const CState& ts)
  // ----------------------------------------------
  {
    m_xi_r = 0.0;
    m_xi_t = 0.0;
    m_xi_h = cos(ps.f()+ps.omega());
  };

  void lsl_imin(const CState& ps, const CState& ts)
  // ----------------------------------------------
  {
    m_xi_r = 0.0;
    m_xi_t = 0.0;
    m_xi_h = -cos(ps.f()+ps.omega());
  };

  void lsl_Omegamax(const CState& ps, const CState& ts)
  // --------------------------------------------------
  {
    const double i = ps.i();

    m_xi_r = 0.0;
    m_xi_t = 0.0;

    if (i != 0.0) // to avoid singularity
    {
      m_xi_h = sin(ps.f()+ps.omega())/sin(i);
    }
    else // singular case
    {
      m_xi_h = 0.0;
    }
  };

  void lsl_Omegamin(const CState& ps, const CState& ts)
  // --------------------------------------------------
  {
    const double i = ps.i();

    m_xi_r = 0.0;
    m_xi_t = 0.0;

    if (i != 0.0) // to avoid singularity
    {
      m_xi_h = -sin(ps.f()+ps.omega())/sin(i);
    }
    else // singular case
    {
      m_xi_h = 0.0;
    }
  };

  void lsl_omegamax(const CState& ps, const CState& ts)
  // --------------------------------------------------
  {
    const double p = ps.SemilatusRectum();
    const double r = ps.r();
    const double e = ps.e();
    const double i = ps.i();
    const double f = ps.f();

    if (e != 0.0) // to avoid singularity
    {
      m_xi_r = -p*cos(f)/e;
      m_xi_t = (p+r)*sin(f)/e;
    }
    else // singular case
    {
      m_xi_r = 0.0;
      m_xi_t = 0.0;
    }

    if (i != 0.0) // to avoid singularity
    {
      m_xi_h = -r*sin(ps.omega()+f)*cos(i)/sin(i);
    }
    else // singular case
    {
      m_xi_h = 0.0;
    }
  };

  void lsl_omegamin(const CState& ps, const CState& ts)
  // --------------------------------------------------
  {
    const double p = ps.SemilatusRectum();
    const double r = ps.r();
    const double e = ps.e();
    const double i = ps.i();
    const double f = ps.f();

    if (e != 0.0) // to avoid singularity
    {
      m_xi_r = p*cos(f)/e;
      m_xi_t = -(p+r)*sin(f)/e;
    }
    else // singular case
    {
      m_xi_r = 0.0;
      m_xi_t = 0.0;
    }

    if (i != 0.0) // to avoid singularity
    {
      m_xi_h = r*sin(ps.omega()+f)*cos(i)/sin(i);
    }
    else // singular case
    {
      m_xi_h = 0.0;
    }
  };

  void lsl_aopt(const CState& ps, const CState& ts)
  // ----------------------------------------------
  {
    if (ts.a() >= ps.a())
      lsl_amax(ps, ts);
    else
      lsl_amin(ps, ts);
  };

  void lsl_eopt(const CState& ps, const CState& ts)
  // ----------------------------------------------
  {
    if (ts.e() >= ps.e())
      lsl_emax(ps, ts);
    else
      lsl_emin(ps, ts);
  };

  void lsl_iopt(const CState& ps, const CState& ts)
  // ----------------------------------------------
  {
    if (ts.i() >= ps.i())
      lsl_imax(ps, ts);
    else
      lsl_imin(ps, ts);
  };

  void lsl_Omegaopt(const CState& ps, const CState& ts)
  // --------------------------------------------------
  {
    double Omega_diff = DifferenceAngle(ts.Omega(), ps.Omega());

    if (Omega_diff >= 0.0)
      lsl_Omegamax(ps, ts);
    else
      lsl_Omegamin(ps, ts);
  };

  void lsl_omegaopt(const CState& ps, const CState& ts)
  // --------------------------------------------------
  {
    double omega_diff = DifferenceAngle(ts.omega(), ps.omega());

    if (omega_diff >= 0.0)
      lsl_omegamax(ps, ts);
    else
      lsl_omegamin(ps, ts);
  };

  //
  // pointer array to local steering laws
  //

  void (CSCSteering::*lsl_ptr[N_LSL])(const CState& ps, const CState& ts);

  //
  // functions
  //

  // __constructor__

  CSCSteering()
  {
    m_ga_steering = false;

    m_number_of_lsl   = N_LSL;
    m_number_of_stst  = N_STST;

    lsl_ptr[amax]     = &CSCSteering::lsl_amax;
    lsl_ptr[amin]     = &CSCSteering::lsl_amin;
    lsl_ptr[emax]     = &CSCSteering::lsl_emax;
    lsl_ptr[emin]     = &CSCSteering::lsl_emin;
    lsl_ptr[imax]     = &CSCSteering::lsl_imax;
    lsl_ptr[imin]     = &CSCSteering::lsl_imin;
    lsl_ptr[Omegamax] = &CSCSteering::lsl_Omegamax;
    lsl_ptr[Omegamin] = &CSCSteering::lsl_Omegamin;
    lsl_ptr[omegamax] = &CSCSteering::lsl_omegamax;
    lsl_ptr[omegamin] = &CSCSteering::lsl_omegamin;
    lsl_ptr[aopt]     = &CSCSteering::lsl_aopt;
    lsl_ptr[eopt]     = &CSCSteering::lsl_eopt;
    lsl_ptr[iopt]     = &CSCSteering::lsl_iopt;
    lsl_ptr[Omegaopt] = &CSCSteering::lsl_Omegaopt;
    lsl_ptr[omegaopt] = &CSCSteering::lsl_omegaopt;
  }

  bool operator==(const CSCSteering& rhs)
  {
    return (
      (m_ga_steering == rhs.m_ga_steering) &&
      (m_number_of_lsl == rhs.m_number_of_lsl) &&
      (m_number_of_stst == rhs.m_number_of_stst) &&
      (lsl_ptr[amax]     = rhs.lsl_ptr[amax]) &&
      (lsl_ptr[amin]     = rhs.lsl_ptr[amin]) &&
      (lsl_ptr[emax]     = rhs.lsl_ptr[emax]) &&
      (lsl_ptr[emin]     = rhs.lsl_ptr[emin]) &&
      (lsl_ptr[imax]     = rhs.lsl_ptr[imax]) &&
      (lsl_ptr[imin]     = rhs.lsl_ptr[imin]) &&
      (lsl_ptr[Omegamax] = rhs.lsl_ptr[Omegamax]) &&
      (lsl_ptr[Omegamin] = rhs.lsl_ptr[Omegamin]) &&
      (lsl_ptr[omegamax] = rhs.lsl_ptr[omegamax]) &&
      (lsl_ptr[omegamin] = rhs.lsl_ptr[omegamin]) &&
      (lsl_ptr[aopt]     = rhs.lsl_ptr[aopt]) &&
      (lsl_ptr[eopt]     = rhs.lsl_ptr[eopt]) &&
      (lsl_ptr[iopt]     = rhs.lsl_ptr[iopt]) &&
      (lsl_ptr[Omegaopt] = rhs.lsl_ptr[Omegaopt]) &&
      (lsl_ptr[omegaopt] = rhs.lsl_ptr[omegaopt])
    );

  };

  // __accessors__

  int getNumberOfLSL(void) const;

  int getNumberOfLStSt(void) const;

//  const char *getSteeringStrategyName(const StStIndex stst) const;

  StStIndex getSteeringStrategy(void) const;

  // __mutators__

  void setSteeringStrategy(const StStIndex stst);
                           // steering strategy index

  // purpose: set actual steering strategy

  // __class functionality__: perform spacecraft steering

  void setNCInputCartesian(int first_neuron,
                           int last_neuron,
                           const CState& ps,
                           const CState& ts,
                           CSC* sc,
                           const CState* ga = NULL);

  void setNCInputPolar(int first_neuron,
                       int last_neuron,
                       const CState& ps,
                       const CState& ts,
                       CSC* sc,
                       const CState* ga = NULL);

  void setNCInputKepler(int first_neuron,
                        int last_neuron,
                        const CState& ps,
                        const CState& ts,
                        CSC* sc,
                        const CState* ga = NULL);

  void setNCInputPropellant(int neuron, CSC* sc);

  void setNCInput(CNC*          nc,
                  const CState& ps,
                  const CState& ts,
                  CSC*          sc,
                  const CState* ga);

  void performSteering(CSC* sc,
                       // pointer to the spacecraft object
                       CNC* nc,
                       // pointer to the neurocontroller used for steering
                       const CState& ps,
                       // actual state of spacecraft (neurocontroller input)
                       const CState& ts,
                       // actual state of target body (neurocontroller input)
                       const CState* ga = 0);

  // purpose: do spacecraft steering with neurocontroller

  void performSteering(CSC* sc,
                       // pointer to the spacecraft object
                       const Vector& lsl_vector,
                       // local steering law vector
                       const CState& ps,
                       // actual state of spacecraft
                       const CState& ts);
                       // actual state of target body

  // purpose: do spacecraft steering with local steering law vector

  void performSteering(CSC* sc,
                       // pointer to the spacecraft object
                       const double clock_angle,
                       // spacecraft clock angle
                       const double cone_angle);
                       // spacecraft cone angle

  /*! \brief perfom a steering operation with given neurocontroller and
   *stored input neuron values.
   * \param sc neurocontroller to use for steering*/
  void performSteering(CSC* sc, CNC* nc);

  // purpose: do spacecraft steering with direct input of attitude angles

}; // end class CSCSteering

#endif // CSPACECRAFTSTEERING_H_INCLUDED