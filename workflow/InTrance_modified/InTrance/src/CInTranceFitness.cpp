/*!
\file CInTranceFitness.cpp
\version 1.7
\brief CInTranceFitness class implementation.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 14.05.2006 (AO)

\b Purpose:
   This file provides the interface for the class CInTranceFitness.
   CInTranceFitness implements the multi-objective fitness function of InTrance.

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

// InTrance
#include <configInTrance.h>

#include <math.h>

// SpaceLib
#include <SpaceLib.h>

// InTrance
#include <CInTranceFitness.h>
#include <InTranceConst.h>
                                
using namespace std;
using namespace MathLib;
using namespace Global;

TrajectoryType       tt;   // trajectory type must be global
OptimizationGoalType goal; // optimization goal must be global
FitnessFunctionType  jt;   // fitness function type must be global
int jt2;                   // secondary fitnesss function type
int jt3;                   // third fitnesss function type

const double CRIT_RENDEZVOUS_FITNESS[N_FFT] = 
{                                         
   1.0/sqrt(2.0) // J_AND
  ,1.0/sqrt(2.0) // J_OR
  ,1.0/sqrt(3.0) // J_H_AND
  ,1.0/sqrt(3.0) // J_H_OR
  ,1.0/sqrt(4.0) // J_EH_AND
  ,1.0/sqrt(4.0) // J_EH_OR
  ,1.0/sqrt(5.0) // J_AEH_AND
  ,1.0/sqrt(5.0) // J_AEH_OR
  ,0.0           // J_RIMPACT
  ,0.0           // J_RIMPACT_AND
  ,0.0           // J_RIMPACT_OR
  ,0.0           // J_ACT
  ,0.0           // J_SPI
};


double GetCriticalFitnessValue(void)
{
  switch (tt)
  {
  case RENDEZVOUS:
    return CRIT_RENDEZVOUS_FITNESS[jt];
    break;

  case FLYBY:
    return 1.0;
    break;

  case RIMPACT:
    switch (goal)
    {
    case RETROGRADE_IMPACT_CRANKING:
      return 99.0;
      break;
    case MAXIMUM_DEFLECTION:
      return 40.0;
      break;
    default:
      return 0.0;
      break;
    }
    break; // RIMPACT

  case ACTPROBLEM:
    return 1000.0; // to be changed
    break;

  case SOLAR_POLAR_IMAGER:
    return 1000.0; // to be changed
    break;

  default:
    return 0.0;
    break;
  }
}


//-----------------------------------------------------------------------------
double CInTranceFitness::getMixedFitness(void)
{
  switch (tt)
  {
    case GTOC2TRAJ: return m_fitness[MIXED_FITNESS]; break;
    default:break;
  } // swtich (tt)

  if (tt==RENDEZVOUS && jt==J_ORB_ELEMENTS) return m_fitness[MIXED_FITNESS];

  if ( (tt==RIMPACT) && (goal == RETROGRADE_IMPACT_CRANKING) )
    return (100.0 + m_fitness(RIMPH_FITNESS));

  double base_fitness;

  bool   pos_ok = (m_fitness(  POS_FITNESS) >=  0.0) ? true : false;
  bool   vel_ok = (m_fitness(  VEL_FITNESS) >=  0.0) ? true : false;
  bool rimpe_ok = (m_fitness(RIMPE_FITNESS) >= -3.0) ? true : false;
  bool rimph_ok = (m_fitness(RIMPH_FITNESS) >= -1.0) ? true : false;

  bool pos_and_vel_ok             =   pos_ok &&   vel_ok;
  bool rimpe_and_rimph_fitness_ok = rimpe_ok && rimph_ok;
  
  if (goal == TRANSFER_TIME)
    base_fitness = pos_and_vel_ok ? m_fitness(TIME_FITNESS) : 0.0;
  
  else if (goal == PROPELLANT_MASS)
    base_fitness = pos_and_vel_ok ? m_fitness(MPROP_FITNESS) : 0.0;
  
  else if (goal == MAXIMUM_DEFLECTION)
  {
    base_fitness = 10.0 + m_fitness(RIMPH_FITNESS) + m_fitness(RIMPE_FITNESS); // optimize rimpe and rimph
    if (rimpe_and_rimph_fitness_ok)
    {    
      base_fitness = 30.0 + m_fitness(POS_FITNESS); // optimize pos
      if (pos_ok)
      {
        base_fitness=50.0 + m_fitness(DEFL_FITNESS); // optimize deflection
      //base_fitness=50.0 + m_fitness(RIMPV_FITNESS); // optimize retrograde impact velocity
      } 
    }
  }

  else if (goal == SOLAR_POLAR_IMAGER)
  {
    base_fitness = m_fitness(SPI_I_FITNESS); // optimize semimajor axis
    if (m_fitness(SPI_I_FITNESS) >= 0.0)
    {    
      base_fitness = 1.0 + m_fitness(SPI_A_FITNESS); // optimize inclination
      if (m_fitness(SPI_A_FITNESS) >= 0.0)
      {
        base_fitness = 5.0 + m_fitness(SPI_E_FITNESS); // optimize eccentricity
        if (m_fitness(SPI_E_FITNESS) >= 0.0)
        {
          base_fitness = 10.0 + 100.0 * m_fitness(TIME_FITNESS); // optimize time
        } 
      } 
    }
  }
  else if (goal == ACT_DEFLECTION/* && m_fitness(TIME_FITNESS) > -1.0e-20*/)
  {
    base_fitness = m_fitness[MIXED_FITNESS];
  } else base_fitness = -20000;                                              
  
  double r = 1.0/(m_fitness(POS_FITNESS) < -10.0 ? 1.0E-10 : pow(10.0,m_fitness(POS_FITNESS)));
  if (r < 1.0) r = 0.99 + r/100.0;

  double v = 1.0/(m_fitness(VEL_FITNESS) < -10.0 ? 1.0E-10 : pow(10.0,m_fitness(VEL_FITNESS)));
  if (v < 1.0) v = 0.99 + v/100.0;

  double a = 1.0/(m_fitness(A_FITNESS) < -10.0 ? 1.0E-10 : pow(10.0,m_fitness(A_FITNESS)));
  if (a < 1.0) a = 0.99 + a/100.0;

  double e = 1.0/(m_fitness(E_FITNESS) < -10.0 ? 1.0E-10 : pow(10.0,m_fitness(E_FITNESS)));
  if (e < 1.0) e = 0.99 + e/100.0;

  double h = 1.0/(m_fitness(H_FITNESS) < -10.0 ? 1.0E-10 : pow(10.0,m_fitness(H_FITNESS)));
  if (h < 1.0) h = 0.99 + h/100.0;

  if (tt == RENDEZVOUS) // rendezvous
  {
    const double offset2 = 1.0/sqrt(2.0);
    const double offset3 = 1.0/sqrt(3.0);
    const double offset4 = 1.0/sqrt(4.0);
    const double offset5 = 1.0/sqrt(5.0);
    if (pos_and_vel_ok)
    {
      double m = 0.0;     
      switch (jt)
      {
      case J_AND:
      case J_OR:
        m = std::max(r,v);
        //  return offset2 + 10.0*base_fitness + 1.0/sqrt(2.0*m*m); // 2.0 due to m*m+m*m          
        return 10.0 * base_fitness + 1.0/sqrt(2.0*m*m); // 2.0 due to m*m+m*m
        break;

      case J_H_AND:
      case J_H_OR:
        m = std::max(r,v); // h does not play a role
        return 10.0 * base_fitness + 1.0/sqrt(3.0*m*m); // 3.0 due to m*m+m*m+m*m, although h is droped
        break;

      case J_EH_AND:
      case J_EH_OR:
        m = std::max(r,v); // e and h do not play a role
        return 10.0 * base_fitness + 1.0/sqrt(4.0*m*m); // 4.0 due to m*m+m*m+m*m+m*m, although e and h are droped
        break;

      case J_AEH_AND:
      case J_AEH_OR:
        m = std::max(r,v); // a, e, and h do not play a role
        return 10.0 * base_fitness + 1.0/sqrt(5.0*m*m); // 5.0 due to m*m+m*m+m*m+m*m+m*m, although a, e and h are droped
        break;

      case J_ORB_ELEMENTS:
        return m_fitness[MIXED_FITNESS];
        break;

      default:
        FatalError("CEAFitness:getMixedFitness(): unknown fitness function type");
        return 0.0; // to prevent compiler warning
        break;
      } // switch (jt)
    } // if (pos_and_vel_ok)
    else //if (!pos_and_vel_ok)
    {
      double m=0.0, m1=0.0, m2=0.0, m3=0.0;
      switch (jt)
      {
      case J_AND:
        m = std::max(r,v);
        return 10.0 * base_fitness + 1.0/sqrt(2.0*m*m); // 2.0 due to m*m+m*m
        break;

      case J_H_AND:
        m1 = std::max(r,v);
        m = std::max(m1,h); // h plays a role
        return 10.0 * base_fitness + 1.0/sqrt(3.0*m*m); // 3.0 due to m*m+m*m+m*m, although h is droped
        break;

      case J_EH_AND:
        m1 = std::max(r,v);
        m2 = std::max(e,h); // e and h play a role
        m  = std::max(m1,m2);
        return 10.0 * base_fitness + 1.0/sqrt(4.0*m*m); // 4.0 due to m*m+m*m+m*m+m*m, although e and h are droped
        break;

      case J_AEH_AND:
        m1 = std::max(r,v);
        m2 = std::max(a,e); // a and e play a role
        m3 = std::max(m1,m2);
        m  = std::max(h,m3); // and h also
        return 10.0 * base_fitness + 1.0/sqrt(5.0*m*m); // 5.0 due to m*m+m*m+m*m+m*m+m*m, although a, e and h are droped
        break;

      case J_OR:
        //r = Global::max(1.0,r);v = Global::max(1.0,v);
        return 1.0/sqrt(r*r+v*v)-offset2;
        //return 10.0 * base_fitness + 1.0/sqrt(r*r+v*v);
        break;

      case J_H_OR:        
        r = std::max(1.0,r);v = std::max(1.0,v);h= std::max(1.0,h);
        return 10.0 * base_fitness + 1.0/sqrt(r*r+v*v+h*h) - offset3;
        //return 10.0 * base_fitness + 1.0/sqrt(r*r+v*v+h*h);
        break;

      case J_EH_OR:
        r = std::max(1.0,r);v = std::max(1.0,v);h= std::max(1.0,h);
        e = std::max(1.0,e);
        return 10.0 * base_fitness + 1.0/sqrt(r*r+v*v+e*e+h*h) - offset4;
        //return 10.0 * base_fitness + 1.0/sqrt(r*r+v*v+e*e+h*h);
        break;

      case J_AEH_OR:
        r = std::max(1.0,r);v = std::max(1.0,v);h= std::max(1.0,h);
        e = std::max(1.0,e);a = std::max(1.0,a);
        return 10.0 * base_fitness + 1.0/sqrt(r*r+v*v+a*a+e*e+h*h) - offset5;
        break;

      case J_ORB_ELEMENTS:
        return m_fitness[MIXED_FITNESS];
        break;
      default:
        FatalError("CEAFitness:getMixedFitness(): unknown fitness function type");
        return 0.0; // to prevent compiler warning
        break;
      } // switch (jt)
    } // else //if (!pos_and_vel_ok)
  } // if (tt == RENDEZVOUS)
  
  else if (tt == FLYBY) // flyby
  {
    return 10.0 * base_fitness + 1.0/r;
  }
  // retrograde impact or ACT problem or SPI problem
  else if (tt == RIMPACT || tt == ACTPROBLEM || tt == SPIPROBLEM)
  {
    return base_fitness;
  }
  
  else
  {
    FatalError("CEAFitness::getMixedFitness(): unknown trajectory type");
    return -HUGE_VAL; // dummy to prevent compiler warning
  }

}; // end CInTranceFitness::getMixedFitness

//------------------------------------------------------------------------------
CEAFitness* CInTranceFitness::Duplicate(void) const 
{
  CInTranceFitness* fitness = 0;
  fitness = new CInTranceFitness;

#ifndef no_check_of_memory_allocation
  if (fitness == 0) FatalError("CInTranceFitness::Duplicate(void) => Error during memory allocating!");
#endif

  fitness->m_fitness = m_fitness;

  return fitness;

} // end CInTranceFitness::Duplicate(void)
