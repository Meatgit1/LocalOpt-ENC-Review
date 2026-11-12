/*!
\file CApproxFitness.h
\version 1.7
\brief CApproxFitness class interface.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 02.05.2005 (AO)

\b Purpose:


\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.

\b Copyright:
   (c) 2001-2005 Bernd Dachwald, Andreas Ohndorf
*/

#ifndef CAPPROXFITNESS_H_INCLUDED
#define CAPPROXFITNESS_H_INCLUDED

// STL
#include <cassert>

// InTrance
#include <InTranceConst.h>

// ENCLib
#include <CEAFitness.h>

// SpaceLib
#include <Global.h>

extern OptimizationGoalType goal; // global in CInTranceFitness.cpp
extern FitnessFunctionType  jt;   // global in CInTranceFitness.cpp
extern int jt2;                   // global in CInTranceFitness.cpp
extern int jt3;                   // global in CInTranceFitness.cpp

namespace approx
{
  /// application specific secondary fitness function type
  enum FitnessFunctionType2   /// fitness function type 2 (secondary)
  {
    NO_FFT1,          /// after having matched bound. cond. don't use fft1 anymore
    PLUS_FFT1,  /// continue assessing fft1 even after having matched bound. cond.

    FFT2_END
  };

  enum FitnessFunctionType3
  {
    SEGM_SIGNCHANGE1,   // use segmentation and sign changes of first derivative
    SIGNCHANGE,           // use only number of sign changes in first derivative
    SEGM_SIGNCHANGE2,   // use segmentation and sign changes of first derivative

    FFT3_END
  };

  enum ApproxFitnessIndex
  {
    CLOCK_SEG       //= 0
   ,CONE_SEG        //= 1
   ,SEGMENTS        //= 2
   ,CLOCK_SC        //= 3
   ,CONE_SC         //= 4
   ,SIGN_CHANGES    //= 5
   ,OPT_FITNESS     //= 8 (original CInTranceFitness object fitness)
   ,MAX_CLOCK       //= 9
   ,MIN_CLOCK       //= 10
   ,MAX_CONE        //= 11
   ,MIN_CONE        //= 12
   ,MIXED_FITNESS   //= 13
  };

  const int FITNESS_DIM(MIXED_FITNESS+1);


  /// Class that implements the multi-objective fitness function for trajectory
  /// approximation.
  class CApproxFitness: virtual public CEAFitness
  {

  public:

    /*!
     * \brief Constructor.
     */
    CApproxFitness(void) : CEAFitness(FITNESS_DIM) {};


    /*!
     * \brief Destructor.
     */
    virtual ~CApproxFitness() {};

    /*!
     * \brief Returns overall fitness.
     */
    virtual double getMixedFitness(void)
    {
      const double clock_seg_fit(m_fitness[CLOCK_SEG]);
      const double cone_seg_fit(m_fitness[CONE_SEG]);
      setFitness(SEGMENTS, clock_seg_fit*cone_seg_fit);

      const double clock_sc_fit(m_fitness[CLOCK_SC]);
      const double cone_sc_fit(m_fitness[CONE_SC]);
      setFitness(SIGN_CHANGES, 0.5*(clock_sc_fit+cone_sc_fit));

      const double sc_fitness(m_fitness[SIGN_CHANGES]);
      const double opt_fit(m_fitness[OPT_FITNESS]);
      const double macl_fit(m_fitness[MAX_CLOCK]*m_fitness[MAX_CLOCK]);
      const double micl_fit(m_fitness[MIN_CLOCK]*m_fitness[MIN_CLOCK]);
      const double maco_fit(m_fitness[MAX_CONE]*m_fitness[MAX_CONE]);
      const double mico_fit(m_fitness[MIN_CONE]*m_fitness[MIN_CONE]);

      // combined steering angles fitness (bound. cond. are matched
      // when this value reaches 1.0)
      const double ang_fit((macl_fit + micl_fit + maco_fit + mico_fit)/4.0);

      double cl1 = 1.0/(clock_seg_fit<-10.0 ? 1.0E-10 : pow(10.0,clock_seg_fit));
      // in case of having satisfied this boundary condition (a1==1.0)
      // modify subfitness value calculation
      if (cl1 < 1.0) cl1 = 0.99 + cl1/100.0;

      double co1 = 1.0/(cone_seg_fit<-10.0 ? 1.0E-10 : pow(10.0,cone_seg_fit));
      // in case of having satisfied this boundary condition (a1==1.0)
      // modify subfitness value calculation
      if (co1 < 1.0) co1 = 0.99 + co1/100.0;

      double cl2 = 1.0/(clock_sc_fit<-10.0 ? 1.0E-10 : pow(10.0,clock_sc_fit));
      // in case of having satisfied this boundary condition (a1==1.0)
      // modify subfitness value calculation
      if (cl2 < 1.0) cl2 = 0.99 + cl2/100.0;

      double co2 = 1.0/(cone_sc_fit<-10.0 ? 1.0E-10 : pow(10.0,cone_sc_fit));
      // in case of having satisfied this boundary condition (a1==1.0)
      // modify subfitness value calculation
      if (co2 < 1.0) co2 = 0.99 + co2/100.0;

      double crit_val(GetCriticalFitnessValue());
      double o( opt_fit/crit_val );

      double ret_val(0.0);
      if (ang_fit == 1.0)
      {
        if (o<1.0) // boundary conditions are not satisfied
        {
          switch (jt3)  // switch to determine composition of final fitness
          {             // depending on segmentation, sign changes of first derivative
                        // and possibly original fitness function value
          case SEGM_SIGNCHANGE1:               // just use orig. fitness value
            ret_val = o;
            break;                             
          case SIGNCHANGE:             // use sign changes, orig. fitness value
            ret_val = clock_sc_fit*cone_sc_fit*o;
              break;
          case SEGM_SIGNCHANGE2: // use segmentation, sign changes, orig. fitness value
            ret_val = sqrt(21.0/(5.0*cl1*cl1 + 5.0*co1*co1+ 5.0*cl2*cl2+5.0*co2*co2+1.0/(o*o)));
            break;
          default:
            ret_val = 0.0;
            break;
          }
        } 
        else // else branch of if (o<1.0) and therefore boundary conditions
        {    // are satisfied
          switch (jt2)  // switch to determine whether original fitness value
          {             // is to be assessed for final fitness value composition
                        // even after having satisfied boundary conditions
          case NO_FFT1:    // no further assessment of original fitness function
            switch (jt3)
            {
            case SEGM_SIGNCHANGE1:  // segmentation, sign changes
              ret_val = 1000.0*m_fitness[SEGMENTS]*m_fitness[SIGN_CHANGES];
              break;
            case SIGNCHANGE:                      // sign changes
              ret_val = 5.0*(clock_sc_fit*cone_sc_fit);
              break;
            case SEGM_SIGNCHANGE2:  // segmentation, sign changes
              ret_val = 1.0/(cl1*co1+cl2*co2);
              break;
            default:
              ret_val = 0.0;
              break;
            }
            break;


          case PLUS_FFT1: // original fitness function will still be assessed
            switch (jt3)
            {
            case SEGM_SIGNCHANGE1: // segmentation, sign changes, orig. fitness
              ret_val = 1000.0* m_fitness[SEGMENTS]*m_fitness[SIGN_CHANGES] 
                        + 0.01*opt_fit;
              break;
            case SIGNCHANGE:                     // sign changes, orig. fitness
              ret_val = 5.0*(clock_sc_fit*cone_sc_fit) + 0.01 * opt_fit;
              break;
            case SEGM_SIGNCHANGE2: // segmentation, sign changes, orig. fitness
              ret_val = 1.0/(cl1*co1+cl2*co2);
              break;
            default:
              ret_val = 0.0;
              break;
            }
          break;
          default:
            ret_val = 0;
            break;
          }
        }
      } else ret_val = (ang_fit-1);
      return ret_val;
                                              
    } // CApproxFitness::getMixedFitness(void) const
                                   
                  
    /*!   
     * \brief Returns pointer to a copy of itself.
     */                         
    virtual CEAFitness* Duplicate(void) const
    {
      CApproxFitness* f = 0;      // pointer initialization
      f = new CApproxFitness;     // new fitness object
      assert(f != 0);             // memory allocation assertation
      f->m_fitness = m_fitness;   // copy fitness data
      return f;
    }

  }; // end class CEAFitness

} // namespace approx
#endif
