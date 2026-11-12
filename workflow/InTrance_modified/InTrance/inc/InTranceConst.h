/*!
\file InTranceConst.h
\version 1.7
\brief InTrance constants and enums.
\author Bernd Dachwald <bernd.dachwald@dlr.de> (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 12.10.2006 (AO)

\b Purpose:
   This file provides the constants and enums for InTrance.

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

#ifndef INTRANCECONST_H_INCLUDED
#define INTRANCECONST_H_INCLUDED

#include <string>

const std::string COMMENT_SEPERATOR("//");  // C++ like double slash
const std::string TRAILING_CHARS(" \t");    // blanks and tabs

enum InTranceCommand     ///< available InTrance commands
{ 
  OPTIMIZE_NC       = 0, ///< optimize neurocontroller
  EVALUATE_NC       = 1, ///< evaluate neurocontroller
  EVALUATE_CONTROL  = 2, ///< evaluate control history
  LSL_MINMAX        = 3, ///< fly with minimize-maximize-steering strategy
  LSL_OPT           = 4, ///< fly with optimize-steering strategy
  ANGLES            = 5, ///< fly with constant steering angles
  OPTIMIZE_APPROX   = 6, ///< approximate given steering angles
  EVALUATE_APPROX   = 7, ///< evaluate steering angles approximation
  N_INTRANCE_COMMANDS
};

/// trajectory type
enum TrajectoryType
{
   RENDEZVOUS                                       ///< rendezvous problem type
  ,FLYBY                                                 ///< flyby problem type
  ,RIMPACT                        ///< retrograde impact trajectory problem type
  ,ACTPROBLEM        ///< 1st ESA ACT Global Trajectory Optimization Competition
  ,SPIPROBLEM                       ///< Solar Polar Imager mission problem type
  ,GTOC2TRAJ          ///< second Global Trajectory Optimization Problem type
  ,TRAJ_TEST1
  ,TRAJ_TEST2
  ,TRAJ_TEST3
  ,TRAJ_TEST4
  ,TRAJ_TEST5
  ,TRAJ_TEST6
  ,TRAJ_TEST7
  ,TRAJ_TEST8
  ,TRAJ_TEST9
  ,N_TRAJ_TYPES                                       ///< # of trajectory types
};

/// literal description of possible trajectory types
const std::string TRAJ_TYP_DESC[N_TRAJ_TYPES]=
{
   "rendezvous trajectory"
  ,"fly-by trajectory"
  ,"retrograde impact trajectory"
  ,"trajectory for 1st ESA ACT Global Trajectory Optimization Competition"
  ,"Solar Polar Imager mission problem trajectory"
  ,"2nd Global Trajectory Optimization Competition Problem"
  ,"Test Trajectory 1"
  ,"Test Trajectory 2"
  ,"Test Trajectory 3"
  ,"Test Trajectory 4"
  ,"Test Trajectory 5"
  ,"Test Trajectory 6"
  ,"Test Trajectory 7"
  ,"Test Trajectory 8"
  ,"Test Trajectory 9"
};

/// optimization goal type
enum OptimizationGoalType
{
   TRANSFER_TIME                                     ///< minimize transfer time
  ,PROPELLANT_MASS                          ///< minimize propellant consumption
  ,MAXIMUM_DEFLECTION 
                     ///< max. deflection w.r.t. formerly undisturbed trajectory
  ,RETROGRADE_IMPACT_CRANKING 
                         ///< max. orbit cranking for subsequent retrogr. impact
  ,ACT_DEFLECTION
             ///< for the 1st ESA ACT Global Trajectory Optimization Competition
  ,SOLAR_POLAR_IMAGER       ///< optimize solar polar imager transfer trajectory
  ,GTOCOMP2                    ///< Global Trajectory Optimization Competition 2
  ,OPT_GOAL_TEST1
  ,OPT_GOAL_TEST2
  ,OPT_GOAL_TEST3
  ,OPT_GOAL_TEST4
  ,OPT_GOAL_TEST5
  ,OPT_GOAL_TEST6
  ,OPT_GOAL_TEST7
  ,OPT_GOAL_TEST8
  ,OPT_GOAL_TEST9
  ,N_OPT_GOAL_TYPES                            ///< # of optimization goal types
};

/// optimization goal types (for program logic purposes)
const OptimizationGoalType OPT_GOAL_TYPES[N_OPT_GOAL_TYPES] =
{
   TRANSFER_TIME
  ,PROPELLANT_MASS
  ,MAXIMUM_DEFLECTION
  ,RETROGRADE_IMPACT_CRANKING
  ,ACT_DEFLECTION
  ,SOLAR_POLAR_IMAGER
  ,GTOCOMP2
  ,OPT_GOAL_TEST1
  ,OPT_GOAL_TEST2
  ,OPT_GOAL_TEST3
  ,OPT_GOAL_TEST4
  ,OPT_GOAL_TEST5
  ,OPT_GOAL_TEST6
  ,OPT_GOAL_TEST7
  ,OPT_GOAL_TEST8
  ,OPT_GOAL_TEST9
};

/// optimization goal parameter name
const std::string OPT_GOAL_TYPE_PARAM_VALUE[N_OPT_GOAL_TYPES] =
{
   "MINIMUMTRANSFERTIME"
  ,"MINIMUMPROPELLANTMASS"
  ,"DEFLECTION"
  ,"RETROGRADEIMPACTORBITCRANKING"
  ,"ACTDEFLECTION"
  ,"SOLARPOLARIMAGER"
  ,"GTOC2"
  ,"OPT_GOAL_TEST1"
  ,"OPT_GOAL_TEST2"
  ,"OPT_GOAL_TEST3"
  ,"OPT_GOAL_TEST4"
  ,"OPT_GOAL_TEST5"
  ,"OPT_GOAL_TEST6"
  ,"OPT_GOAL_TEST7"
  ,"OPT_GOAL_TEST8"
  ,"OPT_GOAL_TEST9"
};

//  text description of optimization goals
const std::string OPT_GOAL_DESC[N_OPT_GOAL_TYPES] =
{
   "minimum transfer time"
  ,"minimum consumed propellant mass"
  ,"maximum deflection"
  ,"maximum orbit cranking"
  ,"1st ESA ACT Global Trajectory Optimization Competition"
  ,"Solar Polar Imager Mission"
  ,"2nd Global Trajectory Optimization Competition"
  ,"not named yet"
  ,"not named yet"
  ,"not named yet"
  ,"not named yet"
  ,"not named yet"
  ,"not named yet"
  ,"not named yet"
  ,"not named yet"
  ,"not named yet"
};

// Application specific RDC constants.
const double LIMIT_HYPERSPACE_LARGE                 (1e-1); ///< magic
const double LIMIT_FINAL_BOUNDARY_CONDITION_MET     (6e-1); ///< magic
const double EPOCH_CONVERGENCE_LIMIT_FINAL_BOUNDARY_CONDITION_MET     (1e-5); ///< magic
const double EPOCH_CONVERGENCE_LIMIT_FINAL_BOUNDARY_CONDITION_NOT_MET (1e-3); ///< magic
const double EA_CONVERGENCE_LIMIT_FINAL_BOUNDARY_CONDITION_MET     (1e-6); ///< magic
const double EA_CONVERGENCE_LIMIT_FINAL_BOUNDARY_CONDITION_NOT_MET (1e-2); ///< magic

// application specific convergence criteria (steering angles approximation)
const double APPROX_HCS_UPPER_LIMIT(0.05);      // min. hypercube size upper limit
const double APPROX_FBC_MET(0.0);               // min. fitness value to satisfy final boundary conditions (FBC)
const double APPROX_POP_CONV_FBC_MET(1E-5);     // max. relative improvement of population convergence (FBC met)
const double APPROX_POP_CONV_FBC_NOT_MET(1E-4); // max. relative improvement of population convergence (FBC not met)
const double APPROX_EPC_CONV_FBC_MET(1E-4);     // max. relative improvement of epoche convergence (FBC met)
const double APPROX_EPC_CONV_FBC_NOT_MET(1E-6); // max. relative improvement of epoche convergence (FBC not met)

// Application specific fitness function types.
enum FitnessFunctionType /// fitness function type
{
  J_AND         =  0, ///< magic
  J_OR          =  1, ///< magic
  J_H_AND       =  2, ///< magic
  J_H_OR        =  3, ///< magic
  J_EH_AND      =  4, ///< magic
  J_EH_OR       =  5, ///< magic
  J_AEH_AND     =  6, ///< magic
  J_AEH_OR      =  7, ///< magic
  J_RIMPACT     =  8, ///< magic
  J_RIMPACT_AND =  9, ///< magic
  J_RIMPACT_OR  = 10, ///< magic
  J_ACT         = 11, ///< for the 1st ESA ACT Global Trajectory Optimization Competition
  J_SPI         = 12, ///< solar polar imager
  J_ORB_ELEMENTS= 13, ///< orbital elements matching
  J_GTOCOMP2    = 14,
  J_TEST_1      = 15,
  J_TEST_2      = 16,
  J_TEST_3      = 17,
  J_TEST_4      = 18,
  J_TEST_5      = 19,
  J_TEST_6      = 20,
  J_TEST_7      = 21,
  J_TEST_8      = 22,
  J_TEST_9      = 23,
  N_FFT         = 24  ///< number of fitness function types
};

/// fitness function type description
const std::string FitnessFunctionTypeDesc[] =
{
  "J_AND",             //  0
  "J_OR",              //  1
  "J_H_AND",           //  2
  "J_H_OR",            //  3
  "J_EH_AND",          //  4
  "J_EH_OR",           //  5
  "J_AEH_AND",         //  6
  "J_AEH_OR",          //  7
  "J_RIMPACT",         //  8
  "J_RIMPACT_AND",     //  9
  "J_RIMPACT_OR",      // 10
  "J_ACT",             // 11
  "J_SPI",             // 12
  "J_ORB_ELEMENTS",    // 13
  "J_GTOCOMP2",        // 14
  "J_TEST1",           // 15
  "J_TEST2",           // 16
  "J_TEST3",           // 17
  "J_TEST4",           // 18
  "J_TEST5",           // 19
  "J_TEST6",           // 20
  "J_TEST7",           // 21
  "J_TEST8",           // 22
  "J_TEST9"            // 23
};


// Application specific subfitness types.
enum SubfitnessType /// subfitness type \todo numbers should be changed somewhen
{
  MIXED_FITNESS  =  0, ///< final mixed fitness
  TIME_FITNESS   =  1, ///< transfer time
  POS_FITNESS    =  2, ///< final distance
  VEL_FITNESS    =  3, ///< final relative velocity
  A_FITNESS      =  4, ///< final difference in semi-major axis
  E_FITNESS      =  5, ///< final difference in eccentricity
  I_FITNESS      =  6, ///< final difference in inclination
  NODE_FITNESS   =  7, ///< final difference in longitude of ascending node
  PERI_FITNESS   =  8, ///< final difference in argument of pericenter
  M_FITNESS      =  9, ///< final difference in mean anomaly
  H_FITNESS      = 10, ///< final difference in orientation of the orbital plane
  IMPV_FITNESS   = 11, ///< impact velocity (including orientation)
  RIMPH_FITNESS  = 12, ///< final difference in orientation of the orbital plane (for retrograde impact)
  RIMPE_FITNESS  = 13, ///< final difference in orientation of perihelion (for retrograde impact)
  ACT_FITNESS    = 14, ///< fitness for the 1st ESA ACT Global Trajectory Optimization Competition
  MPROP_FITNESS  = 15, ///< propellant consumption
  GA_FITNESS     = 16, ///< gravity assist fitness
  DEFL_FITNESS   = 17, ///< deflection fitness
  SF_ACCURACY    = 18, ///< dynamic accuracy fitness
  SPI_A_FITNESS  = 19, ///< solar polar imager semi-major axis
  SPI_I_FITNESS  = 20, ///< solar polar imager inclination
  SPI_E_FITNESS  = 21, ///< solar polar imager eccentricity
  N_SUB_FITNESS  = 22  ///< total # of sub fitnesses
};

/// subfitness description text array
const std::string NAMES_SUBFITNESS[] =
{
  "final mixed fitness"                 //  0
 ,"t (time)"                            //  1
 ,"r (position)"                        //  2
 ,"v (velocity)"                        //  3
 ,"a (semimajor axis)"                  //  4
 ,"e (excentricity)"                    //  5
 ,"i (inclination)"                     //  6
 ,"O (node longitude)"                  //  7
 ,"o (pericenter argument)"             //  8
 ,"M (mean anomaly)"                    //  9
 ,"h (angular momentum)"                // 10
 ,"impact velocity"                     // 11
 ,"h_r (retro. impact)"                 // 12
 ,"e_r (retro. impact)"                 // 13
 ,"act (special)"                       // 14
 ,"m (mass. consumption)"               // 15
 ,"ga (gravity assist)"                 // 16
 ,"defl (deflection)"                   // 17
 ,"precision"                           // 18
 ,"Solar Polar Imager semimajor axis"   // 19
 ,"Solar Polar Imager inclination"      // 20
 ,"Solar Polar Imager eccentricity"     // 21
};

#endif // INTRANCECONST_H_INCLUDED
