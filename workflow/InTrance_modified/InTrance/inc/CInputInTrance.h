/*!
\file CInputInTrance.h
\version 1.7
\brief CInputInTrance class interface.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 21.09.2006 (AO)

\b Purpose:
   This file provides the interface for the class CInputInTrance.
   CInputInTrance is an auxiliary class that contains the values of the InTrance input file.

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
#ifndef CINPUTINTRANCE_H_INCLUDED
#define CINPUTINTRANCE_H_INCLUDED

//InTrance
#include <configInTrance.h>

// SpaceLib
#include <Global.h>

// InTrance
#include <InTranceConst.h>

const unsigned int DFLT_RANDOM_SEED = 1;

// possible InTrance commands
const InTranceCommand COMMANDS[] =
{
   OPTIMIZE_NC
  ,EVALUATE_NC
  ,EVALUATE_CONTROL
  ,LSL_MINMAX
  ,LSL_OPT
  ,ANGLES
  ,OPTIMIZE_APPROX
  ,EVALUATE_APPROX
};

// possible command parameter values associated to certain commands
const std::string COMMAND_PARAM_VALUE[] = 
{
   "OPTIMIZENC"
  ,"EVALUATENC"
  ,"EVALUATECONTROL"
  ,"LSLMINIMIZEMAXIMIZE"
  ,"LSLOPTIMIZE"
  ,"STEERINGANGLES"
  ,"OPTIMIZESTEERINGAPPROXIMATION"
  ,"EVALUATESTEERINGAPPROXIMATION"
};

/// Auxiliary class that contains the values of the InTrance input file.
class CInputInTrance
{
public:

  std::string path;                            ///< possible path of input files

  bool random_seed_provided;
  unsigned int random_seed; ///< seed for (pseudo-)random number generation

  InTranceCommand what_to_do; ///< what should the program do?

  // input files

  std::string simulation_parameter_file;             ///< name of the simulation parameter input file
  std::string spacecraft_parameter_file;             ///< name of the spacecraft parameter input file
  std::string nc_topology_file;                      ///< name of the neurocontroller topology file
  std::string snc_topology_file;                     ///< name of steering nc describing configuration file
  std::string evolutionary_algorithm_parameter_file; ///< name of the evolutionary algorithm parameter input file
  std::string evolutionary_algorithm_warmstart_file; ///< name of the evolutionary algorithm file containing the start solution
  std::string simulation_input_file_ang;             ///< name of the simulation input file containing the steering angles
  std::string clock_approx_data_file;                ///<
  std::string cone_approx_data_file;                 ///<
  std::string reference_trajectory_file;             ///<
  std::string rot_state_traj_file;                   ///<

  // input/output files

  std::string evolutionary_algorithm_best_solution_file;
        // name of the file containing the best solution found by
        // the evolutionary algorithm
        // =  input file, if what_to_do = EVALUATE_NC
        // = output file, if what_to_do = OPTIMIZE_NC

  // output files

  std::string simulation_output_file_csv;   ///< name of the simulation CSV output file
  std::string simulation_output_file_dat;   ///< name of the simulation DAT output file
  std::string simulation_output_file_gesop; ///< name of the GESOP output file
  std::string simulation_output_file_wrl;   ///< name of the simulation VRML output file
  std::string simulation_output_file_ang;   ///< name of the simulation output file containing the steering angles
  std::string report_file;                  ///< name of the report file

  // local steering law vector components

  double clock; ///< clock angle
  double cone;  ///< cone angle

  double     a_max; ///< component of LSL vector for the maximization of a
  double     a_min; ///< component of LSL vector for the minimization of a
  double     a_opt; ///< component of LSL vector for the optimization of a
  double     e_max; ///< component of LSL vector for the maximization of e
  double     e_min; ///< component of LSL vector for the minimization of e
  double     e_opt; ///< component of LSL vector for the optimization of e
  double     i_max; ///< component of LSL vector for the maximization of i
  double     i_min; ///< component of LSL vector for the minimization of i
  double     i_opt; ///< component of LSL vector for the optimization of i
  double Omega_max; ///< component of LSL vector for the maximization of Omega
  double Omega_min; ///< component of LSL vector for the minimization of Omega
  double Omega_opt; ///< component of LSL vector for the optimization of Omega
  double omega_max; ///< component of LSL vector for the maximization of omega
  double omega_min; ///< component of LSL vector for the minimization of omega
  double omega_opt; ///< component of LSL vector for the optimization of omega


  /*!
   * \brief Constructor.
   */
  CInputInTrance(void);

  /// destructor
  ~CInputInTrance() {};

  /// read parameters for local steering law minimization/maximization
  void readLslMinMaxParams(Global::ParamList&)
    throw (std::invalid_argument);

  /// read parameters for local steering law optimization
  void readLslOptParams(Global::ParamList&) throw (std::invalid_argument);

  /// read neurocontroller parameter input filename
  void readNcParamFilename(Global::ParamList&) throw (std::invalid_argument);

  /*!
   *\brief Reads parameters for approximation optimization.
   *\param params A list of parameter string pairs
   */
  void readOptimizeApproxParams(ParamList& parameters)
    throw (std::invalid_argument&);

  /// read parameters for "OPTIMIZE NC" & "EVALUATE NC" commands
  void readOptimizeNcParams(Global::ParamList&) throw (std::invalid_argument);

  /// read "ANGLES"-command parameters
  void readAngleParams(Global::ParamList&) throw (std::invalid_argument);

  /// read filename of best chromosome
  void readBestChromFilename(Global::ParamList&,
    bool file_existence_check = false) throw (std::invalid_argument);

  /*! \brief read and parse InTrance command */
  void readCommand(Global::ParamList&) throw (std::invalid_argument);

  /// read filename of desired csv-output file
  void readCsvFilename(Global::ParamList&) throw (std::invalid_argument&);

  /// read control values output filename
  void readCtrlFilename(Global::ParamList&) throw (std::invalid_argument);

  /// read evolutionary algorithm parameter filename
  void readEaParamFilename(Global::ParamList&) throw (std::invalid_argument);

  /*!
   *\brief Reads parameters for evaluation of approximation ncs.
   *\param params A list of parameter string pairs
   */
  void readEvaluateApproxOptimizationParams(
    ParamList& parameters) throw (std::invalid_argument&);

  /// read parameters need for evaluation of given control history
  void readEvaluateControlParams(Global::ParamList&)
    throw (std::invalid_argument&);

  /*! \brief input file identifier string and check for validity */
  void readFileIdentifier(Global::ParamList&) throw (std::invalid_argument&);

  /// read name of GESOP data output file
  void readGESOPFilename(Global::ParamList& p) throw (std::invalid_argument&);
            
  /// Read InTrance input file.
  /*! \param parameters tbd */
  void readParams(ParamList& parameters)
    throw (std::invalid_argument&);

  /// read input files pathname
  void readPath(Global::ParamList& p) throw (std::invalid_argument);

  /// read random number generator initialization value
  void readRandomSeedInitValue(Global::ParamList&) throw (std::invalid_argument&);

  /// read report filename
  void readReportFilename(Global::ParamList&) throw (std::invalid_argument&);

  /// reads name of spacecraft parameter file and check availibility
  void readScParamFilename(Global::ParamList&) throw (std::invalid_argument);

  /// reads name of simulation parameter file and check availibility
  void readSimParamFilename(Global::ParamList&) throw (std::invalid_argument);

  /// reads name of TecPlot trajectory data file
  void readTecPlotFilename(Global::ParamList&) throw (std::invalid_argument);

  /// read VRML output filename
  void readVrmlFilename(Global::ParamList&) throw (std::invalid_argument);

  /*! reads warmstart solution file
   *
   */
  void readWarmstartSolutionFilename(ParamList& p)
    throw (std::invalid_argument);

  /// sets a new value for pseudo random number generator initialization
  void setRandomSeed(const unsigned int seed) {random_seed=seed;};

  /// returns currently set value for pseudo random number generator initialization
  unsigned int getRandomSeed(void) const {return random_seed;};

  /// prints stored actual settings
  void printSettings(void) const;

  /// reset input object
  void reset(void);

}; // end class CInputInTrance

#endif // CINPUTINTRANCE_H_INCLUDED
