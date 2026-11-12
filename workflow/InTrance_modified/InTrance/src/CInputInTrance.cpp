/*!
\file CInputInTrance.cpp
\version 1.7
\brief CInputInTrance class implementation.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 16.05.2006 (AO)

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

// InTrance
#include <configInTrance.h>

// STL
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cassert>

// SpaceLib
#include <Global.h>

// InTrance
#include <CInputInTrance.h>

using namespace std;

//--------------------------------------------------------------------------
CInputInTrance::CInputInTrance(void)
: simulation_parameter_file(),
  spacecraft_parameter_file(),
  nc_topology_file(),
  evolutionary_algorithm_parameter_file(),
  evolutionary_algorithm_warmstart_file(),
  simulation_input_file_ang(),
  clock_approx_data_file(),
  cone_approx_data_file(),
  evolutionary_algorithm_best_solution_file(),
  random_seed_provided(false),
  random_seed(DFLT_RANDOM_SEED),
  simulation_output_file_csv(),
  simulation_output_file_dat(),
  simulation_output_file_gesop(),
  simulation_output_file_wrl(),
  simulation_output_file_ang(),
  report_file(),
  rot_state_traj_file(),
  clock(0.0),
  cone(0.0),
  a_max(0.0),
  e_max(0.0),
  i_max(0.0),
  Omega_max(0.0),
  omega_max(0.0),
  a_min(0.0),
  e_min(0.0),
  i_min(0.0),
  Omega_min(0.0),
  omega_min(0.0),
  a_opt(0.0),
  e_opt(0.0),
  i_opt(0.0),
  Omega_opt(0.0),
  omega_opt(0.0)
{} // constructor

//------------------------------------------------------------------------------
void CInputInTrance::printSettings(void) const
{
  const unsigned int w(37);
  const unsigned int p(10);

  // name of best solution file
  ostringstream ss;
  Output("\n");
  ss << setw(w) << left << "input: best solution file"<<": " 
    << evolutionary_algorithm_best_solution_file << endl;
  Output(ss.str());

  // name of evolutionary algorithm configuration file
  ss.str("");
  ss << setw(w) << left << "input: ea configuration"<<": "
    << evolutionary_algorithm_parameter_file << endl;
  Output(ss.str());

  // evolutionary algorithm warmstart file
  ss.str("");
  ss << setw(w) << left << "input: warmstart file"<<": " 
    << evolutionary_algorithm_warmstart_file << endl;
  Output(ss.str());

  // random number generator seed value
  ss.str("");
  ss << setw(w) << left << "input: random seed given"<<": "
    << (random_seed_provided==true?"YES":"NO") << endl;
  Output(ss.str());

  // used random number generator seed value (given or generated)
  ss.str("");
  ss << setw(w) << left << "input: random seed value"<<": "
    << random_seed << endl;
  Output(ss.str());

  // simulation steering angle output file
  ss.str("");
  ss << setw(w) << left << "input: steering angle output file"<<": " 
    << simulation_output_file_ang << endl;
  Output(ss.str());

  // comma seperated value output file
  ss.str("");
  ss<< setw(w) << left << "input: simulation data file(CSV)"<<": " 
    << simulation_output_file_csv << endl;
  Output(ss.str());

  // TecPlot trajectory data
  ss.str("");
  ss << setw(w) << left << "input: trajectory data (TecPlot)"<<": " 
    << simulation_output_file_dat << endl;
  Output(ss.str());

  // GESOP data
  ss.str("");
  ss << setw(w) << left << "input: GESOP data file"<<": " 
    << simulation_output_file_gesop << endl;
  Output(ss.str());

  // VRML data
  ss.str("");
  ss << setw(w) << left << "input: VRML file"<<": " 
    << simulation_output_file_wrl << endl;
  Output(ss.str());

} // CInputInTrance::printSettings(void)



//------------------------------------------------------------------------------
void CInputInTrance::readAngleParams(Global::ParamList& p)
throw (std::invalid_argument)
{
  const std::string method_name =
    "void CInputInTrance::readAngleParams(Global::ParamList& p)";

  InputParameter("CLOCK_ANGLE", clock, p, method_name, TRAILING_CHARS);
  InputParameter("CONE_ANGLE", cone, p, method_name, TRAILING_CHARS);

  readCsvFilename(p);                                 // data csv file
  readTecPlotFilename(p);         // trajectory data file (for TECPLOT)
  readGESOPFilename(p);

} // void CInputInTrance::readAngleParams(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputInTrance::readBestChromFilename(Global::ParamList& p
                                           ,bool file_existence_check)
  throw (std::invalid_argument)
{
  const std::string method_name = 
    "void readBestChromFilename(Global::ParamList& p)";

  InputFilenameParameter("BEST_CHROM_FILE"
    , evolutionary_algorithm_best_solution_file, p
    , method_name+": missing parameter CTRL_FILE=", TRAILING_CHARS, path
    , true, file_existence_check);
} // void CInputInTrance::readBestChromFilename(Global::ParamList&)

//------------------------------------------------------------------------------
void CInputInTrance::readCommand(Global::ParamList& p)
throw (std::invalid_argument)
{
  const std::string method_name =
    "void CInputInTrance::readCommand(Global::ParamList& p)";

  const std::string err_msg1 = method_name + 
    std::string(": command not found\n");

  std::string s_command;
  InputStringParameter("COMMAND", s_command, p, method_name, TRAILING_CHARS);
  removeChars(s_command, " \"");
  upper(s_command);

  int i= -1, sz1 = (sizeof(COMMANDS)/(sizeof(COMMANDS[0]))),
   sz2 = (sizeof(COMMAND_PARAM_VALUE)/(sizeof(COMMAND_PARAM_VALUE[0])));

  assert(sz1==sz2);

  while ( (++i) < sz1)
  {
    if (s_command==COMMAND_PARAM_VALUE[i])
    {
      what_to_do = COMMANDS[i];
      break;
    }
  }
  if (i==sz1) throw std::invalid_argument(err_msg1);

} // void CInputInTrance::readCommand(Global::ParamList& p)


//-----------------------------------------------------------------------------
void CInputInTrance::readCsvFilename(Global::ParamList& p)
throw (std::invalid_argument&)
{
  const std::string method_name =
    "void CInputInTrance::readCsvFileName(Global::ParamList& p)";
  InputFilenameParameter("SIM_DATA_FILE", simulation_output_file_csv, p
    , method_name + ": SIM_DATA_FILE=", TRAILING_CHARS, path, true, false);
} // void CInputInTrance::readCsvFileName(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputInTrance::readCtrlFilename(Global::ParamList& p) 
throw (std::invalid_argument)
{
  const std::string method_name =
    "void CInputInTrance::readCtrlFilename(Global::ParamList& p)";
  InputFilenameParameter("CTRL_FILE", simulation_output_file_ang, p
    , method_name + ": CTRL_FILE=", TRAILING_CHARS, path, true, false);
} // void CInputInTrance::readTecPlotFilename(Global::ParamList&) 


//-----------------------------------------------------------------------------
void CInputInTrance::readEaParamFilename(Global::ParamList& p)
throw (std::invalid_argument)
{
  const std::string method_name = 
    "void CInputInTrance::readEaParamFilename(Global::ParamList& p)";

  InputFilenameParameter("EA_PARAM_FILE"
    , evolutionary_algorithm_parameter_file, p
    , method_name + ": EA_PARAM_FILE=", TRAILING_CHARS, path, true, true);
} // void CInputInTrance::readEaParamFilename(Global::ParamList&)


//-----------------------------------------------------------------------------
void CInputInTrance::readEvaluateApproxOptimizationParams(
  ParamList& p) throw (std::invalid_argument&)
{
  readCsvFilename(p);                                           // data csv file
  readTecPlotFilename(p);                  // trajectory data file (for TECPLOT)
  readVrmlFilename(p);                                              // VRML file
  readCtrlFilename(p);                         // control values output filename
}


//-----------------------------------------------------------------------------
void CInputInTrance::readEvaluateControlParams(Global::ParamList& p)
throw (std::invalid_argument&)
{
  const std::string method_name = 
    "void CInputInTrance::readEvaluateControlParams(Global::ParamList& p)";

  InputFilenameParameter("CTRL_EVAL_FILE", simulation_input_file_ang,
    p, method_name+": CTRL_EVAL_FILE=", TRAILING_CHARS, path, true, true);

  readCsvFilename(p);                                           // data csv file
  readTecPlotFilename(p);                  // trajectory data file (for TECPLOT)
  readGESOPFilename(p);
  readVrmlFilename(p);                                              // VRML file
  readCtrlFilename(p);                         // control values output filename
} // void CInputInTrance::readEvaluateControlParams(Global::ParamList&)



//------------------------------------------------------------------------------
void CInputInTrance::readFileIdentifier(Global::ParamList& p) 
throw (std::invalid_argument&)
{
  const std::string method_name = 
    "void CInputInTrance::readFileIdentifier(Global::ParamList& p)";

  const std::string err_msg1 = method_name + 
    std::string(": not an InTrance 1.7 input file\n");

  std::string s_identifier;

  InputStringParameter("FILE_TYPE", s_identifier, p, method_name,
    TRAILING_CHARS);
  removeChars(s_identifier, string(" \""));          // remove quotes and spaces
  upper(s_identifier);

  if (s_identifier!="INTRANCE1.7INPUTFILE") throw
    std::invalid_argument(err_msg1);

} // void CInputInTrance::readFileIdentifier(Global::ParamList&)


//-----------------------------------------------------------------------------
void CInputInTrance::readGESOPFilename(Global::ParamList& p)
throw (std::invalid_argument&)
{
  const std::string method_name = 
    "void readGESOPFilename(Global::ParamList& p)";

  InputFilenameParameter("GESOP_FILE", simulation_output_file_gesop,
    p, method_name+": GESOP_FILE=", TRAILING_CHARS, path, true, false);
} // void CInputInTrance::readGESOPFilename(Global::ParamList& p)


//-----------------------------------------------------------------------------
void CInputInTrance::readLslMinMaxParams(Global::ParamList& p)
 throw (std::invalid_argument)
{
  const std::string method_name =
    "void CInputInTrance::readLslMinMaxParams(Global::ParamList& p)";

  InputParameter("A_MAX", a_max, p, method_name, TRAILING_CHARS);
  InputParameter("A_MIN", a_min, p, method_name, TRAILING_CHARS);
  InputParameter("E_MAX", e_max, p, method_name, TRAILING_CHARS);
  InputParameter("E_MIN", e_min, p, method_name, TRAILING_CHARS);
  InputParameter("I_MAX", i_max, p, method_name, TRAILING_CHARS);
  InputParameter("I_MIN", i_min, p, method_name, TRAILING_CHARS);
  InputParameter("OMEGA_MAX", Omega_max, p, method_name, TRAILING_CHARS);
  InputParameter("OMEGA_MIN", Omega_min, p, method_name, TRAILING_CHARS);
  InputParameter("OM_MAX", omega_max, p, method_name, TRAILING_CHARS);
  InputParameter("OM_MIN", omega_min, p, method_name, TRAILING_CHARS);

  readCsvFilename(p);                                           // data csv file
  readTecPlotFilename(p);                  // trajectory data file (for TECPLOT)
  readGESOPFilename(p);
  readVrmlFilename(p);                                              // VRML file
  readCtrlFilename(p);                         // control values output filename

} // void CInputInTrance::readLslMinMaxParams(Global::ParamList&)


//-----------------------------------------------------------------------------
void CInputInTrance::readLslOptParams(Global::ParamList& p)
throw (std::invalid_argument)
{
  const std::string method_name =
    "void CInputInTrance::readLslOptParams(Global::ParamList& p)";

  InputParameter("A_OPT", a_opt, p, method_name, TRAILING_CHARS);
  InputParameter("E_OPT", e_opt, p, method_name, TRAILING_CHARS);
  InputParameter("I_OPT", i_opt, p, method_name, TRAILING_CHARS);
  InputParameter("OMEGA_OPT", Omega_opt, p, method_name, TRAILING_CHARS);
  InputParameter("OM_OPT", omega_opt, p, method_name, TRAILING_CHARS);

  readCsvFilename(p);                                           // data csv file
  readTecPlotFilename(p);                  // trajectory data file (for TECPLOT)
  readGESOPFilename(p);
  readVrmlFilename(p);                                              // VRML file
  readCtrlFilename(p);                         // control values output filename

} // void CInputInTrance::readLslOptParams(Global::ParamList&)


//-----------------------------------------------------------------------------
void CInputInTrance::readNcParamFilename(Global::ParamList& p)
throw (std::invalid_argument)
{
  const std::string method_name = 
    "void readNcParamFilename(Global::ParamList& p)";

  InputFilenameParameter("NC_PARAM_FILE", nc_topology_file,
    p, method_name+": NC_PARAM_FILE=", TRAILING_CHARS, path, true, true);

} // void CInputInTrance::readNcParamFilename(Global::ParamList& p)


//-----------------------------------------------------------------------------
void CInputInTrance::readOptimizeApproxParams(ParamList& p)
  throw (std::invalid_argument&)
{
  const string prefix("CInputInTrance::ReadApproxOptimizationParams(ParamList&)");

  // best chromosome file
  InputFilenameParameter("BEST_CHROM_FILE", evolutionary_algorithm_best_solution_file,
    p, prefix+"BEST_CHROM_FILE=", TRAILING_CHARS, path, true, false);

  // evolutionary algorithm parameter file
  InputFilenameParameter("EA_PARAM_FILE", evolutionary_algorithm_parameter_file,
    p, prefix+"EA_PARAM_FILE=", TRAILING_CHARS, path, true, true);

  // neurocontroller file (steering angle approximation)  
  InputFilenameParameter("NC_PARAM_FILE", nc_topology_file,
    p, prefix+"NC_PARAM_FILE=", TRAILING_CHARS, path, true, true);

  // input steering angles file (steering angle approximation)
  InputFilenameParameter("CTRL_EVAL_FILE", simulation_input_file_ang,
    p, prefix+"CTRL_EVAL_FILE=", TRAILING_CHARS, path, true, true);

  // input reference trajectory data file (steering angle approximation)
  InputFilenameParameter("INPUT_TRAJ_DATA_FILE", reference_trajectory_file,
    p, prefix+"INPUT_TRAJ_DATA_FILE=", TRAILING_CHARS, path, true, true);

  // possible warmstart solution
  readWarmstartSolutionFilename(p);

  // input clock angle approx. parameters
  InputFilenameParameter("CLOCK_APPROX_DATA", clock_approx_data_file,
    p, prefix+"CLOCK_APPROX_DATA=", TRAILING_CHARS, path, true, false);

  // input cone angle approx. parameters
  InputFilenameParameter("CONE_APPROX_DATA", cone_approx_data_file,
    p, prefix+"CONE_APPROX_DATA=", TRAILING_CHARS, path, true, false);

  // rotational state trajectory data file
  InputFilenameParameter("ROT_TRAJ_FILE", rot_state_traj_file,
    p, prefix+"ROT_TRAJ_FILE=", TRAILING_CHARS, path, true, false);

  // GESOP file
  InputFilenameParameter("GESOP_FILE", simulation_output_file_gesop,
    p, prefix+"GESOP_FILE=", TRAILING_CHARS, path, false, false);

  readCsvFilename(p);                                           // data csv file
  readTecPlotFilename(p);                  // trajectory data file (for TECPLOT)
  readVrmlFilename(p);                                              // VRML file
  readCtrlFilename(p);                         // control values output filename
} // void ReadApproxOptimizationParams(ParamList& params)


//------------------------------------------------------------------------------
void CInputInTrance::readOptimizeNcParams(Global::ParamList& p) 
throw (std::invalid_argument)
{
  const std::string method_name = 
    "void CInputInTrance::readOptimizeNcParams(Global::ParamList& p)";

  readNcParamFilename(p);                 // read neurocontroller parameter file
  readEaParamFilename(p);      // read evolutionary algorithm parameter filename

  switch (what_to_do)
  {
  case OPTIMIZE_NC:
    readBestChromFilename(p);
    readWarmstartSolutionFilename(p); // read possible warmstart filename
    break;
  case EVALUATE_NC:
    readBestChromFilename(p, true); // read name of best chromosome file
    // and check for existence of file
    break;
  default:
    throw std::invalid_argument(method_name + 
      ": invalid InTrance command\n");
    break;
  }

  readCsvFilename(p);                                           // data csv file
  readTecPlotFilename(p);                  // trajectory data file (for TECPLOT)
  readGESOPFilename(p);
  readVrmlFilename(p);                                              // VRML file
  readCtrlFilename(p);                         // control values output filename
} // void CInputInTrance::readOptimizeNcParams(Global::Param&)


//--------------------------------------------------------------------------
void CInputInTrance::readParams(Global::ParamList& p)
{
  const string method_name = 
    "CInputInTrance::Read(Global::ParamList& parameters):";

  readPath(p);
  readFileIdentifier(p);
  readRandomSeedInitValue(p);
  readCommand(p);
  // read input files and check their existence
  readScParamFilename(p);                                  // s/c parameter file
  readSimParamFilename(p);                          // simulation parameter file
  readReportFilename(p);                                     // report file name

  switch (what_to_do)
  {
  case LSL_MINMAX:
    try {readLslMinMaxParams(p);}
    catch (std::invalid_argument& exc)
    {throw std::invalid_argument(method_name+exc.what());}
    break;

  case LSL_OPT:
    try {readLslOptParams(p);}
    catch (std::invalid_argument& exc)
    {throw std::invalid_argument(method_name+exc.what());}
    break;

  case ANGLES:
    try {readAngleParams(p);}
    catch (std::invalid_argument& exc)
    {throw std::invalid_argument(method_name+exc.what());}
    break;

  case OPTIMIZE_NC:
  case EVALUATE_NC:
    try {readOptimizeNcParams(p);}
    catch (std::invalid_argument& exc)
    {throw std::invalid_argument(method_name+exc.what());}
    break;

  case EVALUATE_CONTROL:
    try {readEvaluateControlParams(p);}
    catch (std::invalid_argument& exc)
    {throw std::invalid_argument(method_name+exc.what());}
    break;

  case OPTIMIZE_APPROX:
  case EVALUATE_APPROX:
    try {readOptimizeApproxParams(p);}
    catch (std::invalid_argument& exc)
    {throw std::invalid_argument(method_name+exc.what());}
    break;

  default:break;
  } // switch(what_to_do)

} // end CInputInTrance::Read(ParamList& parameters)


//------------------------------------------------------------------------------
void CInputInTrance::readPath(Global::ParamList& p)
throw (std::invalid_argument)
{
  try
  {
    InputStringParameter("PATH", path, p, "", "");
  }
  catch (std::invalid_argument) {};
} // void CInputInTrance::inputPath(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputInTrance::readRandomSeedInitValue(Global::ParamList& p)
throw (std::invalid_argument&)
{
  const std::string method_name =
    "void CInputInTrance::readRandomSeedInitValue(Global::ParamList& p)";
  try 
  {
    random_seed_provided = true;
    InputParameter("RANDOM_SEED", random_seed, p, method_name
      , TRAILING_CHARS, true);
  }
  catch (std::invalid_argument)
  {
    random_seed_provided = false;
  }
} // void CInputInTrance::readRandomSeedValue(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputInTrance::readReportFilename(Global::ParamList& p)
throw (std::invalid_argument&)
{
  const std::string method_name =
    "void CInputInTrance::readReportFilename(Global::ParamList&)";

  // report file name
  try
  {
    InputStringParameter("REPORT_FILE", report_file, p, "", TRAILING_CHARS);
    report_file=path+report_file;
  }
  catch (std::invalid_argument) {Warning("No report file parameter!\n");};
  if (report_file.empty()==true) Warning("No report file name!\n");

} // void CInputInTrance::readReportFilename(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputInTrance::readScParamFilename(Global::ParamList& p)
 throw (std::invalid_argument)
{
  const std::string method_name = "void readScParamFile(Global::ParamList& p)";

  InputFilenameParameter("SC_PARAM_FILE", spacecraft_parameter_file, p,
    method_name+": SC_PARAM_FILE=", TRAILING_CHARS, path, true, true);

} // void readScParamFile(Global::ParamList& p)


//------------------------------------------------------------------------------
void CInputInTrance::readSimParamFilename(Global::ParamList& p) 
throw (std::invalid_argument)
{
  const std::string method_name = 
    "void readSimParamFilename(Global::ParamList& p)";

  InputFilenameParameter("SIM_PARAM_FILE", simulation_parameter_file, p,
    method_name+": SIM_PARAM_FILE=", TRAILING_CHARS, path, true, true);

} // void CInputInTrance::readSimParamFilename(Global::ParamList&)


//------------------------------------------------------------------------------
void CInputInTrance::readTecPlotFilename(Global::ParamList& p) 
throw (std::invalid_argument)
{
  const std::string method_name =
    "void CInputInTrance::readTecPlotFilename(Global::ParamList& p)";
  InputFilenameParameter("TRAJ_DATA_FILE", simulation_output_file_dat, p
    , method_name + ": TRAJ_DATA_FILE=", TRAILING_CHARS, path, true, false);
} // void CInputInTrance::readTecPlotFilename(Global::ParamList&) 


//------------------------------------------------------------------------------
void CInputInTrance::readVrmlFilename(Global::ParamList& p) 
throw (std::invalid_argument)
{
  const std::string method_name =
    "void CInputInTrance::readVrmlFilename(Global::ParamList& p)";
  InputFilenameParameter("VRML_FILE", simulation_output_file_wrl, p
    , method_name + ": VRML_FILE=", TRAILING_CHARS, path, true, false);
} // void CInputInTrance::readTecPlotFilename(Global::ParamList&) 


//------------------------------------------------------------------------------
void CInputInTrance::readWarmstartSolutionFilename(ParamList& p) 
throw (std::invalid_argument)
{
  const std::string method_name = 
    "void CInputInTrance::readWarmstartSolutionFilename(ParamList& p)";

  std::string s_warmstart;
  try
  {
    InputStringParameter("COLDSTART", s_warmstart, p, "", TRAILING_CHARS);
  }
  catch (std::invalid_argument) {} 
  upper(s_warmstart);
  if (s_warmstart=="NO")
  {
    InputFilenameParameter("WARMSTART_CHROM_FILE"
      , evolutionary_algorithm_warmstart_file, p
      , method_name+": WARMSTART_CHROM_FILE=", TRAILING_CHARS, path, true, true);
  }
} // void CInputInTrance::readWarmstartSolution(ParamList& p, const std::string& path)


//------------------------------------------------------------------------------
void CInputInTrance::reset(void)
{
  path.erase();
  random_seed_provided = false;
  random_seed = DFLT_RANDOM_SEED;
  simulation_parameter_file.erase();
  spacecraft_parameter_file.erase();
  nc_topology_file.erase();
  snc_topology_file.erase();
  evolutionary_algorithm_parameter_file.erase();
  evolutionary_algorithm_warmstart_file.erase();
  simulation_input_file_ang.erase();
  clock_approx_data_file.erase();
  cone_approx_data_file.erase();
  evolutionary_algorithm_best_solution_file.erase();
  simulation_output_file_csv.erase();
  simulation_output_file_dat.erase();
  simulation_output_file_gesop.erase();
  simulation_output_file_wrl.erase();
  simulation_output_file_ang.erase();
  reference_trajectory_file.erase();
  report_file.erase();
  rot_state_traj_file.erase();
  clock = 0.0;
  cone = 0.0;
  a_max = 0.0;
  e_max = 0.0;
  i_max = 0.0;
  Omega_max = 0.0;
  omega_max = 0.0;
  a_min = 0.0;
  e_min = 0.0;
  i_min = 0.0;
  Omega_min = 0.0;
  omega_min = 0.0;
  a_opt = 0.0;
  e_opt = 0.0;
  i_opt = 0.0;
  Omega_opt = 0.0;
  omega_opt = 0.0;
}
