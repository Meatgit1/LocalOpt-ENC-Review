/*!
\file CInTrance.cpp
\version 1.7
\brief CInTrance class implementation.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 23.10.2006 (AO)

\b Purpose:
   This file provides the interface for the class CInputInTrance.


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

#include <conio.h>

// STL
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdlib>

// Java Native Interface
#if defined JAVA_GUI
 #include <jni.h>
#endif

// SpaceLib
#include <Global.h>

// ENCLib
#include <CInputEA.h>
#include <CEAFitness.h>
#include <CIndividual.h>
#include <CEA.h>
#include <CNC.h>

// SCLib
#include <CSCAttitudeHistory.h>
#include <CSC.h>
#include <NativeBridgeProbe.h>

//InTrance
#include <InTranceConst.h>
#include <CInputInTrance.h>
#include <CInputSimulation.h>
#include <CInTranceFitness.h>
#include <CApproxFitness.h>
#include <CSimulation.h>
//#include <CGTOC2Sim.h>
//#include <CCONSEPSim.h>
#include <CInTrance.h>
#include <CInputNC.h>

#if defined (COMMONCPP)
 #include <TQueue.h>
#endif

using namespace std;
using namespace Global;
using namespace MathLib;
using namespace SpaceLib;

#undef CONSEP // activate for CONSEP-problem specific code

#if defined JAVA_GUI
extern "C" {

 #if defined COMMONCPP

  extern TransferDataItemBuffer* dataBuffer;
  extern EAStatusRec EA_status;
  extern SimStatusRec SIM_status;

 #else #error CInTrance.cpp => No multithreading library for transfer buffer purposes
 #endif //COMMONCPP
}
#endif //JAVA_GUI


const int ADDITIONAL_CHROMOSOME_DATA(5);

extern FitnessFunctionType  jt;  // global in CInTranceFitness.cpp
extern int jt2;                  // global in CInTranceFitness.cpp
extern int jt3;                  // global in CInTranceFitness.cpp

// time variable to estimate performance parameters like
// - reproductions per second (rep/sec)
static clock_t time_of_last_call = 0;
static int n_last_reproductions = 0;

// Isolated smoke accounting: monitor/report evaluations are not search calls.
static bool enc_smoke_accounting=false;
static bool enc_smoke_final_report=false;
static unsigned long enc_smoke_search_monitors=0;
static unsigned long enc_smoke_final_monitors=0;
static unsigned long enc_smoke_final_evaluations=0;

bool CInTrance::m_break;
CSimulation* CInTrance::sim;
CEA* CInTrance::ea;
CInputInTrance CInTrance::input_inp;
CInputEA CInTrance::input_eap;

CSimulation* simulation;
double accuracy_factor;


//------------------------------------------------------------------------------
void SaveBestSolution(CEA* ea, const char* filename)
{
  if (ea!=0) ea->SaveBestSolution(filename); 
} // void SaveTempSolution(CEA* ea, const char* filename)


//------------------------------------------------------------------------------
void SaveBestSolution(CEA* ea, const std::string filename)
{
  SaveBestSolution(ea, filename.c_str());
} // void SaveBestSolution(CEA* ea, const std::string filename)


//------------------------------------------------------------------------------
void SaveTempSolution(CEA* ea)
{
  SaveBestSolution(ea, "tmp.eac");
} // void SaveTmpSolution(CEA* ea)


//------------------------------------------------------------------------------
/// handler of evolutionary algorithm events
void CInTrance::handleEaEvents(const EAEvent ea_event)
{
  double best_val,limit,up_v,lo_v, spread;
  int cl;
  CEA* ea;
  CSimulation* sim = getCurrentSimObj();

  switch (ea_event)
  {
  case NOEVENT:break;
  case UNSUCCESSFUL_EPOCH:
    break;
  case LAST_UNSUCCESSFUL_EPOCH:
    break;
  case UNSUCCESSFUL_SSS:
    break;
  case SUCCESSFUL_EPOCH:
  case SUCCESSFUL_SSS:
    SaveBestSolution(CInTrance::ea, 
      CInTrance::input_inp.evolutionary_algorithm_best_solution_file);
    break;
  case HYPERCUBE_SIZE_CHANGE:
    break;
  case NEW_POPULATION_INITIALIZED:

    ea = getCurrentEaObj();
    cl = ea->getChromosomeLength();
    spread = (ea->getPopulationSize()-1)*ea->getCurrentHypercubeSize();

    if (getCurrentSimObj()->getMinStartTime() < 
      getCurrentSimObj()->getMaxStartTime())
    {
      best_val = ea->getBestSolutionEver()->getAllele(cl-5);
      up_v = best_val+spread/2.0;
      lo_v = best_val-spread/2.0;
      limit = ea->getHypercubeInitialSize();

      if (up_v > limit)
      {
        up_v = limit;
        lo_v = std::max(-limit, (up_v-spread) );
      } // if (up_v > limit)

      if (lo_v < (-limit))
      {
        lo_v = (-limit);
        up_v = std::min(limit, (lo_v + spread) );
      } // if (lo_v < (-limit))

      ea->uniformSpreadSolutionAlleles(cl-5, lo_v, up_v);
    } // if (getCurrentSimObj()->getMinStartTime() <...

    if (!sim->getCurrentSC().isSolarSail())
    {
      if (sim->getCurrentSC().getPropSys()->getMinInitialPropellantMass() < sim->getCurrentSC().getPropSys()->getMaxInitialPropellantMass())
      {
        best_val = ea->getBestSolutionEver()->getAllele(cl-1);
        up_v = best_val+spread/2.0;
        lo_v = best_val-spread/2.0;
        limit = ea->getHypercubeInitialSize();

        if (up_v > limit)
        {
          up_v = limit;
          lo_v = std::max(-limit, (up_v-spread) );
        } // if (up_v > limit)

        if (lo_v < (-limit))
        {
          lo_v = (-limit);
          up_v = std::min(limit, (lo_v + spread) );
        } // if (lo_v < (-limit))

        ea->uniformSpreadSolutionAlleles(cl-1, lo_v, up_v);
      } // if (sim->getCurrentSC().getPropSys()->getMinInitialPropellantMass()
    }

    break; // case NEW_POPULATION_INITIALIZED:
  default:
    break;
  }
} // void EaHandler(const EAEvent event)


const std::string MENU_SCREEN[] =
{//1-------10--------20--------30--------40--------50--------60--------70--------80
  "-----------------------------  keyboard commands -------------------------------",
  "                                    ee   : end current epoch",
  "a : end current optimization        ei   : display ea status",
  "m : menu                            em   : change mutation rate",
  "s : save best solution              eh   : change hypercube shrinking factor",
  "p : pause                           ebe  : change max. # of bad epochs",
  "",
  "--------------------------------------------------------------------------------"
};




//------------------------------------------------------------------------------
void Pause(void)
{
  cout << "\n\npausing....";
  do {} while (_kbhit()==0);
  cout << "ended\n";
} // Pause(void)

// ----------------------------------------------------------------------------
void ModifyHypercubeShrinkingFactor(CEA* ea, double inc1, double inc2,int ow)
{
  int key = 0;
  cout << "\nmodify hypercube shrinking factor\n";
  cout << "use +/- ("<<inc1<<" step size) or PgUp/PgDn ("<<inc2<<" step size)"
    << "and ESC if finished\n";
  cout << "-------------------------------------------------------------------------\n";
  do 
  {
    key = _getch();
    switch (key)
    {
    case 43:          // "+"
      ea->incHypercubeShrinkFactor(inc1);
      ea->printHypercubeShrinkFactorSetting(ow);
      break;
    case 45:          // "-"
      ea->decHypercubeShrinkFactor(inc1);
      ea->printHypercubeShrinkFactorSetting(ow);
      break;

    case 224:         // special key
      key = _getch();
      switch (key)
      {
        //------------------------------// "eh + page up"
      case 73:
        ea->incHypercubeShrinkFactor(inc2);
        ea->printHypercubeShrinkFactorSetting(ow);
        break;
        //------------------------------// "eh + page down"
      case 81:
        ea->decHypercubeShrinkFactor(inc2);
        ea->printHypercubeShrinkFactorSetting(ow);
        break;
      default:
        break;
      }
    default:
      break;
    }
  } while (key!=KEY_ESCAPE);
}

// ----------------------------------------------------------------------------
void ModifyMutationRate(CEA* ea, double inc1, double inc2,int ow)
{
  int key = 0;
  cout << "\nmodify mutation rate\n";
  cout << "use +/- ("<<inc1<<" step size) or PgUp/PgDn ("<<inc2<<" step size)"
    << "and ESC if finished\n";
  cout << "-------------------------------------------------------------------------\n";
  do 
  {
    key = _getch();
    switch (key)
    {
    case KEY_PLUS:          // "+"
      ea->incChromMutationRate(inc1);
      ea->printMutationRateSetting(ow);
      break;
    case KEY_MINUS:          // "-"
      ea->decChromMutationRate(inc1);
      ea->printMutationRateSetting(ow);
      break;

    case 224:         // special key
      key = _getch();
      switch (key)
      {
        //------------------------------// "eh + page up"
      case 73:
        ea->incChromMutationRate(inc2);
        ea->printMutationRateSetting(ow);
        break;
        //------------------------------// "eh + page down"
      case 81:
        ea->decChromMutationRate(inc2);
        ea->printMutationRateSetting(ow);
        break;
      default:
        break;
      }
    default:
      break;
    }
  } while (key!=KEY_ESCAPE);
}

// ----------------------------------------------------------------------------
void ModifyBadEpochCount(CEA* ea, int ow)
{
  int key = 0;
  cout << "\nmodify max. bad epoch count\n";
  cout << "use +/- (1 step size)" << " and ESC if finished\n";
  cout << "-------------------------------------------------------------------------\n";
  do 
  {
    key = _getch();
    switch (key)
    {
    case 43:          // "+"
      ea->setMaxBadEpochsCount(ea->getMaxBadEpochsCount()+1);
      ea->printMaxBadEpochsSetting(ow);
      break;
    case 45:          // "-"
      ea->setMaxBadEpochsCount(ea->getMaxBadEpochsCount()-1);
      ea->printMaxBadEpochsSetting(ow);
      break;
    default:
      break;
    }
  } while (key!=KEY_ESCAPE);
}

// ----------------------------------------------------------------------------
void CInTrance::handleKeyboardInputs(const int first_key, const int second_key)
{
  const int OW = 50;
  int key=0;
  switch (second_key)
  {
  case 0:
    switch (first_key)
    {

    case KEY_A:                   // "a"
      setAbort();
      break;

    case KEY_E:                   // "e"
      key = _getch();
      switch (key)
      {
      //------------------------------// "ee"
      case KEY_E:               
        ea->abortCurrentEpoch();
        break;
      //------------------------------// "eh"
      case KEY_H:
        ModifyHypercubeShrinkingFactor(ea, HCSF_INC1, HCSF_INC2, OW);
        break;

      //------------------------------// "ei"
      case KEY_I:                           
        PrintEaSettings(ea);
        break;

      //------------------------------// "em"
      case KEY_M:
        ModifyMutationRate(ea, MR_INC1, MR_INC2, OW);
        break;

      //------------------------------// "eb"
      case KEY_B:                        
        key = _getch();
        switch (key)
        {
        //------------------------------// "ebe"
        case KEY_E:
          ModifyBadEpochCount(ea, OW);
          break;
        default:
          break;
        }
        break;
      default:
        break;
      }
      break;

    case KEY_H:                   // "h"
      break;
    case KEY_M:                   // "m"
      printHelpScreen();
      break;
    case KEY_P:                 // "p"
      Pause();
      break;
    case KEY_S:                 // "s"
      if (sim!=0) 
      {
        SaveTempSolution(ea);
        switch (input_inp.what_to_do)
        {
        case OPTIMIZE_NC:
        case EVALUATE_NC:
          if (ea != 0) 
          {
            CIndividual solution(ea,"tmp.eac");
            // final solution evaluation
            sim->evaluate(solution, false, true);
          }        
          sim->writeOutputFiles("tmp.csv", "temp.dat", "tmp.gesop.txt"
            , "tmp.wrl", "tmp.ang");
          break;
        default:break;
        } // switch (input_inp.what_to_do)

      }
      break;
    case 117:               // "u"
      break;
    }
    break;
  default:
    break;
  }
}

// ----------------------------------------------------------------------------
void CInTrance::printHelpScreen(void)
{  
  cout << endl;
  for (int i=0;i<sizeof(MENU_SCREEN)/sizeof(MENU_SCREEN[0]);++i)
    cout << MENU_SCREEN[i] << endl;

}


//------------------------------------------------------------------------------
void CalcAccuracyFactor(double hc_init, double hc_curr, double hc_min)
{
  if (simulation != 0)
  {
    if (simulation->useDAC())
    {
      const double hc_max = simulation->getDACMaxHcSize();
      const double hc_min = simulation->getDACFinalHcSize();
      const double f_max  = simulation->getDACMaxFactor();
      const double f_min  = simulation->getDACMinFactor();
      const double f = 
        std::max(0.0, (f_max-f_min)*(std::min(hc_curr,hc_max)-hc_min)/(hc_max-hc_min))+f_min;
      simulation->setAccuracyFactor(f/*/(hc_max-hc_min)*/);
    } // if (simulation->useDAC())
  } // if (simulation != 0)
} // void CalcAccuracyFactor(double hc_init, double hc_curr, double hc_min)


//------------------------------------------------------------------------------
#if defined COMMONCPP

SLOutput& GUIOutput::operator<<(const char* message)
{
  if (message != 0)
  {
    if (dataBuffer != 0)
    {
      dataBuffer->put(TransferDataItem(0, 0, 0, 0, 0, 0, 0, 0, 0, string(message), 0.0));
    }
    return (SLOutput::operator <<(message));
  }
  return *this;
}

GUIOutput java_output;
#endif //COMMONCPP

#ifdef COMMONCPP
/*! \brief set subfitness string fields of transfer data structure
 *
 * You may specifiy the field width and precision of subfitness value output,
 * while field width of the particular subfitness name will be determined from
 * the subfitness names string lengths, so the longest subfitness name sets
 * the field width for all subfitness names.
 * \param ind pointer to individual and its subfitness values
 * \param item exchange data structure
 * \param w2 desired field width of subfitness value
 * \param prec desired output precision of subfitness value
 */
void StoreSubfitnessValues(CIndividual* ind, TransferDataItem& item, int w2
                           , int prec)
{
  std::ostringstream text("t");

  // subfitness name field width (w1), its value (w2) and precision (prec)
  int w1(0);
  for (int i=0;i<N_SUB_FITNESS;++i)                    // determine field width
    w1=Global::max(w1,int(NAMES_SUBFITNESS[i].length()));

  // set output format
  text << setfill(' ') << left << scientific << showpos << setprecision(prec);

  for (int i=0;i<N_SUB_FITNESS;++i)           // loop over subfiness fitness values
  {
    text.seekp(0);text.str("");             // reset string stream
    text << setw(w1) << NAMES_SUBFITNESS[i];                  // subfitness name
    text << setw(1)  << ":";                             // seperator character
    text << setw(10) << ind->getFitnessObject()->getFitness(i);        // value
    item.m_fitness_values.push_back(text.str());        // store created string
  }
}

//-----------------------------------------------------------------------------
void SubmitTransferData(CIndividual* ind, CSimulation* sim, std::string& msg,
                        float rep_per_sec)
{
  Trajectory traj                  = 0;
  Trajectory target_traj           = 0;
  DataArray  clock_angle           = 0;
  DataArray  cone_angle            = 0;
  DataArray  thrust                = 0;

  if (sim != 0)
  { 
    sim->evaluate(*ind, false, true);                      // monitor evaluation
    sim->saveLastTrajectory(*ind);                        // save
    traj               = sim->getLastTrajectoryData();    // s/c trajectory data
    target_traj        = sim->getLastTrajectoryTargetData(); // s/c target
    clock_angle        = sim->getLastTrajectoryClockAngleData();
    cone_angle         = sim->getLastTrajectoryConeAngleData();
    thrust             = sim->getLastTrajectoryThrustData();
  }

  // All the above standing calls to "get..."-methods of CSimulation class may
  // return zero but there's no special error handling code needed
  // because the data transfer routines of "InTranceLib.cpp" ignore data array
  // pointers set to zero. So in case of errors or no data at all there
  // is no serious threat to program execution.

  TransferDataItem item(traj
                       ,0
                       ,target_traj
                       ,clock_angle
                       ,0
                       ,cone_angle
                       ,0
                       ,thrust
                       ,0
                       ,msg
                       ,rep_per_sec);

  // transfer subfitness values (fitness value field width=10; precision=5)
  StoreSubfitnessValues(ind, item, 10, 5); 
  if (sim != 0) sim->writeMessageToReportFile(msg);    // write log file entry
  dataBuffer->put(item);   // put completed transfer item into storage structure
}
#endif // COMMONCPP

//--------------------------------------------------------------------------------
inline void FitnessText_a(CIndividual* i, CSimulation* sim, const int width, 
                     const int prec, std::ostringstream& oss)
{
  if (sim->useTargetOrbitalElement(SEMIMAJORAXIS)==true)
  {
    oss<<setw(width)<<setfill(' ')<< "a=" << scientific
      <<setprecision(prec)<<i->getFitnessObject()->getFitness(A_FITNESS);
  }
}
//--------------------------------------------------------------------------------
inline void FitnessText_e(CIndividual* i, CSimulation* sim, const int width, 
                     const int prec, std::ostringstream& oss)
{
  if (sim->useTargetOrbitalElement(ECCENTRICITY)==true)
  {
    oss<<setw(width)<<setfill(' ')<< "e=" << scientific
      <<setprecision(prec)<<i->getFitnessObject()->getFitness(E_FITNESS);
  }
}

//--------------------------------------------------------------------------------
inline void FitnessText_O(CIndividual* i, CSimulation* sim, const int width, 
                     const int prec, std::ostringstream& oss)
{
  if (sim->useTargetOrbitalElement(NODE)==true)
  {
    oss<<setw(width)<<setfill(' ')<< "Om=" << scientific
      <<setprecision(prec)<<i->getFitnessObject()->getFitness(NODE_FITNESS);
  }
}

//--------------------------------------------------------------------------------
inline void FitnessText_i(CIndividual* i, CSimulation* sim, const int width, 
                     const int prec, std::ostringstream& oss)
{
  if (sim->useTargetOrbitalElement(INCLINATION)==true)
  {
    oss<<setw(width)<<setfill(' ')<< "i=" << scientific
      <<setprecision(prec)<<i->getFitnessObject()->getFitness(I_FITNESS);
  }
}

//--------------------------------------------------------------------------------
inline void FitnessText_o(CIndividual* i, CSimulation* sim, const int width, 
                     const int prec, std::ostringstream& oss)
{
  if (sim->useTargetOrbitalElement(PERIAPSIS)==true)
  {
    oss<<setw(width)<<setfill(' ')<< "om=" << scientific
      <<setprecision(prec)<<i->getFitnessObject()->getFitness(PERI_FITNESS);
  }
}

//--------------------------------------------------------------------------------
inline void FitnessText_M(CIndividual* i, CSimulation* sim, const int width, 
                     const int prec, std::ostringstream& oss)
{
  if (sim->useTargetOrbitalElement(ANOMALY)==true)
  {
    oss<<setw(width)<<setfill(' ')<< "M=" << scientific
      <<setprecision(prec)<<i->getFitnessObject()->getFitness(M_FITNESS);
  }
}

//--------------------------------------------------------------------------------
void PrintReproduction(      CIndividual* ind,
                       const double improvement,
                       const int n_epc_repr,         // # of epoch reproductions
                       const int n_abs_repr,      // absolute # of reproductions
                       void* data)
{
  // get system time and calculate time since last call
  clock_t time_of_new_call = clock();
  float diff_time = (float(time_of_new_call)-float(time_of_last_call))/float(CLOCKS_PER_SEC);

  CSimulation* sim = static_cast<CSimulation*>(data);
#ifndef COMMONCPP  
  if (sim!=0)
  {
    sim->evaluate(*ind, false, true);
    if (enc_smoke_accounting)
    {
      if (enc_smoke_final_report) ++enc_smoke_final_monitors;
      else ++enc_smoke_search_monitors;
    }
  }
#endif

  ostringstream text;

  // estimate the reproductions per second ratio 
  float repro_per_sec = 0.0;
  if (diff_time > 0.0) repro_per_sec = ((n_abs_repr-n_last_reproductions)/diff_time);
  if (repro_per_sec < 0.0) repro_per_sec = 0.0;

  // save absolute # of reproductions
  n_last_reproductions = n_abs_repr;
  
  const double t_fitness=ind->getFitnessObject()->getFitness(TIME_FITNESS);
  const double r_fitness=ind->getFitnessObject()->getFitness(POS_FITNESS);
  const double v_fitness=ind->getFitnessObject()->getFitness(VEL_FITNESS);
  const double a_fitness=ind->getFitnessObject()->getFitness(A_FITNESS);
  const double e_fitness=ind->getFitnessObject()->getFitness(E_FITNESS);
  const double i_fitness=ind->getFitnessObject()->getFitness(I_FITNESS);
  const double O_fitness=ind->getFitnessObject()->getFitness(NODE_FITNESS);
  const double o_fitness=ind->getFitnessObject()->getFitness(PERI_FITNESS);
  const double M_fitness=ind->getFitnessObject()->getFitness(M_FITNESS);
  const double h_fitness=ind->getFitnessObject()->getFitness(H_FITNESS);
  const double impv_fitness=ind->getFitnessObject()->getFitness(IMPV_FITNESS);
  const double rimph_fitness=ind->getFitnessObject()->getFitness(RIMPH_FITNESS);
  const double rimpe_fitness=ind->getFitnessObject()->getFitness(RIMPE_FITNESS);
  const double act_fitness=ind->getFitnessObject()->getFitness(ACT_FITNESS);
  const double mprop_fitness=ind->getFitnessObject()->getFitness(MPROP_FITNESS);
  const double ga_fitness=ind->getFitnessObject()->getFitness(GA_FITNESS);
  const double defl_fitness=ind->getFitnessObject()->getFitness(DEFL_FITNESS);
  const double spi_a_fitness=ind->getFitnessObject()->getFitness(SPI_A_FITNESS);
  const double spi_i_fitness=ind->getFitnessObject()->getFitness(SPI_I_FITNESS);
  const double spi_e_fitness=ind->getFitnessObject()->getFitness(SPI_E_FITNESS);

  const double replaced_fitness=ind->getFitnessObject()->getMixedFitness();

  text<<setiosflags(ios::fixed|ios::showpos) << setprecision(8);
  text<<setw(10)<<setfill(' ')<<scientific<<setprecision(5) << replaced_fitness;
  text<<setw(3)<<setfill(' ')<<"[" <<scientific<<setprecision(5)<< improvement <<"]";


  switch (jt)
  {
  case J_ORB_ELEMENTS:
    text<<setw(4)<<setfill(' ')<<" t=" <<scientific<<setprecision(5)
      <<t_fitness/100.0;
    FitnessText_a(ind, sim, 4, 5, text);
    FitnessText_e(ind, sim, 4, 5, text);
    FitnessText_O(ind, sim, 4, 5, text);
    FitnessText_i(ind, sim, 4, 5, text);
    FitnessText_o(ind, sim, 4, 5, text);
    FitnessText_M(ind, sim, 4, 5, text);
    break;
  default:
    text<<setw(4)<<setfill(' ')<<" t=" <<scientific<<setprecision(5)<<t_fitness/100.0;
    text<<setw(4)<<setfill(' ')<<" r=" <<scientific<<setprecision(5)<<r_fitness;
    break;
  }

  if (sim)
  {
    switch (sim->getTrajectoryType())
    {
    case GTOC2TRAJ:
      text<<setw(4)<<setfill(' ')<<" v=" <<v_fitness;
      break;

    case RENDEZVOUS:
      switch (jt)
      {
      case J_ORB_ELEMENTS:
        break;
      default:
        text<<setw(4)<<setfill(' ')<<" v=" <<v_fitness;
        break;
      }

      switch (jt)
      {
      case J_AEH_AND:
      case J_AEH_OR : text<<setw(4)<<setfill(' ')<<" a=" <<a_fitness;
      case J_EH_AND :
      case J_EH_OR  : text<<setw(4)<<setfill(' ')<<" e=" <<e_fitness;
      case J_H_AND  :
      case J_H_OR   : text<<setw(4)<<setfill(' ')<<" h=" <<h_fitness;
        break;
      default:
        break;
      } // switch (jt)
      break;
    case RIMPACT:
      text<<setw(4)<<setfill(' ')<<"  impv=" <<impv_fitness;
      text<<setw(4)<<setfill(' ')<<"  imph=" <<rimph_fitness;
      text<<setw(4)<<setfill(' ')<<"  impe=" <<rimpe_fitness;
      text<<setw(4)<<setfill(' ')<<"  defl=" <<setprecision(0)<<defl_fitness;
      text<< " km"<<setprecision(8);
      break;
    case ACTPROBLEM:
      text<<setw(4)<<setfill(' ')<<" act=" << act_fitness;
      break;
    case SPIPROBLEM:
      text<<setw(4)<<setfill(' ')<<"  SPIa=" <<spi_a_fitness;
      text<<setw(4)<<setfill(' ')<<"  SPIi=" <<spi_i_fitness;
      text<<setw(4)<<setfill(' ')<<"  SPIe=" <<spi_e_fitness;
      break;
    default:break;
    }
  } // if (sim)

  if (sim!=0)
  {
    if (sim->getCurrentSC().isSolarSail() != true)
      text<<setw(4)<<setfill(' ')<<" m=" 
          << scientific << setprecision(2) << mprop_fitness;
  }
  else text<<setw(4)<<setfill(' ')<<" m=" 
           << scientific << setprecision(2) << mprop_fitness;
  text<<resetiosflags(ios::fixed|ios::showpos);

  // print reproduction rate (in reproductions per second)
  text<<setw(3)<<setfill(' ')<<" (" <<setw(6)<<setfill('0')<<n_epc_repr
      <<"/"<<setw(6)<<setfill('0')<<n_abs_repr<<")";
  text << " (" << setw(5) << fixed<< setprecision(1) 
    << repro_per_sec << " rep/sec)" << endl;

#if defined JAVA_GUI
  SubmitTransferData(ind, sim, text.str(), repro_per_sec);
  if (sim!=0) sim->fillStatusRec();
  EA_status.reproduction = n_abs_repr;
  EA_status.epoch_reproduction = n_epc_repr;
  EA_status.improvement = float(improvement);
  EA_status.reproduction_rate = repro_per_sec;
#else  // JAVA_GUI
  if (sim != 0) sim->printMessage(text);
#endif // JAVA_GUI

  time_of_last_call = clock();
} // end CInTrance::PrintReproduction()


//-----------------------------------------------------------------------------
void PrintReproduction2(     CIndividual* ind,
                       const double improvement,
                       const int n_epc_repr,         // # of epoch reproductions
                       const int n_abs_repr,      // absolute # of reproductions
                       void* data)
{
  using namespace approx;

  // get system time and calculate time since last call
  clock_t time_of_new_call = clock();
  float diff_time        = (float(time_of_new_call)-float(time_of_last_call))/float(CLOCKS_PER_SEC);

  // estimate the reproductions per second ratio 
  float repro_per_sec = 0.0;
  if (diff_time > 0.0) repro_per_sec = ((n_abs_repr-n_last_reproductions)/diff_time);
  if (repro_per_sec < 0.0) repro_per_sec = 0.0;

  // save absolute # of reproductions
  n_last_reproductions = n_abs_repr;

  CSimulation* sim = 0;
  try
  { 
    sim = static_cast<CSimulation*>(data);
  }
  catch (std::bad_cast&)
  {
    Warning("No valid CSimulation pointer in PrintReproduction2!");
  }

  CApproxFitness* f(0);
  ostringstream mess("");
  if (ind != 0)
  {
    f = dynamic_cast<CApproxFitness*>(ind->getFitnessObject());
    if (f != 0)
    {
      const unsigned int w(12), p(9);
      mess << setw(w) << showpos << left << setprecision(p);
      mess << f->getMixedFitness() << setw(1) <<"(";
      mess << setw(w) 
           << setfill('0') << showpos << showpoint << left << setprecision(p)<<improvement
           << ") ";
      mess << ", cl1:" << setw(w/2) << setfill('0') << showpos << left << setprecision(p/2) 
           << f->getFitness(CLOCK_SEG);
      mess << ", co1:" << setw(w/2) << setfill('0') << showpos << left << setprecision(p/2) 
           << f->getFitness(CONE_SEG);
      mess << ", seg:" << setw(w/2) << setfill('0') << showpos << left << setprecision(p/2) 
           << f->getFitness(SEGMENTS);
      mess << ", cl2:" << setw(w/2) << setfill('0') << showpos << left << setprecision(p/2) 
           << f->getFitness(CLOCK_SC);
      mess << ", co2:" << setw(w/2) << setfill('0') << showpos << left << setprecision(p/2) 
           << f->getFitness(CONE_SC);
      mess << ", sc:" << setw(w/2) << showpos << left << setprecision(p/2) 
           << f->getFitness(SIGN_CHANGES);
      mess << ", opt:" << setw(w) << showpos << left << setprecision(p) 
           << f->getFitness(OPT_FITNESS);
      mess << resetiosflags(ios::fixed|ios::showpos);
      mess << setw(1)<<setfill(' ')<<" (" 
           << setfill('0') << setw(6) << right <<n_epc_repr<<"/"<< setw(6)<<n_abs_repr;
      mess << ", " << setfill(' ') 
           << right << showpoint<< setw(6) << setprecision(2) 
           << fixed << repro_per_sec << " rep/sec)"<< endl;
    }
  }
  else 
  {
    mess << "Pointer of individual points to zero!\n";
  }

#if defined JAVA_GUI
  string text_stl               = mess.str();

  Trajectory ref_traj              = 0;
  Trajectory traj                  = 0;
  Trajectory target_traj           = 0;
  Trajectory approx_traj           = 0;
  DataArray  clock_angle           = 0;
  DataArray  cone_angle            = 0;
  DataArray  clock_angle_approx    = 0;
  DataArray  cone_angle_approx     = 0;
  DataArray  thrust                = 0;
  SCTraj     sc_traj;

  if (sim!=0 && f!=0)
  {
    try 
    {      
      sim->evaluate(*ind, false, true);
    }
    catch (std::invalid_argument)
    {
      Warning("Desired fitness could not be found in last \
 evaluated approximations!\n");
    }
    approx_traj        = sim->getLastTrajectoryApproxData();
    clock_angle_approx = sim->getLastTrajectoryClockAngleApproximation();
    cone_angle_approx  = sim->getLastTrajectoryConeAngleApproximation();
    ref_traj           = sim->getReferenceTrajectoryData();
    target_traj        = sim->getReferenceTrajectoryTargetData();
    clock_angle        = sim->getReferenceTrajectoryClockAngleData();
    cone_angle         = sim->getReferenceTrajectoryConeAngleData();
    thrust             = sim->getReferenceTrajectoryThrustData();
  }

  // All the above standing calls to "get..."-methods of CSimulation class may
  // return zero but there's no special error handling code needed
  // because the data transfer routines of "InTranceLib.cpp" ignore data array
  // pointers set to zero. So in case of errors or no data at all there
  // is no serious threat to program execution.

  TransferDataItem item(ref_traj
                       ,approx_traj
                       ,target_traj
                       ,clock_angle
                       ,clock_angle_approx
                       ,cone_angle
                       ,cone_angle_approx
                       ,thrust
                       ,0
                       ,text_stl
                       ,repro_per_sec);


  dataBuffer->put(item);   
  sim->writeMessageToReportFile(mess);
#else // #if defined JAVA_GUI
  if (sim != 0) sim->printMessage(mess);
#endif // JAVA_GUI

  time_of_last_call = clock();
} // void PrintReproduction2(...)


//-----------------------------------------------------------------------------
void CInTrance::reset(void) 
{
  m_initialized = false;
  input_inp.reset();
  input_ncp.reset();
  input_eap.reset();
  input_scp.reset();
  input_sim.reset();
  SafeDelete(m_solsys);
  m_rep.erase();
}


//-----------------------------------------------------------------------------
void CInTrance::initializeFitnessFunctionType2(const int fft2)
throw (std::invalid_argument&)
{
  const string method_name("void CInTrance::initializeFitnessFunctionType2(const int fft3):");
  const string err_msg1("Invalid InTrance command!\n");
  const string warn_msg1("Not applicable for given InTrance command!\n");
  const string warn_msg2("Invalid fitness function type!\n");

  switch (input_inp.what_to_do )
  {
  case OPTIMIZE_NC:
  case EVALUATE_NC:
  case EVALUATE_CONTROL:
  case LSL_MINMAX:
  case LSL_OPT:
  case ANGLES:
  case OPTIMIZE_APPROX:
  case EVALUATE_APPROX:
    if ((fft2<int(approx::NO_FFT1)) || (fft2>=approx::FFT2_END)) 
    {
      Warning(method_name + warn_msg2);
      jt2 = approx::NO_FFT1;
    }
    else jt2 = static_cast<approx::FitnessFunctionType2>(fft2);
    break;
  default:
    throw std::invalid_argument(method_name + err_msg1);
    break;
  }
}


//-----------------------------------------------------------------------------
void CInTrance::initializeFitnessFunctionType3(const int fft3)
throw (std::invalid_argument&)
{
  const string method_name("void CInTrance::initializeFitnessFunctionType3(const int fft3):");
  const string err_msg1("Invalid InTrance command!\n");
  const string warn_msg1("Not applicable for given InTrance command!\n");
  const string warn_msg2("Invalid fitness function type!\n");

  switch (input_inp.what_to_do )
  {
  case OPTIMIZE_NC:
  case EVALUATE_NC:
  case EVALUATE_CONTROL:
  case LSL_MINMAX:
  case LSL_OPT:
  case ANGLES:
  case OPTIMIZE_APPROX:
  case EVALUATE_APPROX:
    if ((fft3<int(approx::SEGM_SIGNCHANGE1)) || (fft3>=approx::FFT3_END)) 
    {
      Warning(method_name + warn_msg2);
      jt3 = approx::SIGNCHANGE;
    }
    else jt3 = fft3;
    break;
  default:
    throw std::invalid_argument(method_name + err_msg1);
    break;
  }
}


//-----------------------------------------------------------------------------
void CInTrance::initializeFitnessFunctionType(void)
throw(std::invalid_argument&)
{
  const std::string method_name = "CInTrance::\
initializeFitnessFunctionType(void):";

  // read fitness function type 1
  string fft(input_eap.getFitnessFunctionType());

  jt=N_FFT;
  int i;
  for (i = 0;i<N_FFT;++i)
  {
    if (fft == FitnessFunctionTypeDesc[i])
    {
      jt = FitnessFunctionType(i);
      i=N_FFT;
    }
  }
  if (jt==N_FFT) throw std::invalid_argument(method_name + 
    " unknown fitness function type " + fft);

  // read fitness function type 2
  int t = input_eap.getFitnessFunctionType2();
  initializeFitnessFunctionType2(t);

  // read fitness function type 3
  t = input_eap.getFitnessFunctionType3();
  initializeFitnessFunctionType3(t);
}

//-----------------------------------------------------------------------------
void CInTrance::readConfigParams(
                                     ParamList  & parameters,
                                     const std::string& filename,
                                     const std::string& proc_name
                                     )
                                     throw (std::invalid_argument&, std::logic_error&)
{
  parameters.clear();
  ifstream config_file(filename.c_str());
  if (config_file.is_open()) 
  {
    try 
    { 
      parameters=GetParameters(config_file, COMMENT_SEPERATOR); 
    }
    catch (syntax_error_in_config_file& exc) 
    {
      throw logic_error(proc_name+ string(exc.what()));
    }
    config_file.close();
  } else throw invalid_argument(proc_name
    +string("could not open ea input file (")
    +input_inp.evolutionary_algorithm_parameter_file+string(")\n"));
}


//-----------------------------------------------------------------------------
void CInTrance::readSimConfigParams(void) 
  throw (std::invalid_argument&, std::logic_error&)
{
  const string proc_name("void CInTrance::readSimConfigParams(void) ");
  ParamList p;
  readConfigParams(p, input_inp.simulation_parameter_file, proc_name);
  input_sim = CInputSimulation(p, m_solsys);
}

//-----------------------------------------------------------------------------
void CInTrance::readScConfigParams(void) 
  throw (std::invalid_argument&, std::logic_error&)
{
  const string proc_name("void CInTrance::readScConfigParams(void) ");
  ParamList p;
  readConfigParams(p, input_inp.spacecraft_parameter_file, proc_name);
  input_scp = CInputSC(p);
}

//-----------------------------------------------------------------------------
void CInTrance::readEaConfigParams(void)
  throw (std::invalid_argument&, std::logic_error&)
{
  const string proc_name("void CInTrance::readEaConfigParams(void) ");
  ParamList p;
  readConfigParams(p, input_inp.evolutionary_algorithm_parameter_file,
    proc_name);
  input_eap = CInputEA(p);
}

//-----------------------------------------------------------------------------
void CInTrance::readNcConfigParams(const bool extra_params)
  throw (std::invalid_argument&, std::logic_error&)
{
  const string proc_name("void CInTrance::readNcConfigParams(void) ");
  ParamList p;
  readConfigParams(p, input_inp.nc_topology_file, proc_name);
  input_ncp = CInputNC(p,extra_params);
}

//-----------------------------------------------------------------------------
void CInTrance::initSolarSystem(Global::ParamList& p)
{
  const std::string method_name = 
    "void CInTrance::initSolarSystem(Global::ParamList& p)";
  const std::string err_msg_01 = "could not allocate solar system object!\n";

  // allocate solar system object by reading DE405 ephemerides parameter file 
  // and creating a CSolarSystem object
  std::string fn;
  InputFilenameParameter("DE405_DATA_FILE", fn, p, method_name
    , TRAILING_CHARS, "", false, false);

  try
  {
    SafeDelete(m_solsys);
    if (fn.empty() == false)
    {
      removeChar(fn, '\"');
      m_solsys = new CSolarSystem(fn);
    }
    else
    {
      SafeDelete(m_solsys);
      m_solsys = new CSolarSystem;
    }
  }
  catch (std::bad_alloc)
  {
    throw invalid_argument(method_name+err_msg_01);
  }
} // void initSolarSystem(Global::ParamList& p)

//-----------------------------------------------------------------------------
void CInTrance::initAsteroidDataFile(Global::ParamList& p)
{
  const std::string method_name = 
    "void CInTrance::initAsteroidDataFile(const Global::ParamList& p)";
  std::string fn;
  ostringstream t("");
  int idx=0;
  if (m_solsys!=0)
  {
    try
    {
      do
      {
        t.str("");t << (++idx);
        InputFilenameParameter("ASTEROID_DATA_FILE"+t.str(), fn
          , p, method_name, TRAILING_CHARS, "", true, false);
        removeChar(fn, '\"');
        m_solsys->setAsteroidDataFile(fn);
      } while (1==1);
    }
    catch (std::invalid_argument) {}
  }
} // void initAsteroidDataFile(const Global::ParamList& p)

//-----------------------------------------------------------------------------
void CInTrance::initCometDataFile(Global::ParamList& p)
{
  const std::string method_name = 
    "void CInTrance::initAsteroidDataFile(Global::ParamList& p)";
  std::string fn;
  try
  {
    InputFilenameParameter("COMET_DATA_FILE", fn, p, method_name
      , TRAILING_CHARS,"",true, false);
    removeChar(fn, '\"');
    if (m_solsys!=0)
    {
      try
      {
        m_solsys->setCometDataFile(fn);
      }
      catch (invalid_argument& exc) 
      {
        Output(method_name+std::string(exc.what()));
      }
    } // if (m_solsys!=0)
  } // try
  catch (std::invalid_argument&) 
  {
    Output("\nnot using comet ephemeris file\n");
  }

} // void initAsteroidDataFile(Global::ParamList& p)

//-----------------------------------------------------------------------------
void CInTrance::initialize(const std::string& inp_filename)
throw (std::logic_error&)
{
  std::string PROC_NAME("CInTrance::initialize(...):");
  bool error=false;

  // check name of given input file for empty string and throw exception if so
  if (inp_filename.empty() == true) throw invalid_argument(PROC_NAME+
    string("empty input file name string!"));

  time_of_last_call = 0;      // init static time variable
  n_last_reproductions = 0;   // also the number of last reproductions
  ResetReproCtr() ;           // reset counter of absolute # of reproductions

  // get parameters from config file
  ifstream config_file;
  config_file.open(inp_filename.c_str());
  ParamList parameters;
  try 
  { 
    parameters=GetParameters(config_file, COMMENT_SEPERATOR); 
  }
  catch (syntax_error_in_config_file&exc) 
  {
    error=true;
    throw logic_error(PROC_NAME+string(exc.what()));
  }
  config_file.close();

  std::string path = inp_filename;

  StripPathFromFilename(path);        // strip leading path from input file name

  if (!error) initSolarSystem(parameters);// create internal solar system object
  if (m_solsys == 0) error = true;

  // asteroid orbital elements data file
  if (!error) initAsteroidDataFile(parameters);

  // comet orbital elements data file
  if (!error) initCometDataFile(parameters);

  // read InTrance input data file
  //-----------------------------------
  if (error == false)
  {
    try 
    {
      if (path.empty()==false) parameters.push_back(Param("PATH",path));
      if (!error) input_inp.readParams(parameters);
    }
    catch (invalid_argument& exc) 
    {
      throw invalid_argument(PROC_NAME + string(exc.what()));
    }
  }

  m_rep=input_inp.report_file;    // set InTrance report file name

  // initialize random number generator
  if (error == false)
  {
    // with given seed value
    if (input_inp.random_seed_provided==true) srand(input_inp.random_seed);
    else 
    {
      // or with system time, if no seed value is given
      input_inp.setRandomSeed(unsigned int(time(0)));
      srand(input_inp.getRandomSeed());
    }

    if (m_rep.size()!=0)
    {
      ofstream report(m_rep.c_str(), ios::out);
      if (report.is_open() == true) flush(report<<"InTrance 1.7 Report File\n");
    }
  }

  // initialize optimization depending on given task
  try 
  {
    switch (input_inp.what_to_do)
    {
    case OPTIMIZE_NC:
    case EVALUATE_NC:initializeNcOptimization();break;
    case EVALUATE_CONTROL:initializeControlEvaluation();break;
    case LSL_OPT:
    case LSL_MINMAX: initializeLslOptimization();break;
    case ANGLES:initializeLslOptimization();break;
    case OPTIMIZE_APPROX:
    case EVALUATE_APPROX: initializeApproximationOptimization();break;
    default: throw invalid_argument("\n Don't know what to do!\n");
    }
  }
  catch (std::invalid_argument& exc)
  {
    throw logic_error(PROC_NAME+string(exc.what()));
  }

  if (error == false) m_initialized = true;

} // initialize()


//-----------------------------------------------------------------------------
void CInTrance::initializeNcOptimization(void) throw(std::invalid_argument&)
{
  // the InTrance input object defines, what has to be done
  // ============================================
  if ((input_inp.what_to_do == OPTIMIZE_NC) ||
    (input_inp.what_to_do == EVALUATE_NC)) 
  {
    readNcConfigParams();
    input_ncp.setExtraParameters(ADDITIONAL_CHROMOSOME_DATA);

    readScConfigParams();
    readSimConfigParams();
    readEaConfigParams();

    // set default values for possibly not given necessary parameters

    if (input_eap.getHypercubeUpperLimit()==DEFAULT_HYPERCUBE_UPPER_LIMIT)
      input_eap.setHypercubeUpperLimit(LIMIT_HYPERSPACE_LARGE);

    // final boundary condition met limit
    if (input_eap.isDefaultConvergenceLimit(FBC_MET)==true) 
      input_eap.setConvergenceLimit(FBC_MET,LIMIT_FINAL_BOUNDARY_CONDITION_MET);

    // population convergence limit 1 (FBC met) 
    if (input_eap.isDefaultConvergenceLimit(POP_CL_FBC_MET)==true) 
      input_eap.setConvergenceLimit(POP_CL_FBC_MET,EPOCH_CONVERGENCE_LIMIT_FINAL_BOUNDARY_CONDITION_MET);

    // population convergence limit 2 (FBC not met) 
    if (input_eap.isDefaultConvergenceLimit(POP_CL_FBC_NMET)==true) 
      input_eap.setConvergenceLimit(POP_CL_FBC_NMET,EPOCH_CONVERGENCE_LIMIT_FINAL_BOUNDARY_CONDITION_NOT_MET);

    // evol. alg. convergence limit 1 (FBC met) 
    if (input_eap.isDefaultConvergenceLimit(EA_CL_FBC_MET)==true) 
      input_eap.setConvergenceLimit(EA_CL_FBC_MET,EA_CONVERGENCE_LIMIT_FINAL_BOUNDARY_CONDITION_MET);

    // evol. alg. convergence limit 2 (FBC not met) 
    if (input_eap.isDefaultConvergenceLimit(EA_CL_FBC_NMET)==true) 
      input_eap.setConvergenceLimit(EA_CL_FBC_NMET,EA_CONVERGENCE_LIMIT_FINAL_BOUNDARY_CONDITION_NOT_MET);

    initializeFitnessFunctionType(); // set fitness function type

  } // end "OPTIMIZE_NC" || "EVALUATE_NC"
}

//-----------------------------------------------------------------------------
void CInTrance::initializeControlEvaluation(void) throw(std::invalid_argument&)
{
  readScConfigParams();
  readSimConfigParams();
}

//-----------------------------------------------------------------------------
void CInTrance::initializeLslOptimization(void) throw(std::invalid_argument&)
{
  readScConfigParams();
  readSimConfigParams();
}

//-----------------------------------------------------------------------------
void CInTrance::initializeApproximationOptimization(void) 
throw (std::invalid_argument&)
{
  ParamList parameters;
  readScConfigParams();
  readNcConfigParams(false);
  readSimConfigParams();
  readEaConfigParams();
  initializeFitnessFunctionType(); // set fitness function type

  // set default values for possibly not given necessary parameters
  // final boundary condition met limit
  if (input_eap.isDefaultConvergenceLimit(FBC_MET)==true) 
    input_eap.setConvergenceLimit(FBC_MET,APPROX_FBC_MET);

  // population convergence limit 1 (FBC met) 
  if (input_eap.isDefaultConvergenceLimit(POP_CL_FBC_MET)==true) 
    input_eap.setConvergenceLimit(POP_CL_FBC_MET,APPROX_POP_CONV_FBC_MET);

  // population convergence limit 2 (FBC not met) 
  if (input_eap.isDefaultConvergenceLimit(POP_CL_FBC_NMET)==true) 
    input_eap.setConvergenceLimit(POP_CL_FBC_NMET,APPROX_POP_CONV_FBC_NOT_MET);

  // evol. alg. convergence limit 1 (FBC met) 
  if (input_eap.isDefaultConvergenceLimit(EA_CL_FBC_MET)==true) 
    input_eap.setConvergenceLimit(EA_CL_FBC_MET,APPROX_EPC_CONV_FBC_MET);

  // evol. alg. convergence limit 2 (FBC not met) 
  if (input_eap.isDefaultConvergenceLimit(EA_CL_FBC_NMET)==true) 
    input_eap.setConvergenceLimit(EA_CL_FBC_NMET,APPROX_EPC_CONV_FBC_NOT_MET);
}


//-----------------------------------------------------------------------------
bool IsBatchFile(const std::string& filename)
{
  if (filename.rfind(".ibf") == filename.length()-4) return true;
  else return false;
} // bool CIntrance::isBatchFile(const std::string& filename)


//-----------------------------------------------------------------------------
void CInTrance::insertConfigFilename(const std::string& filename)
{
  const std::string method_name = 
    "void CInTrance::insertConfigFilename(const std::string& filename)";

  const std::string err_msg1 = method_name +
    std::string(": couldn't find \"") + filename + std::string("\"\n");

  std::ifstream test_file(filename.c_str());

  // throw exception if file doesn't exist
  if (test_file.is_open() == false) throw std::runtime_error(err_msg1);

  // determine file type (single input file or batch file ?)
  if ( IsBatchFile(filename) == true )
  {
    std::string path(filename);
    const char* seperator = "/\\";       // use both types (WINDOWS,UNIX)
    std::string::size_type pos = path.find_last_of(seperator);
    if (pos != std::string::npos) path.erase(pos+1, path.length()-pos-1);
    else path.erase();

    std::string file_content;
    while (test_file.eof() != true)
    {
      file_content.erase();              // clear string for actual file name
      getline(test_file, file_content); // get file name
      removeComments(file_content, COMMENT_SEPERATOR);  // remove comments
      // and remaining characters at the end
      removeTrailingChars(file_content, TRAILING_CHARS); 
      if (file_content.length() > 0)
      {
        // remove leading spaces
        while (file_content[0] == ' ') file_content.erase(0, 1);
        if (file_content.length()>2)
        {
          int i = 0;
          if (file_content[0]=='"') ++i;
          if ( (file_content[i]=='/') || (file_content[i+1]==':') )
            insertConfigFilename(file_content);
          else insertConfigFilename(path+file_content);
        }
      }
    }
  } // if (isBatchFile(filename) == true )
  else // single input file
  {
    m_input_files.push_back(filename);
  } // else

} // void CInTrance::insertConfigFilename(const std::string& filename)

//-----------------------------------------------------------------------------
CInTrance::CInTrance(const char* inp_filename)
: m_initialized(false)
 ,m_is_running(false)
 ,m_solsys(0)
 ,m_ea_status(0)
#if defined COMMONCPP
 ,m_gui_output(0)
 ,m_clock_angle_history(false)
 ,m_clock_angle_history_approximation(false)
 ,m_clock_angle_approx_poly_degree(0)
 ,m_clock_angle_approx_seg_size(0)
 ,m_cone_angle_history(false)
 ,m_cone_angle_history_approximation(false)
 ,m_cone_angle_approx_poly_degree(0)
 ,m_cone_angle_approx_seg_size(0)
 ,m_thrust_history(false)
#endif //COMMONCPP
 ,input_ncp()
 ,input_scp()
{
  m_break = false;
  sim = 0;
  ea = 0;
  input_inp = CInputInTrance();
  input_sim = CInputSimulation();
  input_eap = CInputEA();

#if defined COMMONCPP
  setOutput(&java_output);
#endif

  insertConfigFilename(inp_filename);

}; // end CInTrance::CInTrance(const char* inp_filename)

//------------------------------------------------------------------------------
void CInTrance::stop(void) 
{
  m_break = true;
}

//------------------------------------------------------------------------------
void CInTrance::evaluateControl(void)
{
  input_scp.printSettings();  // print s/c settings to inform user
  input_sim.printSettings();  // print simulation settings to inform user

  // create simulation object from simulation parameter input file
  auto_ptr<CSimulation> sim(new CSimulation(input_sim
    , input_scp, NULL, input_inp.report_file));

  // create spacecraft attitude history
  CSCAttitudeHistory att_history(input_inp.simulation_input_file_ang.c_str());

  // evaluate attitude history
  if (sim.get() != 0) 
  {
    sim->evaluate(att_history);

    // write simulation output files
    sim->writeOutputFiles(
      input_inp.simulation_output_file_csv
      ,input_inp.simulation_output_file_dat
      ,input_inp.simulation_output_file_gesop
      ,input_inp.simulation_output_file_wrl
      ,input_inp.simulation_output_file_ang
      );
  } // if (sim.get() != 0) 
} // void CInTrance::evaluateControl(void)


//-----------------------------------------------------------------------------
void CInTrance::run(void)
{
  const char* bridge_mode = std::getenv("INTRANCE_BRIDGE_MODE");
  const char* observer_config = std::getenv("INTRANCE_ENC_OBSERVER_CONFIG");
  const bool bridge_requested = bridge_mode && *bridge_mode;
  const bool controlled_run = bridge_requested || (observer_config && *observer_config);
  m_is_running = true;
  m_break = false;
  std::ostringstream mess;

  while( (m_input_files.empty()==false) && (abort()==false) ) 
  {
    CInTrance::reset();                                // reset static variables
    time_of_last_call = clock();
    resetTimer();
    
    try                          // try initializing of optimization environment
    {
      initialize(m_input_files.front());
      setOutputLogFilename(input_inp.report_file);
      input_inp.printSettings();
    }
    catch (logic_error& exc)  // error handling
    {
      if (bridge_requested) NativeBridgeProbe::mark_failed();
      // A failed controlled stage must not be reported as a completed flight.
      if (controlled_run) { m_is_running = false; throw; }
      Output(string(exc.what()));
      input_inp = CInputInTrance();
    }
    catch (std::exception&)
    {
      if (bridge_requested) NativeBridgeProbe::mark_failed();
      m_is_running = false;
      throw; // Other initialisation failures already escaped the legacy logic_error handler.
    }
    m_input_files.pop_front();  // remove processed input file from queue

    try // try evaluation depending on input task
    {
      const char* observer_request=std::getenv("INTRANCE_ENC_OBSERVER_CONFIG");
      if(observer_request && std::string(observer_request)!="" && input_inp.what_to_do!=OPTIMIZE_NC)
        throw std::runtime_error("ENC_OBSERVER_REFUSE_NON_SEARCH_MODE");
      switch (input_inp.what_to_do)
      {
      case OPTIMIZE_NC      :optimizeNc();break;
      case EVALUATE_NC      :optimizeNc();break;
      case LSL_MINMAX       :optimizeNc();break;
      case LSL_OPT          :optimizeNc();break;
      case ANGLES           :optimizeNc();break;
      case EVALUATE_CONTROL :evaluateControl();break;
      case OPTIMIZE_APPROX  :optimizeApproximation();break;
      case EVALUATE_APPROX  :optimizeApproximation();break;
      default:break;
      } // switch (input_inp.what_to_do)
    } // try
    // error handling
    catch (std::exception& exc)
    {
      if (bridge_requested) NativeBridgeProbe::mark_failed();
      if (controlled_run) { m_is_running = false; throw; }
      Output(exc.what());
    }
  
    stopTimer();    
    mess.str("");                                       //print calculation time
    mess <<"total running time: "<<getOptimizationTimespan() << " seconds\n\n";
    Output(mess.str());

    if (m_input_files.size()>0)      // are there still input files to process ?
    {
      Output("\ninput files left to process:\n");         // do some info output
      std::list<std::string>::iterator i;
      for (i = m_input_files.begin(); i != m_input_files.end(); ++i) 
      {
        Output(*i+"\n");
      }
      Output("\n");
    } // if (m_input_files.size()>0)      
    else { 
        //_CrtDumpMemoryLeaks();
        Output("\nInTrance run finished!\n"); };
  } // while(m_input_files.empty != true

  m_is_running = false;
} // end of CInTrance::run(void)


//-----------------------------------------------------------------------------
void CInTrance::setEAStatusObject(EAStatusRec* status)
{
  assert(status);
  m_ea_status=status;
} // void CInTrance::setEAStatusObject(EAStatusRec* status)


//-----------------------------------------------------------------------------
void CInTrance::optimizeNc(void) throw(std::exception&)
{
  const std::string method_name("CInTrance::optimizeNc(void)");
  std::auto_ptr<CNC> nc;

  // the InTrance input object defines, what has to be done
  if ( (input_inp.what_to_do==OPTIMIZE_NC) 
    || (input_inp.what_to_do==EVALUATE_NC))
  {
    if (abort() == false)                       // create neurocontroller object
    {
      // first, create neurocontroller
      try 
      {
        input_ncp.printSettings();
        nc = std::auto_ptr<CNC>(new CNC(input_ncp)); 
      }
      catch (runtime_error& exc)
      {
        Output(method_name+string("\n")+string(exc.what()));
      }
    } // if (abort() == false)


    // create simulation object    
    if ((abort() == false) && (nc.get()!=0) )
    { // try creation of object
      try 
      {
        input_scp.printSettings();
        input_sim.printSettings();

#ifdef CONSEP
        sim = new CCONSEPSim(input_sim,input_scp, nc.get(),input_inp.report_file);
#else
        sim = 
          new CSimulation(input_sim,input_scp, nc.get(),input_inp.report_file);
#endif // CONSEP

        if (sim==0) throw runtime_error(" (simulation object = 0");
      }
      // handling of possibly occured error
      catch (runtime_error& exc)
      {
        Output(method_name + string(" : Could not create simulation object!\n")+
          string(exc.what()));
      }

      // create evolutionary algorithm object 
      if ((abort() == false) && (sim!=0))
      {
        simulation = sim;
        sim->setAbortBoolPointer(&m_break);
#if defined COMMONCPP
        sim->setStatusObject(&SIM_status);

        if (m_clock_angle_history == true) 
          sim->activateClockAngleHistory();

        if (m_clock_angle_history_approximation == true) 
          sim->activateClockAngleApproximation();

        sim->setClockAngleApproximationParameters(
          m_clock_angle_approx_seg_size,
          m_clock_angle_approx_poly_degree);

        if (m_cone_angle_history == true) 
          sim->activateConeAngleHistory();

        if (m_cone_angle_history_approximation == true) 
          sim->activateConeAngleApproximation();

        sim->setConeAngleApproximationParameters(
          m_cone_angle_approx_seg_size,
          m_cone_angle_approx_poly_degree);

        if (m_thrust_history) sim->activateThrustHistory();
#endif //COMMONCPP

        input_eap.setChromosomeLength(nc->getNumberOfParameters()
          + nc->getNumberOfExtraParameters());
        input_eap.printSettings();
        SafeDelete(ea);
        ea = new CEA( input_eap,nc.get(),sim
          ,m_rep,PrintReproduction,CalcAccuracyFactor,sim);

        if ( (abort()==false) && (ea!=0) ) 
        {
          ea->setStopBoolPointer(&m_break);
          ea->setStatusObject(m_ea_status);
          ea->setKeyboardHandler(&(handleKeyboardInputs));
          ea->setEventHandler(&(handleEaEvents));
          // The count-bank observer is opt-in and only valid for cold-start NC search.
          const char* observer_config=std::getenv("INTRANCE_ENC_OBSERVER_CONFIG");
          if(observer_config && std::string(observer_config)!="")
          {
            const char* smoke=std::getenv("INTRANCE_ENC_ORDERLY_STOP_SMOKE");
            if(input_inp.what_to_do!=OPTIMIZE_NC || NativeBridgeProbe::configured() ||
               (smoke && std::string(smoke)!="") ||
               !input_inp.evolutionary_algorithm_warmstart_file.empty() ||
               !input_inp.random_seed_provided || enc_smoke_accounting ||
               sim->getMinStartTime()!=52654.5 || sim->getMaxStartTime()!=52654.5 ||
               sim->getInitialStateType()!=BODY_STATE || sim->getInitialBody()!=EARTH ||
               sim->getTargetStateType()!=BODY_STATE || sim->getTargetBodyId()!=MERCURY ||
               sim->getTrajectoryType()!=RENDEZVOUS ||
               sim->getMinLaunchVelocity()!=0.0 || sim->getMaxLaunchVelocity()!=0.0 ||
               sim->getIntegrationInterval()!=600.0 || sim->useDAC() || sim->getStatusDSSC() ||
               sim->getDynamicIntegrationIntervalReduction() ||
               !sim->getCurrentSC().isSolarSail() || !input_scp.ideal_or_real_sail)
              throw std::runtime_error("ENC_OBSERVER_REFUSE_MODE");
#if defined(COMMONCPP)
            throw std::runtime_error("ENC_OBSERVER_REFUSE_COMMONCPP");
#endif
            ea->enableObserver(observer_config,input_inp.random_seed);
            Output("ENC_OBSERVER_ENABLED count_bank_only; capture is not an optimal handover policy\n");
          }
          // ------------------------------------------
          if (input_inp.what_to_do == OPTIMIZE_NC)
          {
#if defined(INTRANCE_ENC_ORDERLY_STOP_SMOKE)
            // Explicit opt-in in this isolated test build, never in EVALUATE_NC.
            const char* smoke_token=std::getenv("INTRANCE_ENC_ORDERLY_STOP_SMOKE");
            if (smoke_token!=0 && std::string(smoke_token)!="")
            {
              if (std::string(smoke_token)!="42_SEARCH_CALLS")
                throw std::runtime_error("ENC_SMOKE_REFUSE_TOKEN");
#if defined(COMMONCPP)
              throw std::runtime_error("ENC_SMOKE_REFUSE_COMMONCPP");
#endif
              if (enc_smoke_accounting || NativeBridgeProbe::configured() ||
                  !input_inp.evolutionary_algorithm_warmstart_file.empty() ||
                  !input_inp.random_seed_provided || input_inp.random_seed!=1 ||
                  ea->getPopulationSize()!=40 || ea->getMaxSSSEpochsCount()!=1 ||
                  sim->getMinStartTime()!=52654.5 || sim->getMaxStartTime()!=52654.5 ||
                  sim->getInitialStateType()!=BODY_STATE || sim->getInitialBody()!=EARTH ||
                  sim->getTargetStateType()!=BODY_STATE || sim->getTargetBodyId()!=MERCURY ||
                  sim->getTrajectoryType()!=RENDEZVOUS ||
                  sim->getMinLaunchVelocity()!=0.0 || sim->getMaxLaunchVelocity()!=0.0 ||
                  sim->getIntegrationInterval()!=600.0 || sim->useDAC() || sim->getStatusDSSC() ||
                  sim->getDynamicIntegrationIntervalReduction() ||
                  !sim->getCurrentSC().isSolarSail() || !input_scp.ideal_or_real_sail)
                throw std::runtime_error("ENC_SMOKE_REFUSE_CHANGED_FIXTURE_OR_MODE");
              ea->enablePreparationSmoke();
              enc_smoke_accounting=true;
              Output("ENC_SMOKE_ENABLED search_limit=42 population=40 scan_epochs=1 seed=1\n");
            }
#endif

            // try to run ev. alg.
            if ((abort()==false) && (ea!=0) ) // should abort or invalid ea object?
            {
              ea->Optimize(input_inp.evolutionary_algorithm_warmstart_file.c_str());
              if (ea->preparationSmokeEnabled())
              {
                ostringstream smoke;
                smoke << "ENC_SMOKE_SEARCH_RETURN search_calls=" << ea->preparationSearchCalls()
                  << " completed_pairs=" << ea->preparationPairs()
                  << " preparation_stop=" << ea->preparationStopRequested()
                  << " global_abort=" << m_break << " search_monitors=" << enc_smoke_search_monitors
                  << " stop_reason=" << (m_break ? "user_or_simulation_abort" :
                    (ea->preparationStopRequested() ? "fixed_preparation_limit" : "other")) << "\n";
                Output(smoke.str());
                if (!m_break && (!ea->preparationStopRequested() ||
                    ea->preparationSearchCalls()!=42 || ea->preparationPairs()!=1))
                  throw std::runtime_error("ENC_SMOKE_REFUSE_INCOMPLETE_NORMAL_RETURN");
              }
              // save best found solution to file
              ea->SaveBestSolution(
                input_inp.evolutionary_algorithm_best_solution_file.c_str());
              if (enc_smoke_accounting)
                Output("ENC_SMOKE_AUTHORITATIVE_SAVE_RETURNED\n");
            }
          } // input_inp.what_to_do == OPTIMIZE_NC

          // ------------------------------------------
          switch (input_inp.what_to_do)
          {
          case OPTIMIZE_NC:
          case EVALUATE_NC:
            if (ea != 0) 
            {
              CIndividual solution(ea
                , input_inp.evolutionary_algorithm_best_solution_file.c_str());
              if (enc_smoke_accounting)
              {
                enc_smoke_final_report=true;
                ostringstream smoke;
                smoke << "ENC_SMOKE_RELOAD_RETURNED parameters=" << solution.getLength() << "\n";
                Output(smoke.str());
              }
              // final solution evaluation
              sim->evaluate(solution, true,true);
              if (enc_smoke_accounting) ++enc_smoke_final_evaluations;
              // In the normal non-COMMONCPP build PrintReproduction performs
              // a second monitor evaluation and replaces the stored trajectory
              // subsequently written by writeOutputFiles. Capture/replay only
              // that whole reporting pass; leave the prior NN evaluation intact.
              if (NativeBridgeProbe::configured())
              {
                NativeBridgeProbe::require(input_inp.what_to_do == EVALUATE_NC,
                  "diagnostic supports saved-controller EVALUATE_NC only");
#if defined(COMMONCPP)
                NativeBridgeProbe::require(false,
                  "reporting-pass diagnostic does not support COMMONCPP builds");
#endif
              }
              {
                NativeBridgeProbe::ReportingPass bridge_pass;
                PrintReproduction(&solution, 1.0, 1, 1, getCurrentSimObj());
                bridge_pass.complete();
                if (enc_smoke_accounting)
                {
                  ostringstream smoke;
                  smoke << "ENC_SMOKE_REPORT_RETURNED search_monitors=" << enc_smoke_search_monitors
                    << " final_evaluations=" << enc_smoke_final_evaluations
                    << " final_monitors=" << enc_smoke_final_monitors
                    << " global_abort=" << m_break
                    << " best_flight_days=" << sim->getBestFlightTime()
                    << " native_mixed_fitness=" << solution.getFitnessObject()->getMixedFitness() << "\n";
                  Output(smoke.str());
                }
              }
            }        
            if (sim != 0) sim->writeOutputFiles(  // write simulation output files
              input_inp.simulation_output_file_csv
              ,input_inp.simulation_output_file_dat
              ,input_inp.simulation_output_file_gesop
              ,input_inp.simulation_output_file_wrl
              ,input_inp.simulation_output_file_ang);
            if (enc_smoke_accounting)
            {
              if (m_break || !ea->preparationStopRequested() ||
                  ea->preparationSearchCalls()!=42 || ea->preparationPairs()!=1 ||
                  enc_smoke_final_evaluations!=1 || enc_smoke_final_monitors!=1 ||
                  enc_smoke_search_monitors<1 || enc_smoke_search_monitors>2)
                throw std::runtime_error("ENC_SMOKE_REFUSE_INCOMPLETE_CONTROL_FLOW");
              Output("ENC_SMOKE_CONTROL_FLOW_COMPLETE search_calls=42 completed_pairs=1 global_abort=0\n");
            }
            break;
          default:break;
          } // switch (input_inp.what_to_do)
        } // end "OPTIMIZE_NC" || "EVALUATE_NC"

        // ================================================
        else if (input_inp.what_to_do == LSL_MINMAX ||
          input_inp.what_to_do == LSL_OPT    ||
          input_inp.what_to_do == ANGLES) 
        {
          input_scp.printSettings();
          input_sim.printSettings();
          // create simulation object from simulation parameter input file and neurocontroller object
          if (sim!=0) delete sim;
          sim = new CSimulation(input_sim, input_scp, NULL, input_inp.report_file);

          // evaluate local steering law vector
          Vector lsl_vector1, lsl_vector2, lsl_vector;

          switch (input_inp.what_to_do)
          {
          case LSL_MINMAX:
            // ------------
            lsl_vector=Vector(
              input_inp.a_max,     input_inp.a_min
              ,input_inp.e_max,     input_inp.e_min
              ,input_inp.i_max,     input_inp.i_min
              ,input_inp.Omega_max, input_inp.Omega_min
              ,input_inp.omega_max, input_inp.omega_min);
            break;

          case LSL_OPT:
            // ---------
            lsl_vector=Vector(
              input_inp.a_opt
              ,input_inp.e_opt
              ,input_inp.i_opt
              ,input_inp.Omega_opt
              ,input_inp.omega_opt);
            break;

          case ANGLES:
            // --------
            lsl_vector=Vector(input_inp.clock, input_inp.cone);
            break;

          default:
            //-----
            FatalError("main(): wrong steering strategy");
            break;
          }

          // evaluate local steering law vector in case of valid simulation object
          if (sim != 0)
          {
#if defined COMMONCPP
            if (m_clock_angle_history == true) sim->activateClockAngleHistory();

            if (m_clock_angle_history_approximation == true) 
              sim->activateClockAngleApproximation();

            sim->setClockAngleApproximationParameters(
              m_clock_angle_approx_seg_size,
              m_clock_angle_approx_poly_degree);

            if (m_cone_angle_history == true) 
              sim->activateConeAngleHistory();

            if (m_cone_angle_history_approximation == true) 
              sim->activateConeAngleApproximation();

            sim->setConeAngleApproximationParameters(
              m_cone_angle_approx_seg_size,
              m_cone_angle_approx_poly_degree);

#endif // COMMONCPP
            sim->evaluate(lsl_vector);
          }

          // write simulation output files
          if (sim != 0) sim->writeOutputFiles(
            input_inp.simulation_output_file_csv
            ,input_inp.simulation_output_file_dat
            ,input_inp.simulation_output_file_gesop
            ,input_inp.simulation_output_file_wrl
            ,input_inp.simulation_output_file_ang);

        } // end "LSL minimize maximize" || "LSL optimize" || "steering angles"
      }
    }
  }
  SafeDelete(sim);
  SafeDelete(ea);
}


//-----------------------------------------------------------------------------
void CInTrance::optimizeApproximation(void) throw(std::exception&)
{
  const unsigned int prec(12); // output precision

  auto_ptr<CNC>         nc(0);   // neurocontroller 1 (clock angle approx.)
  auto_ptr<CNC>         nc2(0);  // neurocontroller 2 (cone angle approx.)
  auto_ptr<CSimulation> sim(0);  // simulation (used for solution evaluating)
  auto_ptr<CEA>         ea(0);   // evolutionary algorithm

#undef NEW_CODE

  // load steering angles history to optimize
  CSCAttitudeHistory attitudes(input_inp.simulation_input_file_ang.c_str());

#if defined NEW_CODE
  // set # of input neurons
  input_ncp.input_neurons  = attitudes.getDimension()+1;
  input_ncp.output_neurons = 2;
#endif

  bool end_current_optimization(false);

  // create approx. neurocontroller objects
	//----------------------------------------------------------------
  if (abort() == false)
  {
    // first, info dump and creation of neurocontrollers
	  try 
    {
      input_ncp.printSettings();
      nc  = auto_ptr<CNC>(new CNC(input_ncp)); 
      nc2 = auto_ptr<CNC>(new CNC(input_ncp));
    }
	  catch (runtime_error& exc)
    {
 	    throw std::runtime_error(string("CInTrance::optimizeApproximation()\n")
        +string(exc.what()));
    }

    // error handling
      if (nc.get() == 0) {
          std::string error_message = "CInTrance::optimizeApproximation(): Memory allocation for neural net controller failed!\n";
          throw std::runtime_error(error_message);
      }
  } // if (abort() == false)

	// create simulation object
  //---------------------------------------------------------------
  if (abort() == false)
  {
    // try to create simulation object
	  try 
    {
      input_scp.printSettings();
      input_sim.printSettings();
	    sim = auto_ptr<CSimulation>(new CSimulation(input_sim
                                                 ,input_scp
                                                 ,nc.get()
                                                 ,input_inp.report_file));
      sim->setNC2(nc2.get());
    }
    // error handling
	  catch (runtime_error& exc)
    {
      throw std::runtime_error(string("CInTrance::optimizeApproximation()\n")
        +string(exc.what()));
    } // catch (runtime_error& exc)

    // again, some error handling is needed
    if (sim.get()==0) {
        std::string error_message = ("CInTrance::optimizeApproximation(void):Memory :\
        allocation for simulation object failed!\n");
        throw std::runtime_error(error_message);
    }

    sim->setSimType(APPROXIMATION);

  } // if (abort() == false)

  if (attitudes.getDimension() > 0)
  {
    sim->evaluate(attitudes);                     // generate reference trajectory
    sim->setReferenceTrajectory();
    sim->extractDataFromReferenceTrajectory();
    std::vector<double> clock_angle = sim->getRefTrajClockAngleArray();
    sim->postprocessAngleData(clock_angle);
    sim->setReferenceTrajectoryClockAngleData(clock_angle);
  } // if (attitudes.getDimension() > 0)
  else throw std::invalid_argument("CInTrance::optimizeApproximation(void):\
 Loaded steering angles history is empty!\n");


  // create evolutionary algorithm object
  //---------------------------------------------------------------
  if (abort() == false)
  {
    // try creation of evolutionary algorithm object
	  try 
    {
      // chromosome size
      input_eap.setChromosomeLength(nc->getNumberOfParameters()
        + nc2->getNumberOfParameters());
      input_eap.printSettings();
      ea = auto_ptr<CEA>(
        new CEA(input_eap
        ,nc.get()
        ,sim.get()
        ,m_rep
        ,PrintReproduction2
        ,0
        ,sim.get()
        )
      );
      if (ea.get() != 0) ea->pushNc(nc2.get());// set second neurocontroller
    }
    // error handling
	  catch (runtime_error& exc)
    {
      throw std::runtime_error(string("CInTrance::optimizeApproximation()\n")
        +string(exc.what()));
    }

    // again, some error handling is needed
    if (ea.get()==0) {
        std::string error_message = ("CInTrance::optimizeApproximation(void):Memory :\
 allocation for evolutionary algorithm object failed!\n");
    throw std::runtime_error(error_message);
  }



  } // if (abort() == false)

  ea->setStopBoolPointer(&m_break);          // inform ev. alg. about break flag pointer

  // run evolutionary algorithm and write possible results to files
  //---------------------------------------------------------------
  if (abort() == false)
  {
    // try optimization run
	  try 
    {
      switch (input_inp.what_to_do)
      {
      case OPTIMIZE_APPROX:
        // do optimization
        ea->Optimize(input_inp.evolutionary_algorithm_warmstart_file.c_str());
        // save found best solution
        ea->SaveBestSolution(input_inp.evolutionary_algorithm_best_solution_file.c_str());
      
        Output("\n best solution:\n\n");
        break;
      case EVALUATE_APPROX:
        Output("\n evaluated solution:\n\n");
        break;
      }

      CIndividual solution(ea.get(),
        input_inp.evolutionary_algorithm_best_solution_file.c_str());
      sim->evaluate(solution, false);
      PrintReproduction2(&solution, 1.0, 1, 1, sim.get());
      Output("\n");
      sim->evaluate(solution, true);


      // write simulation output files
      if (sim.get() != 0) 
      {
        cout.precision(prec);
        sim->writeOutputFiles(
           input_inp.simulation_output_file_csv
          ,input_inp.simulation_output_file_dat
          ,input_inp.simulation_output_file_gesop
          ,input_inp.simulation_output_file_wrl
          ,input_inp.simulation_output_file_ang);

        sim->saveClockApproxParams(input_inp.clock_approx_data_file);
        sim->saveConeApproxParams(input_inp.cone_approx_data_file);
//        sim->PrintRefRotTraj(input_inp.rot_state_traj_file);
      } // if (sim.get() != 0) 
    } // try block
    // error handling
	  catch (runtime_error& exc)
    {
      throw std::runtime_error(string("CInTrance::optimizeApproximation()\n")
        +string(exc.what()));
    } // catch (runtime_error& exc)
  } // if (abort() == false)
} // CInTrance::optimizeApproximation

// end of implementation
