/*!
\file CInTrance.h
\version 1.7
\brief CInTrance class interface.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 23.10.2006 (AO)

\b Purpose:
   This file provides the interface for the class CInputInTrance.
   CInputInTrance implements the functionality of the InTrance program.

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
#ifndef INTRANCE_H_INCLUDED
#define INTRANCE_H_INCLUDED

// InTrance
#include <configInTrance.h>

// STL
#include <string>
#include <list>

// JAVA (for GUI, if desired)
#ifdef JAVA_GUI
 #include <jni.h>    // JAVA Native Interface
#endif

#if defined COMMONCPP
#define _CRT_SECURE_NO_DEPRECATE
 #include <cc++/thread.h>
#endif

// SpaceLib
#include <SpaceLib.h>

// ENCLib
#include <CInputNC.h>
#include <CInputEA.h>
#include <CEATypes.h>

// SCLib
#include <CInputSC.h>

// InTrance
#include <CInputInTrance.h>
//#include <CInputSimulation.h>
#include <CSimulation.h>


#include <time.h>


// forward declarations
class CIndividual;                         // inidividual of ev. alg. population
class SpaceLib::CSolarSystem;                              // solar system class

typedef std::list<std::string> StringList;

// constants
const double HCSF_INC1 = 0.01;  ///< hypercube shrinking factor increment 1 (small)
const double HCSF_INC2 = 0.1;   ///< hypercube shrinking factor increment 2 (bigger)
const double MR_INC1   = 0.01;  ///< mutation rate increment 1 (small)
const double MR_INC2   = 0.1;   ///< mutation rate increment 2 (bigger)

#if defined COMMONCPP
/*!
\brief Class that implements data output in a JAVA GUI. It is used like
cout from C++ standard library.

*/
class GUIOutput: virtual public Global::SLOutput
{
  public:
  GUIOutput(void):SLOutput() {};
  virtual ~GUIOutput(void) {};
  virtual SLOutput& operator<<(const char* message);
};
#endif // COMMONCPP


/// InTrance functionality program
#if defined COMMONCPP
class CInTrance : public ost::Thread
#else
class CInTrance
#endif //COMMONCPP
{
  bool                    m_initialized;    ///< InTrance object is initialized?
  bool                    m_is_running;///< true in case of running optimization
  static bool             m_break;     ///< flag for indicating stopping of opt.
  SpaceLib::CSolarSystem* m_solsys;                     ///< solar system object
  StringList              m_input_files;       ///< input file storage structure
  std::string             m_rep;                        ///< name of report file
  static CSimulation*     sim;                            ///< simulation object
  static CEA*             ea;                 ///< evolutionary algorithm object
  EAStatusRec*            m_ea_status;                     ///< ea status object

  time_t t_start;                                   ///< optimization start time
  time_t t_stop;                                      ///< optimization end time

#if defined (COMMONCPP)
  GUIOutput* m_gui_output;

  /// flag to control clock angle history computation
  bool       m_clock_angle_history;
  /// calculate approximation of steering angle alpha
  bool       m_clock_angle_history_approximation;

  /// degree of polynomials that should approximate steering angle alpha history
  unsigned int m_clock_angle_approx_poly_degree;

  /// number of points per segment to use for approximation
  unsigned int m_clock_angle_approx_seg_size;

  /// flag to control clock angle history computation
  bool       m_cone_angle_history;
  /// calculate approximation of steering angle alpha
  bool         m_cone_angle_history_approximation;

  /// degree of polynomials that should approximate steering angle betahistory
  unsigned int m_cone_angle_approx_poly_degree;

  /// number of points per segment to use for approximation
  unsigned int m_cone_angle_approx_seg_size;

  bool m_thrust_history; ///< thrust history data switch

#endif //COMMONCPP

  static CInputInTrance input_inp;   
               ///< object that contains the contents of the InTrance input file
  CInputSimulation input_sim; 
         ///< object that contains the contents of the simulation parameter file
  CInputSC input_scp;         
         ///< object that contains the contents of the spacecraft parameter file
  CInputNC input_ncp;         
                 ///< object that contains the contents of the NC parameter file
  CInputNC input_sncp;        
        ///< object that contains the contants of steering nc configuration file
  static CInputEA input_eap;
                 ///< object that contains the contents of the EA parameter file

  void reset(void);                                 ///< resets CInTrance object

  void resetTimer(void);         ///< set internal time variable to current time
  void stopTimer(void);              ///< set stop time with current system time

  /// get time span between two subsequent "restTimer", "stopTimer" calls
  double getOptimizationTimespan(void) const;

  /// return currently used simulation object
  static CSimulation* getCurrentSimObj(void) throw();

  /// return currently used evolutionary algorithm  object
  static CEA* getCurrentEaObj(void) throw();

  /// initializes optimization
  void initialize(const std::string& inp_filename) throw (std::logic_error&);

  void initSolarSystem(Global::ParamList& p);
  void initAsteroidDataFile(Global::ParamList& p);
  void initCometDataFile(Global::ParamList& p);

  /// sets fitnessfunction type according to ev. alg. input object
  void initializeFitnessFunctionType(void) throw(std::invalid_argument&);

private:
  static void printHelpScreen(void);
  static void handleEaEvents(const EAEvent);
  static void handleKeyboardInputs(const int first, const int second);
private:

  // method is called by initializeFitnessFunctionType
  void initializeFitnessFunctionType2(const int fft2)
    throw (std::invalid_argument&);
  void initializeFitnessFunctionType3(const int fft3)
    throw (std::invalid_argument&);
public:

  /// reads parameters from input file
  void readConfigParams(
    ParamList  & parameters,
    const std::string& filename,
    const std::string& proc_name
    )
    throw (std::invalid_argument&, std::logic_error&);

  /// read simulation input data
  void readSimConfigParams(void)
    throw (std::invalid_argument&, std::logic_error&);

  /// read spacecraft input data
  void readScConfigParams(void)
    throw (std::invalid_argument&, std::logic_error&);

  /// read evolutionary algorithm input data
  void readEaConfigParams(void)
    throw (std::invalid_argument&, std::logic_error&);

  /// read neurocontroller input data
  void CInTrance::readNcConfigParams(const bool extra_params = true)
    throw (std::invalid_argument&, std::logic_error&);

  /*! \brief returns true if flag "m_break" is true, false if otherwise 
   *or not set
   */
  inline bool abort(void)
  {
    return m_break;
  }

  static inline void setAbort(bool shall = true) {m_break=shall;}

  /*!
   * \brief Output of given message through global output object and report file.
   * The report file output will only be done, when a valid file name has been
   * stored (Through input file parsing).
   * \param mess text message to be put out
   */
  inline void statusMessage(const std::ostringstream& mess)
  {
    Output(mess.str()); // output through global output object

    // and report file
    if (m_rep.size() != 0)
    {
      ofstream report_file(m_rep.c_str(), ios::app);
      if (report_file.is_open() == true) report_file << mess.str();
    }
  };


  /// evaluate given set of steering angles
  void evaluateControl(void);

  /// input of trajectory optimization/evaluation parameters
  void initializeControlEvaluation(void) throw (std::invalid_argument&);

  /// input of local steering laws optimization/evaluation parameters
  void initializeLslOptimization(void) throw (std::invalid_argument&);

  /// input of trajectory optimization/evaluation parameters
  void initializeNcOptimization(void) throw (std::invalid_argument&);
  // optimize NC (trajectory)
  void optimizeNc(void) throw (std::exception&);

  /// input of approximation configuration
  void initializeApproximationOptimization(void) throw (std::invalid_argument&);
  /// optimize steering angles spline approximation
  void optimizeApproximation(void) throw(std::exception&);

  void insertConfigFilename(const std::string& filename);

public:

  /// output method
  inline void printMessage(const std::ostringstream& m) {statusMessage(m);};

  /// InTrance input file constructor
  CInTrance(const char* inp_filename);

  /// destructor
  CInTrance::~CInTrance(void) {
    std::cerr << "EXIT_TRACE CInTrance: before solar system delete "
              << static_cast<const void*>(m_solsys) << std::endl;
    SafeDelete(m_solsys);
    std::cerr << "EXIT_TRACE CInTrance: solar system deleted, member teardown next" << std::endl;
  }

  /// Runs InTrance.
  void run(void);

  /// give status object to underlying EA to continously get ea status data
  void setEAStatusObject(EAStatusRec* status);

  /// is initialized yet?
  bool isInitialized(void) {return m_initialized;}

  /// stops current optimization process
  void stop(void);

  /// returns true if InTrance is running an optimization process
  bool isRunning(void) {return m_is_running;}

#if defined (COMMONCPP)
  /// Tbd.
  void activateClockAngleHistory(void);
  /// Tbd.
  void deactivateClockAngleHistory(void);
  /// Tbd.
  void setClockAngleApproximationParameters(
    const unsigned int segSize
   ,const unsigned int polyDegree);

  /// Tbd.
  void activateClockAngleHistoryApproximation(void);
  /// Tbd.
  void deactivateClockAngleHistoryApproximation(void);

  // cone angle control methods
  void activateConeAngleHistory(void);
  void deactivateConeAngleHistory(void);

  void setConeAngleApproximationParameters(
   const unsigned int segSize
  ,const unsigned int polyDegree);

  void activateConeAngleHistoryApproximation(void);
  void deactivateConeAngleHistoryApproximation(void);

  /// activate thrust data preparation
  void activateThrustHistory(void);
  /// decativate thrust history
  void deactivateThrustHistory(void);
#endif //COMMONCPP

}; // end class CInTrance


#if defined COMMONCPP
inline void CInTrance::activateClockAngleHistory(void)
{
  m_clock_angle_history = true;
  if (sim != 0) sim->activateClockAngleHistory();
} // void CInTrance::activateClockAngleHistory(void)

inline void CInTrance::deactivateClockAngleHistory(void)
{
  m_clock_angle_history = false;
  if (sim != 0) sim->deactivateClockAngleHistory();
} // void CInTrance::deactivateClockAngleHistory(void)


inline void CInTrance::setClockAngleApproximationParameters(
  const unsigned int segSize
  ,const unsigned int polyDegree)
{
  m_clock_angle_approx_poly_degree = polyDegree;
  m_clock_angle_approx_seg_size = segSize;
  if (sim != 0) sim->setClockAngleApproximationParameters(segSize,polyDegree );
} // void CInTrance::setClockAngleApproximationParameters(...)


inline void CInTrance::activateClockAngleHistoryApproximation(void)
{
  m_clock_angle_history_approximation = true;
  if (sim != 0) sim->activateClockAngleApproximation();
} // void CIntrance::activateClockAngleHistoryApproximation(void)

inline void CInTrance::deactivateClockAngleHistoryApproximation(void)
{
  m_clock_angle_history_approximation = false;
  if (sim != 0) sim->deactivateClockAngleApproximation();
} // void CIntrance::deactivateClockAngleHistoryApproximation(void)

inline void CInTrance::activateConeAngleHistory(void)
{
  m_cone_angle_history = true;
  if (sim != 0) sim->activateConeAngleHistory();
} // void CInTrance::activateConeAngleHistory(void)

inline void CInTrance::deactivateConeAngleHistory(void)
{
  m_cone_angle_history = false;
  if (sim != 0) sim->deactivateConeAngleHistory();
} // void CInTrance::deactivateConeAngleHistory(void)

inline void CInTrance::setConeAngleApproximationParameters(
  const unsigned int segSize
  ,const unsigned int polyDegree)
{
  m_cone_angle_approx_poly_degree = polyDegree;
  m_cone_angle_approx_seg_size = segSize;
  if (sim != 0) sim->setConeAngleApproximationParameters(segSize,polyDegree);
} // void CInTrance::setConeAngleApproximationParameters(

inline void CInTrance::activateConeAngleHistoryApproximation(void)
{
  m_cone_angle_history_approximation = true;
  if (sim != 0) sim->activateConeAngleApproximation();
} // void CInTrance::activateConeAngleHistoryApproximation(void)

inline void CInTrance::deactivateConeAngleHistoryApproximation(void)
{
  m_cone_angle_history_approximation = false;
  if (sim != 0) sim->deactivateConeAngleApproximation();
} // void CInTrance::deactivateConeAngleHistoryApproximation(void)

//------------------------------------------------------------------------------
inline void CInTrance::activateThrustHistory(void)
{
  m_thrust_history = true;
  if (sim != 0) sim->activateThrustHistory();
} // void CInTrance::activateThrustHistory(void)

//------------------------------------------------------------------------------
/// decativate thrust history
inline void CInTrance::deactivateThrustHistory(void)
{
  m_thrust_history = false;
  if (sim != 0) sim->deactivateThrustHistory();
} // void CInTrance::deactivateThrustHistory(void)

#endif // COMMONCPP

//------------------------------------------------------------------------------
inline void CInTrance::resetTimer(void)
{
  t_stop = t_start = time(0);
} // void resetTimer(void)

//------------------------------------------------------------------------------
inline void CInTrance::stopTimer(void)
{
  t_stop = time(0);
} // void stopTimer(void)

//------------------------------------------------------------------------------
inline double CInTrance::getOptimizationTimespan(void) const
{
  return difftime(t_stop, t_start);
} // time_t getOptimizationTimespan(void) const

//------------------------------------------------------------------------------
inline CSimulation* CInTrance::getCurrentSimObj(void) throw()
{
  return sim;
} // CSimulation* getCurrentSimObj(void) throw()

//------------------------------------------------------------------------------
inline CEA* CInTrance::getCurrentEaObj(void) throw()
{
  return ea;
}// CEA* CInTrancegetCurrentEaObj(void) throw()

#endif // INTRANCE_H_INCLUDED
