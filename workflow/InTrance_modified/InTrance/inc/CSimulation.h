/*!
\file CSimulation.h
\version 1.7
\brief CSimulation class interface.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 12.11.2006 (AO)

\b Purpose:
This file provides the interface for the class CSimulation.
CSimulation is the class that implements the simulation functionality.

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

#ifndef CSIMULATION_H_INCLUDED
#define CSIMULATION_H_INCLUDED

#define using_constraints
#include <algorithm>
// InTrance
#include <InTranceConst.h>
#include <CInputSimulation.h>
#include <CInTranceFitness.h>
#include <CConstraint.h>

#if defined COMMONCPP
#include <TQueue.h>
#undef min // this is important, because there're templates defined in
#undef max // namespace Global that do this job
#else
typedef std::vector<MathLib::Vector> CTrajectory;
typedef CTrajectory*                 Trajectory;
typedef std::vector<double>          CDataArray;
typedef CDataArray*                  DataArray;
#endif //COMMONCPP

// SpaceLib
#include <MLDE.h>
#include <SpaceLib.h>

// ENCLib
#include <CIndividual.h>
#include <CEAEvalFunc.h>
#include <CEAFitness.h>
#include <CEA.h>

// SCLib
#include <CSC.h>

// forward declarations
class CNC;
class CSCAttitudeHistory;
class CInputSC;
class CInputSimulation;

// declarations

/// exchange struct for transfer of simulation data to Java GUI
struct SimStatusRec
{
  /// list of character strings holding the exchange information of current 
  /// simulation
  std::vector<std::string> data;
}; // struct SimStatusRec


/// simulation type
enum SimulationType
{
   STANDARD      ///< trajectory optimization
  ,APPROXIMATION ///< steering angles optimization
}; // SimulationType

/// trajectory data structure (vector of spacecrafts)
typedef std::vector<CSC> SCTraj;

/// max. # of integrator object calls per evaluation
const int N_MAX_ITGR_CALLS = 50000;
const double INIT_BEST_PROX_SUM = 1.0e20;
const double DFLT_ACCURACY_FACTOR = 1.0;

/*!
* \brief Class that implements simulation functionalities.
* \todo BD: remove outcommented two lines in calcBestValues member function
* \todo AO: evolve template version of evaluateSteeringAngleApprox; then delete outcommented implementation
*/
class CSimulation : virtual public CEAEvalFunc
{
protected:
  ///< all set trajectory integration constraint objects
  std::vector<constraints::CConstraint*> m_constraints;

  SimStatusRec* m_status;  ///< simulation data transfer object (prim. Java GUI)

  CInputSimulation m_data;                 ///< simulation configuration object

  /// dynamic accuracy member fields
  double m_proposed_accuracy_factor;
  double m_accuracy_factor;

  double m_curr_DSSC_max_angle;  ///< current setting of max. control step angle

  // deflection optimization specific variables
  double m_impact_time;                       ///< time of impact on target body
  double m_impact_efficiency;           ///< efficiency of impact on target body
  double m_deflection;                                   ///< deflection measure

  double m_initial_velocity;              ///< absolute excess velocity at start
  double m_initial_azimuth;     ///< launch excess velocity vector azimuth angle
  double m_initial_elevation; ///< launch excess velocity vector elevation angle
  double m_initial_mp;                             ///< propellant mass at start
  double m_initial_proximity;   ///< prox. sum value at start (orbital elements)
  double m_best_initial_proximity;     ///< best value ever of initial prox. sum
  double m_h;                                         ///< current step size [d]
  double m_h_dyn;                             ///< current dynamic step size [d]
  double m_h_min;                             ///< current minimal step size [d]
  double m_t;                                      ///< variable of current time
  double m_time_start;                                ///< trajectory start time
  double m_time_max_dynamic;              ///< max. dynamic integration interval
  double m_time_best;                    ///< time of best objective achievement
  double m_ga_credits;         ///< gravitiy assist credits for fitness function
  double m_best_dr;                  ///< relative distance measure of best step
  double m_best_dv;                  ///< relative velocity measure of best step
  double m_best_da;   ///< rel. deviation measure of semi-major axis (best step)
  double m_best_de;      ///< rel. deviation measure of eccentricity (best step)
  double m_best_di;       ///< rel. deviation measure of inclination (best step)
  double m_best_dO;    ///< rel. deviation measure of ascending node (best step)
  double m_best_do;   ///< rel. deviation measure of periapsis angle (best step)
  double m_best_dM;      ///< rel. deviation measure of mean anomaly (best step)
  double m_best_dh;  ///< rel. deviation measure of orbital momentum (best step)
  double m_best_impv;                          ///< impact velocity of best step
  double m_best_rimph;     ///< rel. orbital momentum measure of retrogr. impact
  double m_best_rimpe;           ///< rel. eccentrity measure of retrogr. impact

  double m_best_act; ///< best steps's ACT fitness (special fitness value)
  ///< 2005, 1st ESA ACT Optimization Competition

  double m_best_mp;                            ///< propellant mass at best step
  double m_best_deflection;     ///< greatest deflection (distance) of best step

  // may be deleted
  double m_best_spi_a;
  double m_best_spi_i;
  double m_best_spi_e;

  double m_launch_prox_sum;
  double m_temp_prox_sum;
  double m_best_prox_sum;     /// integration step of best objective achievement

  unsigned int m_gtoc2_fp;                                      /// flight phase
  std::vector<unsigned int> m_gtoc2_group;
  std::vector<BodyID> m_gtoc2_temp_best_bodies;
  std::vector<BodyID> m_gtoc2_best_bodies;
  CBody* m_gtoc2_temp_asteroid;

  double m_sfa;                        ///< sensitivity factor of semimajor axis
  double m_sfe;                          ///< sensitivity factor of eccentricity
  double m_sfi;                           ///< sensitivity factor of inclination
  double m_sfh;         ///< sensitivity factor of orbit angular momentum vector

  double m_max_dr_km;         ///< maximum final distance to target body [in km]
  double m_max_dr_AU;                                           ///< ... in [AU]
  double m_max_dv_kms;///< max. final relative velocity to target body [in km/s]
  double m_max_dv_AUd;                                        ///< ... in [AU/d]

  std::vector<double> m_ACT_pericenter_times;

  /// initial deviation of Keplerian elements in case of Keplerian elements
  /// target
  std::valarray<double> m_initial_delements;
  /// current deviations of Keplerian elements of date to desired target values
  std::valarray<double> m_delements;
  /// estimated remaining time to change to target element value
  std::valarray<double> m_element_change_rest_times;
  /// best deviation values of Keplerian elements to desired target values
  std::valarray<double> m_best_delements;
  /// best estimated remaining time to change to target element value
  std::valarray<double> m_best_element_change_rest_times;

  SCTraj m_last_saved_trajectory;                 ///< last saved s/c trajectory

  SimulationType m_sim_type;                             ///< type of simulation

  std::string m_rep;                                    ///< name of report file
  std::string m_csv_output_filename;                ///< name of csv output file

  bool* m_should_abort;                               ///< stop variable pointer
  bool* m_ga_body;      ///< pointer to array of potential gravity assist bodies
  bool m_ga_performed;                  ///< gravity assist maneuver performed ?

  SpaceLib::CSolarSystem* m_solsys;                             ///<solar system

  CSC m_launch_sc;                       ///< CSC object at trajectory beginning
  CSC m_spacecraft;                              ///< present spacecraft (state)
  CSC m_best_sc;                   ///< CSC object of best objective achievement

  SCTraj m_spacecraft_trajectory;                     ///< spacecraft trajectory
  SCTraj m_ref_sc_traj;                     ///< reference spacecraft trajectory

#if defined COMMONCPP
  ost::Mutex m_sc_traj_lock;  ///< mutex for data transfer purposes (trajectory)
  ost::Mutex m_ref_sc_traj_lock;    ///< mutual exclusion lock for data transfer
  ost::Mutex m_lock;                          ///< general mutual exclusion lock
  ost::Mutex m_last_saved_trajectory_lock;                  ///< trajectory lock
  ost::Mutex m_last_saved_cone_angle_history_lock;          ///< cone angle lock
  ost::Mutex m_last_saved_clock_angle_history_lock;;       ///< clock angle lock
  ost::Mutex m_last_saved_thrust_history_lock;  ///< thrust setting history lock
#endif // COMMONCPP

  /// data array of thrust settinh history
  CDataArray m_last_saved_thrust_history;
  /// data structure of las saved clock angle setting history
  CDataArray m_last_saved_clock_angle_history;
  /// data structure of las saved cone angle setting history
  CDataArray m_last_saved_cone_angle_history;

  std::ofstream m_csv_stream;             ///< file stream object for csv-output

  int m_itgr_dim;                              ///< integration vector dimension
  int m_steps;                                           ///< current step count
  int m_steps_max_used;        ///< upper boundary for # of steps <= m_steps_max
  int m_steps_best;                                     ///< number of best step

  MathLib::Vector m_ys;                  ///< current state vector (integration)
  MathLib::Vector m_x;                                      ///< position vector
  MathLib::Vector m_v;                                      ///< velocity vector
  MathLib::Vector m_best_s;              ///< state vector of best result so far

  CProximity m_curr_proximity;                            ///< current proximity
  CProximity m_temp_proximity;                       ///< temp. proximity object
  CProximity m_init_proximity;                  ///< initial proximity at launch
  CProximity m_best_proximity;                       ///< best steps's proximity

  RK4 *              m_SpacecraftOrbit_RK4;           ///< Runge-Kutta 4th order
  /// Runge-Kutta 4th order, 5th order predictor
  RKF54 *            m_SpacecraftOrbit_RKF54;
  /// Dormand & Prince 8th order
  DOPRI8 *           m_SpacecraftOrbit_DOPRI8;
  DE *               m_SpacecraftOrbit_DESG;

  CIndividual        m_last_saved_individual;

  CNC* m_nc;                                       ///< steering neurocontroller
  CNC* m_nc2;                                 ///< approximation neurocontroller

  
protected:
  /// initializes the different integration constraints
  virtual void initConstraints(void);

private:
  /// initializes gravitiy assist bodies
  void initGravityAssistBodies(void);

  /// check for correct optimization goal dependend on spacecraft type
  void checkOptimizationGoal(void) throw (std::invalid_argument);

  /// sets internal initial conditions depending on given configuration
  void setInitialState(void) throw (std::invalid_argument);

  /// set initial condition of a given celestial body
  void setBodyInitialState(void);

  /// set intital state of given orbital parameters
  void setOrbitInitialState(void) throw (std::logic_error);

  /// set initial state of given position/velocity vector
  void setVectorInitialState(void);

  /// initialize spacecraft object
  void initSpacecraft(const CInputSC& scp_input) throw (std::logic_error);

  /// sets internal target conditions depending on input configuration
  virtual void setTargetState(void) throw (std::invalid_argument);

  // private methods for setting up different types of target states
  void setTargetStateBody(void);
  void setTargetStateBodyOrbit(void);
  void setTargetStateOrbit(void);
  void setTargetStateElements(void);
  void setTargetStateReserved01(void);

  /// array of deviations in Keplerian elements at trajectory start
  std::valarray<double>& getTOEInitialElementDeviation(void);

  /// set new integration step size
  void setIntegrationStepSize(const double step_size);
  /// return current integration step size
  double getIntegrationStepSize(void);

  /// sets Runge-Kutta 4th order integrator object
  void setIntegrator(MathLib::RK4* itgr);
  /// return pointer currently used Runge-Kutta 4th order integrator object
  MathLib::RK4* getIntegratorRK4(void) const;
  /// sets Runge-Kutta 4th order predictor integrator object
  void setIntegrator(MathLib::RKF54* itgr);
  /// return currently used RKF54 integrator
  MathLib::RKF54* CSimulation::getIntegratorRKF54(void) const;
  /// sets Runge-Kutta 8th order integrator object
  void setIntegrator(MathLib::DOPRI8* itgr);
  /// return currently used 8th order integrator object
  MathLib::DOPRI8* CSimulation::getIntegratorDOPRI8(void) const;
  /// sets Runge-Kutta 4th order multistep integrator object
  void setIntegrator(MathLib::DE* itgr);
  /// return currently used multistep integrator object
  MathLib::DE* CSimulation::getIntegratorDESG(void) const;

  /// help screen of commanline version
  void printHelpScreen(void);
  /// transfer rendezvous simulation data to Java GUI transfer data structure
  void fillStatusRecRendezvous(SimStatusRec* s);
  /// transfer fly-by simulation data to Java GUI transfer data structure
  void fillStatusRecFlyby(SimStatusRec* s);
  /// transfer impact simulation data to Java GUI transfer data structure
  void fillStatusRecRetrogradeImpact(SimStatusRec* s);
  /// transfer ACT-problem simulation data to Java GUI transfer data structure
  void fillStatusRecActProblem(SimStatusRec* s);

  /*! \brief Calculates and prints probably good launch windows for optimizing
  * transfer to an asteroid with repect to the 1st ESA ACT competition held in 
  *november 2005*/
  void PrintACTCompetitionLaunchWindows(void);

  /// Initialize equations of motion
  void initEOM(CSC * sc);

  /// Initialize equations of motion for user provided integrator objects
  void CSimulation::initEOM(
    CSC* sc
    ,RK4*&    integrator_RK4
    ,RKF54*&  integrator_RKF54
    ,DOPRI8*& integrator_DOPRI8
    ,DE*&     integrator_DESG);

  /// free integrator object(s)
  void freeEOM(void);

  /// initialize integrator object(s)
  void CSimulation::initDESSolver( 
    CSC& sc
    ,DESSolverType solver_type
    ,MathLib::Vector& ys
    ,const double t
    ,MathLib::RK4*    integrator_RK4
    ,MathLib::RKF54*  integrator_RKF54
    ,MathLib::DOPRI8* integrator_DOPRI8
    ,MathLib::DE*     integrator_DESG  
    );


  /// initialize numerical equation solver
  void initDESSolver(Vector & ys, const double t);


  /// calculate one integration step
  RK_STATE CSimulation::integrate(
    CSC&          sc
    ,Vector &      ys
    ,double &      t
    ,const double  h
    ,DESSolverType des_solver
    ,RK4*          integr_RK4
    ,RKF54*        integr_RKF54
    ,DOPRI8*       integr_DOPRI8
    ,DE*           integr_DESG
    );


  /// calculate one integration step
  RK_STATE integrate(Vector & ys, double & t, const double h);

  /// pre-integration processing
  void integrationPreProcessing(void);

  /// pre-integration processing for orbital elements target type
  void integrationPreProcessingTOES(void);

  /// initialize sensitivity factors
  void initSensitivityFactors(void);

  void calcBestValuesGTOC2(void);
  void calcBestValuesTest1(void);
  void calcBestValuesTest2(void);
  void calcBestValuesTest3(void);
  void calcBestValuesTest4(void);
  void calcBestValuesTest5(void);
  void calcBestValuesTest6(void);
  void calcBestValuesTest7(void);
  void calcBestValuesTest8(void);
  void calcBestValuesTest9(void);

  /// calc "best" values for rendezvous problem on new best step
  void calcBestValuesRendezvous(void);

  /// calc "best" values for fly-by problem on new best step
  void calcBestValuesFlyby(void);

  /// calc "best" values for ACT problem on new best step
  void calcBestValuesACT(
    BodyID        target_body
    ,CSC&          sc
    ,CProximity&   proximity
    ,CSolarSystem* sol_sys
    ,double&       best_dr
    ,double        max_dr_AU
    ,double&       best_act
    ,double&       best_mp
    );

  /// calc "best" values for SPI problm on new best step
  void calcBestValuesSPI(void);

  /// calc "best" values for retrograde impace problem on new best step
  void calcBestValuesRetrogradeImpact(
    BodyID          target_body
    ,CSC&            sc
    ,CProximity&     proximity
    ,CSolarSystem*   sol_sys
    ,double&         best_dr
    ,double          max_dr_AU
    ,double&         best_dv
    ,double          max_dv_AUd
    ,double&         best_impv
    ,double&         best_rimph
    ,double&         best_rimpe
    ,double&         best_mp
    );

  /// calc "best" values on new best step
  void calcBestValues(
    BodyID          target_body
    ,CSC&            sc
    ,CProximity&     proximity
    ,CSolarSystem*   sol_sys
    ,double&         best_dr
    ,double          max_dr_AU
    ,double&         best_dv
    ,double          max_dv_AUd
    ,double&         best_da
    ,double&         best_de
    ,double&         best_di
    ,double&         best_dO
    ,double&         best_do
    ,double&         best_dM
    ,double&         best_dh
    ,double&         best_impv
    ,double&         best_rimph
    ,double&         best_rimpe
    ,double&         best_act
    ,double&         best_spi_a
    ,double&         best_spi_i
    ,double&         best_spi_e
    ,double&         best_mp
    ,double&         ga_credits
    ,bool            ga_perfomed // gravity assist performed ?
    ,double          sf_a        // sensitivity factor of semi major axis
    ,double          sf_e        // sensitivity factor of eccentricity
    ,double          sf_h        // sensitivity factor of angular momentum
    );


  /// calc "best" values on new best step
  void calcBestValues(void);

  void setFinalTargetStateOrbitalElements(void);

  /*! \brief Sets target state after trajectory integration during evaluation.
  * This is problem-dependend.*/
  void setFinalTargetState(void);


  CEAFitness* setFitnessValuesTest0(CEAFitness* fitness, double time);
  CEAFitness* setFitnessValuesTest1(CEAFitness* fitness, double time);
  CEAFitness* setFitnessValuesTest2(CEAFitness* fitness, double time);
  CEAFitness* setFitnessValuesTest3(CEAFitness* fitness, double time);
  CEAFitness* setFitnessValuesTest4(CEAFitness* fitness, double time);
  CEAFitness* setFitnessValuesTest5(CEAFitness* fitness, double time);
  CEAFitness* setFitnessValuesTest6(CEAFitness* fitness, double time);
  CEAFitness* setFitnessValuesTest7(CEAFitness* fitness, double time);
  CEAFitness* setFitnessValuesTest8(CEAFitness* fitness, double time);
  CEAFitness* setFitnessValuesTest9(CEAFitness* fitness, double time);

  /// set fitness object for ACT problem
  CEAFitness* setFitnessValuesACT(CInTranceFitness*);

  /// set fitness object of 2nd Global Trajectory Optimization Competition
  CEAFitness* setFitnessValuesGTOC2(CInTranceFitness*);
  
  /// set fitness object for orbital elements target problem
  CEAFitness* setFitnessValuesOrbElements(CInTranceFitness*);

  CEAFitness* setFitnessValuesFlyBy(CEAFitness*);

  /*!
  * \brief set fitness object according to best evaluation criteria values
  * so far
  */
  CEAFitness* setFitnessValues(CEAFitness* fitness_input, double time_slice);


  /*! \brief set start state according to given parameters
  \param time start time
  \param velocity start velocity value
  \param azimuth start angle in orbit plane with respect to 
  orbit tangential unit vector
  \param elevation start angle from orbit plane to orbit normal vector*/
  void setLaunchState(const double time,    const double velocity,
    const double azimuth, const double elevation);

  /// set new proximity object of first step
  void setInitProximity(CProximity& prox);

  /// set new proximity object of best step
  void setBestProximity(CProximity& prox);

  /// update proximity object in case of body target
  void updateProximityTargetBody(void);

  void updateProximityGtoc2(void);

  /// update proximity object in case of body orbit target
  void updateProximityTargetBodyOrbit(void);

  /// update proximity object in case of free orbit target
  void updateProximityTargetFreeOrbit(void);

  /// update proximity object in case of orbital elements target
  void updateProximityTargetOrbitalElements(void);

  /// do steering (estimate new s/c orientation)
  void performSteering(void);

  /// add current spacecraft object to trajectory
  void addToSpacecraftTrajectory(const CSC & sc);


  /*!
  * \brief Output of given message through global output object and report file.
  * The report file output will only be done, when a valid file name has been
  * stored (Through input file parsing).
  * \param mess text message to be put out
  */
  void statusMessage(const std::ostringstream& mess);

  /// Tbd.
  void calcProposedAccuracy(void);

  /// Tbd.
  void performGravityAssist(Vector& ys, Vector& x, Vector& v, double& t);

  /*!
  * \brief Special function to obtain the deviation of an asteroid's trajectory
  * from its original path due to impact of a spacecraft (s/c)
  *
  * \param t0 time of s/c impact on asteroid
  * \param t2 end time for evaluation of new trajectory and difference estimating 
  *to old trajectory
  * \param nea_r0 Initial asteroid position at t=t0
  * \param nea_v0 Initial asteroid velocity at t=t0
  * \param sc_v0 Spacecraft velocity at t=t0
  * \param impact_eff impact efficency; measure of how much of the s/c's momentum
  *(because of velocity difference vector) is transferred to asteroid
  * \return Norm of difference vector between asteroid's original position at t=t1 and
  * it's position due to s/c impact at t=t1
  */
  double calcDeflection(
    const double t0
    ,const double t2
    ,const Vector& nea_r0
    ,const Vector& nea_v0
    ,const Vector& sc_v0
    ,const double impact_eff
    );

  /// write simulation data to TECPLOT data file
  void PrintDAT  (const char * filename);

  /// write problem-specific data into GESOP file format file
  void PrintGESOP(const char * filename);

  /// write solution to VRML file
  void PrintVRML (const char * filename);

  /// write steering angles to given output file
  void PrintANG  (const char * filename);

  /*!\brief Extracts position and velocity data from given s/c trajectory.
  * \param sc_traj Spacecraft trajectory
  * \param traj Trajectory consisting of position and velocity data.
  * \return Modified parameter traj*/
  void getTrajectoryData(const SCTraj& sc_traj, CTrajectory& traj
    , const unsigned int sz) const;

  /*! \brief Fills a trajectory data structure with position and velocity
  data of given target body at times from given s/c trajectory.

  \param sc_traj 
  \param target
  \param traj
  \return modified parameter traj */
  void getTrajectoryTargetData(const SCTraj& sc_traj,const BodyID target
    , CTrajectory& traj);

  /// Tbd.
  int getBasicClockAngleData(const SCTraj& traj, double*& angle_values
    ,const int n_max_points);

  /// Clock angle data set from given s/c trajectory
  void getClockAngleData(const SCTraj& traj, CDataArray& data
    ,const int max_n_items = -1);

  /// Tbd.
  int getBasicConeAngleData (const SCTraj& traj, double*& angle_values
    , const int n_max_points);

  /// Cone angle data set from given s/c trajectory
  void getConeAngleData(const SCTraj& traj, CDataArray& data
    ,const int max_n_items = -1);

  /// Tbd.
  void getBasicThrustData(const SCTraj& t, CDataArray& a, const uint n_max_points = 0);

  /// Initializes sensitivity factors depending on current optimization task.
  void checkSensitivityFactorsInitialization(void);

  /*!\brief True if integration loop of standard evaluation shall continue.
  * \param final final evaluation (true)
  * \param monitor for visualization purposes
  * \return true, if integration loop should be continued*/
  bool continueStandardEvaluationIntegrationLoop(bool final, bool monitor);

  /*! \brief reduces maximum integration time (problem-dependent)
  * If certain boundary conditions are satisfied and the primary
  * optimization criterion(s) are to be maximized/minimized, the
  * integration is performed only to the currently known best (shortest) flight
  * time. This behaviour is problem-dependent.*/
  void checkDynamicMaxIntegrationTimeReduction(void);

  /*! \brief check method in case of target state is orbital elements */
  bool checkForNewBestStepTargetOrbitalElements(void);

  /*! Tbd */
  void checkForNewBestStepRendezvous(void);
  /* Tbd. */
  void checkForNewBestStepRetrogradeImpact(void);
  /*! Tbd */
  void checkForNewBestStepACTProblem(void);
  /*! Tbd */
  void checkForNewBestStepSPIProblem(void);
  /*! Tbd */
  bool checkForNewBestStepGTOC2(void);
  /*! Tbd */
  void checkForNewBestStepTest1(void);
  /*! Tbd */
  void checkForNewBestStepTest2(void);
  /*! Tbd */
  void checkForNewBestStepTest3(void);
  /*! Tbd */
  void checkForNewBestStepTest4(void);
  /*! Tbd */
  void checkForNewBestStepTest5(void);
  /*! Tbd */
  void checkForNewBestStepTest6(void);
  /*! Tbd */
  void checkForNewBestStepTest7(void);
  /*! Tbd */
  void checkForNewBestStepTest8(void);
  /*! Tbd */
  void checkForNewBestStepTest9(void);

  /*!
  * \brief Use dynamic step size control to calculate new integration step size.
  * \return new step size
  */
  double calcDSSCStepSize(void);

  /// Tbd. Estimate size of next integration step.
  double calcIntegrationStepWidth(const double h_min, bool final);

  /// dimension of integration vector
  int getIntegrationVectorDimension(void);
  /// set new dimension of integration vector
  void setIntegrationVectorDimension(const int dim);

  /// Sets problem-specific dimension of integration vector
  void initIntegrationVector(void);

  /// set initial conditions of rendezvous problem
  void setRendezvousInitConditions(CIndividual& i, double time);

  /// set initial conditions of fly-by problem
  void setFlyByInitConditions(CIndividual& i, double time);

  /// ste initial conditions of retrograde impact problem
  void setRetrogradeImpactInitConditions(CIndividual& i, double time);

  /// set initial conditions of 1st Global Trajectory Optimization Competition
  void setActProblemInitConditions(CIndividual& i, double time);

  /// set initial conditions of SPI mission
  void setSpiProblemInitConditions(CIndividual& i, double time);

  /// set initial conditions of 2nd Global Trajectory Optimization Competition
  void setGTOC2InitConditions(CIndividual& i, double time);

  /// Tbd.
  void setTrajTest1InitConditions(CIndividual& i, double time);

  /// Tbd.
  void setTrajTest2InitConditions(CIndividual& i, double time);

  /// Tbd.
  void setTrajTest3InitConditions(CIndividual& i, double time);

  /// Tbd.
  void setTrajTest4InitConditions(CIndividual& i, double time);

  /// Tbd.
  void setTrajTest5InitConditions(CIndividual& i, double time);

  /// Tbd.
  void setTrajTest6InitConditions(CIndividual& i, double time);

  /// Tbd.
  void setTrajTest7InitConditions(CIndividual& i, double time);

  /// Tbd.
  void setTrajTest8InitConditions(CIndividual& i, double time);

  /// Tbd.
  void setTrajTest9InitConditions(CIndividual& i, double time);

  /// Tbd.
  void setStandardEvaluationInitConditions(CIndividual& solution, double time);

  /// Tbd.
  void initStandardEvaluation(CIndividual& solution,
    bool final, bool monitor, double time);

  /*!
  * \brief Standard evaluation function for trajectory optimization
  * \param solution candidate for problem solution to evaluation
  * \param final evaluation (true/false) (with min_output_points)
  * \param monitor special, additional evaluation for trajectory evaluation
  * \param time start time [MJD], if 0 the time encoded in solution is taken
  */
  CEAFitness* evaluateStandardSolution(CIndividual& solution
    ,const bool final = false,const bool monitor = false
    ,const double time = 0.0);

  /// Saves an internal copy of actual s/c trajectory.
  void saveLastTrajectory(void); 

  /// postprocess steps after finishing trajectory integration
  void evaluationPostprocess(void);

public:

  /*!
  * \brief 
  * Initializes the simulation object from input objects of simulation
  * and spacecraft which describe the simulation environment and all
  * needed parameters. Additionally the needed neurocontroller object
  * is set with the given object reference. The parameter rep_file
  * is the name of the desired report file into which all messages
  * will be logged.\n If anything goes wrong during the initialization
  * process a standard exeption (std::invalid_argument) will be raised.
  * \param sim_input simulation input object
  * \param scp_input spacecraft describing object
  * \param nc pointer to neurocontroller
  * \param rep_file_name name of report file
  * \throw std::invalid_argument exeption to signal an error
  * \throw std::logic_error if s/c object has no reference solar system
  */
  CSimulation(
    const CInputSimulation & sim_input,              
    const CInputSC &         scp_input,
    CNC *              nc,
    const std::string& rep_file_name
    ) throw (std::invalid_argument, std::logic_error);

  virtual CSimulation::~CSimulation(void);                       ///< destructor

  virtual bool checkForNewBestStep(void);           ///< check for new best step

  virtual void updateProximity(void);        ///< update proximity object values

  SpaceLib::BodyID getInitCentralBody(void) const;    ///< initial cent. body id

  /// currently set (and used if activated) max. orbital control step angle
  double getCurrMaxDSSCOrbAngle(void) const;

  /// new max. orbital control step angle that is to use if DSSC is activated
  void setCurrMaxDSSCOrbAngle(const double angle);

  /// return current spacecraft object
  const CSC& getCurrentSC(void) const;

  /// non-const version of current spacecraft object
  CSC& getCurrentSC(void);

  /// Evaluates spacecraft control history
  /*! \param att_history spacecraft control history */
  void evaluate(const CSCAttitudeHistory & att_history);

  /// Evaluates local steering strategy vector
  /*! \param lsl_vector local steering strategy vector*/
  void evaluate(const Vector & lsl_vector);

  /*! \brief Write simulation output files
  \param csv_file name of the CSV file
  \param dat_file name of the DAT file
  \param gesop_file name of the GESOP file
  \param wrl_file name of the VRML file
  \param ang_file name of the file containing the steering angles */
  void writeOutputFiles(const std::string& csv_file,
    const std::string& dat_file,
    const std::string& gesop_file,
    const std::string& wrl_file,
    const std::string& ang_file);

  /// Tbd.
  void evaluateInterpolated(void);

  /// Tbd.
  void printMessage(const std::ostringstream& m);

  /// Tbd.
  void writeMessageToReportFile(const std::ostringstream& mess);

  /// Tbd.
  void writeMessageToReportFile(const std::string& mess);

  void fillStatusRec(SimStatusRec* s);    ///< fill data exchange data structure

  void fillStatusRec(void);               ///< fill data exchange data structure

  /// set up status object for GTOC2-problem
  void fillStatusRecGTOC2Problem(SimStatusRec* s);
 
  void fillStatusRecTest1(SimStatusRec* s);                   ///< unimplemented
  void fillStatusRecTest2(SimStatusRec* s);                   ///< unimplemented
  void fillStatusRecTest3(SimStatusRec* s);                   ///< unimplemented
  void fillStatusRecTest4(SimStatusRec* s);                   ///< unimplemented
  void fillStatusRecTest5(SimStatusRec* s);                   ///< unimplemented
  void fillStatusRecTest6(SimStatusRec* s);                   ///< unimplemented
  void fillStatusRecTest7(SimStatusRec* s);                   ///< unimplemented
  void fillStatusRecTest8(SimStatusRec* s);                   ///< unimplemented
  void fillStatusRecTest9(SimStatusRec* s);                   ///< unimplemented

  void setStatusObject(SimStatusRec*);         ///< set new data exchange object
  SimStatusRec* getStatusObject(void);            ///< get current status object

  /// Tbd.
  void  setReportFileName(const std::string&);
  /// Tbd.
  const std::string& getReportFileName(void) const;

  /*! \brief set pointer to bool variable as flag to use
  *
  *  The parameter pointer to a boolean type object is checked during every
  *  evaluate operation. If the pointer is not zero and the variable
  *  it points to has a value of true current evaluate operation will be
  *  cancelled and control is given back to caller.\nUp to now it is only
  *  used in context of graphical user interfaces to allow abort operations.
  *
  *  \param abort_pointer Pointer to boolean variable type object
  */
  void setAbortBoolPointer(bool* abort_pointer);
  /// Tbd.
  void abort();
  /// Tbd.
  bool shouldAbort(void) const;

  /// Print simulation data to csv-file
  void PrintCSV  (const char * filename);
  /// Print ACT problem-specific data to csv-file
  void writeACTCsv(const std::string& file_name);
  /// write GTOC2-problem specific data to file
  void writeGTOC2Csv(const std::string& file_name);

  /// use Dynamic Accuracy Control (DAC)
  bool useDAC(void) const;

  /// max. DAC factor
  double getDACMaxFactor(void) const;
  /// min. DAC factor
  double getDACMinFactor(void) const;
  /// max. DAC hypercube size
  double getDACMaxHcSize(void) const;
  /// min. DAC hypercube size
  double getDACFinalHcSize(void) const;
  /// set new dynamic accuracy factor
  void setAccuracyFactor(double factor);
  /// current accuracy factor
  double getAccuracyFactor(void) const;

  /// set simulation type
  void setSimType(const SimulationType type);
  /// current simulation type
  SimulationType getSimType(void) const;

  /// maximum allowed relative distance to target to fulfill boundary conditions
  double getMaxRelDistance(void) const;

  /// maximum allowed relative velocity to target to fulfill boundary conditions
  double getMaxRelVelocity(void) const;

  /// return true, if current simulation is heliocentric and false otherwise
  bool isHeliocentric(void) const;
  /// return true, if current simulation is bodycentric and false otherwise
  bool isBodyCentric(void) const;

  /// status of Dynamic Step Size Control (DSSC), [true/false]
  bool getStatusDSSC(void) const;

  void setStatusDSSC(const bool status);                 ///< switch on/off DSSC
  void enableDSSC(void);                   ///< enable Dynamic Step Size Control
  void disableDSSC(void);                 ///< disable Dynamic Step Size Control

  /// status of max. angle control (part of DSSC)
  bool getStatusDSSCStepAngleControl(void) const;
  /// switch on/off max. angle control (part of DSSC)
  void setStatusDSSCStepAngleControl(const bool status);
  void enableDSSCStepAngleControl(void);  ///< enable max. orbtial angle control
  void disableDSSCStepAngleControl(void);///< disable max. orbtial angle control

  /// set maximum orbital angle during control step
  void setDSSCMaxStepAngle(const double angle);
  /// get maximum orbital angle during control step
  double getDSSCMaxStepAngle(void) const;

  /// status of max. step distance control (part of DSSC)
  bool getStatusDSSCStepDistanceControl(void) const;
  /// switch on/off max. step distance control
  void setStatusDSSCStepDistanceControl(const bool status);
  /// switch on max. step distance control
  void enableDSSCStepDistanceControl(void);
  /// switch off max. step distance control
  void disableDSSCStepDistanceControl(void);
  /// set max. allowed step distance during control step
  void setDSSCMaxStepDistance(const double distance);
  /// maxmimum allowed step distance during control step
  double getDSSCMaxStepDistance(void) const;

  /// status of approach control (part of DSSC) [true/false]
  bool getStatusDSSCApproachControl(void) const;
  /// switch on/off approach control
  void setStatusDSSCApproachControl(const bool status);
  /// switch on approach control
  void enableDSSCApproachControl(void);
  /// switch off apporach control
  void disableDSSCApproachControl(void);
  /// set conversion factor for conversion of target distance into step size [d]
  void setDSSCApproachFactor(const double factor);
  /// get conversion factor
  double getDSSCApproachFactor(void) const;

  /// get current setting of mininum allowed control step size
  double getDSSCMinStepSize(void) const;
  /// set new minimum allowed control step size
  void setDSSCMinStepSize(const double step_size);
  /// get current setting of maximum allowed control step size
  double getDSSCMaxStepSize(void) const;
  /// set new maximum allowed control step size
  void setDSSCMaxStepSize(const double step_size);

  /// Tbd.
  bool getDynamicIntegrationIntervalReduction(void) const;
  /// Tbd.
  void setDynamicIntegrationIntervalReduction(const bool);

  /// Tbd.
  double getStartTime(void) const;
  /// Tbd.
  void setStartTime(const double);

  /// Tbd.
  double getMaxFlightTime(void) const;
  /// current flight time
  double getFlightTime(void) const;
  /// best flight time ever
  double getBestFlightTime(void) const;
  /// remaining time until max. allowed flight time
  double getRemainingFlightTime(void) const;
  /// remaining time of best flight until max. allowed flight time
  double getBestRemainingFlightTime(void) const;

  /// proximity object of best step (const version)
  const CProximity& getBestProximity(void) const;
  /// initial proximity object (const version)
  const CProximity& getInitProximity(void) const;
  /// current proximity object (const version)
  const CProximity& getCurrProximity(void) const;
  /// current proximity object (non-const version)
  CProximity& getCurrProximity(void);


  /// Tbd.
  void resetCurrentStepCounter(void);
  /// Tbd.
  int getCurrentStepCount(void) const;
  /// Tbd.
  void incCurrentStepCounter(void);

  /// Tbd.
  int getMaxStepsCount(void) const;
  /// Tbd.
  void setMaxStepsCount(const int);

  /// Tbd.
  int getMinOutputPoints(void) const;
  /// Tbd.
  void setMinOutputPoints(const int n_min_points);

  /// Tbd.
  StateType getInitialStateType(void) const;
  /// Tbd.
  void setInitialStateType(const StateType state_type);

  /// get identifier of target body
  BodyID getInitialBody(void) const;
  /// set new identifier of target body
  void setInitialBody(const BodyID body_id);

  /// Tbd.
  double getMinLaunchVelocity(void) const;
  /// Tbd.
  void setMinLaunchVelocity(const double vel);

  /// Tbd.
  double getMaxLaunchVelocity(void) const;
  /// Tbd.
  void setMaxLaunchVelocity(const double vel);

  /// Tbd.
  double getMinLaunchAzimuthAngle(void) const;
  /// Tbd.
  void setMinLaunchAzimuthAngle(const double angle);

  /// Tbd.
  double getMaxLaunchAzimuthAngle(void) const;
  /// Tbd.
  void setMaxLaunchAzimuthAngle(const double angle);

  /// Tbd.
  double getMinLaunchElevationAngle(void) const;
  /// Tbd.
  void setMinLaunchElevationAngle(const double angle);

  /// Tbd.
  double getMaxLaunchElevationAngle(void) const;
  /// Tbd.
  void setMaxLaunchElevationAngle(const double angle);

  /// Tbd.
  double getMinStartTime(void) const;
  /// Tbd.
  void setMinStartTime(const double val);
  /// Tbd.
  double getMaxStartTime(void) const;
  /// Tbd.
  void setMaxStartTime(const double val);
  /// Tbd.
  double getIntegrationInterval(void) const;
  /// Tbd.
  void setIntegrationInterval(const double val);

  /// Tbd.
  TrajectoryType getTrajectoryType(void) const;
  /// Tbd.
  OptimizationGoalType getOptGoalType(void) const;

  /// Tbd.
  StateType getTargetStateType(void) const;
  /// Tbd.
  void setTargetStateType(const StateType);

  /// Tbd.
  const std::valarray<double>& getTOEElementDeviation(void) const;
  /// Tbd.
  std::valarray<double>& getTOEElementDeviation(void);
  /// Tbd.
  const std::valarray<double>& getTOEBestElementDeviation(void) const;
  /// Tbd.
  void setTOEBestElementDeviation(const std::valarray<double>&);
  /// Tbd.
  std::valarray<double>& getTOERemChangeTime(void);
  /// Tbd.
  std::valarray<double>& getTOEBestRemChangeTime(void);

  /// Tbd.
  void setTOEBestRemChangeTime(const std::valarray<double>&);

  /// Tbd.
  const std::valarray<double>& getTargetOrbitalElementProxWeights(void) const;

  /// Tbd.
  std::valarray<double>& getTargetOrbitalElementProxWeights(void);

  /// Tbd.
  bool useTargetOrbitalElement(const SpaceLib::KEPLERIAN_ELEMENT element) const;

  /// Tbd.
  DESSolverType getIntegratorType(void) const;
  /// Tbd.
  void setIntegratorType(const DESSolverType);

  /// Tbd.
  double getMaxRelIntegrationError(void) const;
  /// Tbd.
  void setMaxRelIntegrationError(const double rel_error);

  /// Tbd.
  double getMaxAbsIntegrationError(void) const;
  /// Tbd.
  void setMaxAbsIntegrationError(const double abs_error);

  /// Tbd.
  double getMinSolarDistance(void) const;
  /// Tbd.
  void setMinSolarDistance(const double);

  /// Tbd.
  BodyID getTargetBodyId(void) const;
  /// Tbd.
  void   setTargetBodyId(const BodyID target_id);

  /// Tbd.
  CNC* getNC2(void) const {return m_nc2;};
  /// Tbd.
  void setNC2(CNC* nc) {m_nc2=nc;};

  /// Saves a copy of actual s/c trajectory internally.
  void saveLastTrajectory(CIndividual& ind);

  /// Returns trajectory of s/c's target
  Trajectory getTrajectoryTargetData(void);

  /// Returns current trajectory expressed in Mathlib::Vectors
  Trajectory getLastTrajectoryData(void);

  /// Returns trajectory of target bofy of last saved s/c trajectory.
  Trajectory getLastTrajectoryTargetData(void);

  /// Clock angle data of last saved s/c trajectory.
  DataArray getLastTrajectoryClockAngleData(void);

  /// Cone angle data of last saved s/c trajectory.
  DataArray getLastTrajectoryConeAngleData(void);

  /// Tbd.
  DataArray getLastTrajectoryThrustData(void);

  /*!\brief Evaluation function called by evolutionary algorithm.
  * \param solution generated solution chromosome to evaluate
  * \param final final evaluation (with min_output_points)
  * \param time*/
  CEAFitness* evaluate(
    CIndividual& solution, const bool final = false
    , const bool monitor = false, const double time = 0.0);

  //< end of implementation of current InTrance problems



  //------------------------------------------------------------------------------
  // steering angles approximation problem
  //------------------------------------------------------------------------------
private:

  double             m_last_saved_start_time;
  double             m_last_saved_start_velocity;
  double             m_last_saved_start_azimuth;
  double             m_last_saved_start_elevation;

  CDataArray         m_ref_time;
  CDataArray         m_ref_clock_angle;
  CDataArray         m_ref_cone_angle;
  CDataArray         m_ref_thrust;

  CDataArray m_last_saved_time_steps;
  CDataArray m_last_clock_angle_approx_time_steps;
  CDataArray m_last_cone_angle_approx_time_steps;

  ///< last clock angle approximation and it's derivative histories
  CDataArray  m_last_saved_clock_angle_approx_values;
  CDataArray  m_last_saved_clock_angle_approx_values_dot;

  ///< last cone angle approximation and it's derivative histories
  CDataArray  m_last_saved_cone_angle_approx_values;
  CDataArray  m_last_saved_cone_angle_approx_values_dot;

#if defined COMMONCPP
  ost::Mutex m_last_saved_time_steps_lock;
  ost::Mutex m_last_saved_clock_angle_approx_values_dot_lock;
  ost::Mutex m_last_saved_cone_angle_approx_values_dot_lock;
  ost::Mutex m_last_saved_trajectory_approx_lock;
  ost::Mutex m_last_saved_clock_angle_approx_values_lock;
  ost::Mutex m_last_saved_cone_angle_approx_values_lock;
#endif //COMMONCPP

  ///< last approximated s/c trajectories vector
  SCTraj      m_last_saved_trajectory_approx;
  uint  m_last_saved_trajectory_approx_steps;

  ///< spline knot data values of best clock angle history
  std::vector<double> m_best_clock_knot_values;
  std::vector<double> m_best_clock_knot_values_d;
  std::vector<double> m_load_clock_knot_values;

  ///< spline knot time values of best clock angle history
  std::vector<double> m_best_clock_knot_time;
  std::vector<double> m_load_clock_knot_time;

  ///< spline knot time values of best cone angle history
  std::vector<double> m_best_cone_knot_time;
  std::vector<double> m_load_cone_knot_time;

  ///< spline knot data values of best cone angle history
  std::vector<double> m_best_cone_knot_values;
  std::vector<double> m_best_cone_knot_values_d;
  std::vector<double> m_load_cone_knot_values;

  ///< actual index (for write operations) in all data structures concerning
  ///< s/c data due to steering angles approximation
  unsigned int m_approx_put_index;
  unsigned int m_approx_get_index;

  double m_max_clock_jump;   ///< maximum allowed input clock angle difference
  double m_approx_res_factor;        ///< resolution factor of approximation

  CInTranceFitness m_approx_fitness; ///< approximation fitness

  bool m_clock_angle_history;
  bool m_clock_angle_history_approximation;
  bool m_clock_angle_history_near_pi;
  bool m_cone_angle_history;
  bool m_cone_angle_history_approximation;
  bool m_thrust_history;

  int                m_n_approx_sc_traj;

  int  m_clock_angle_history_size;
  int  m_clock_angle_approx_seg_size;
  int  m_clock_angle_approx_poly_degree;
  int  m_cone_angle_approx_seg_size;
  int  m_cone_angle_approx_poly_degree;

  // private methods  

  /// Tbd.
  void getBasicClockAngleData(const SCTraj& t, CDataArray& a, const uint n_max_points = 0);

  /// Tbd.
  void saveApproxParams(
    const std::vector<double>& kt ///< knot time
    , const std::vector<double>& kv ///< knot value
    , const std::string& filename     ///< filename
    , unsigned int sz                 ///< # of values
    ) const;

  /// Tbd.
  void getBasicConeAngleData(const SCTraj& t, CDataArray& a
    , const uint n_max_points = 0);

  /// Tbd.
  CEAFitness* evaluateApproximationSolution(
    CIndividual& solution
    ,const bool final = false
    ,const bool monitor = false
    ,const double time = 0.0);

  /// Sets internal data structures with given data.
  template<class T, class U> void setData(const T& s, U& d);

  /// Tbd.
  void calcLastTrajectoryClockAngleApproximation(void);

  /// Tbd.
  void calcLastTrajectoryConeAngleApproximation(void);

  /*!
  * \brief pointer to vector of trajectory state information
  * \param source trajectory from which the state data should be extracted
  * \return pointer to CTrajectory class object containing the desired data
  */
  Trajectory getTrajectoryData(const SCTraj& source) const;

  /// Tbd.
  template<class T, class U>
  CEAFitness* evaluateSteeringAngleApprox(
    CIndividual& solution
    ,const U&    clock_angle
    ,const U&    cone_angle
    ,const T&    time
    ,SCTraj&        traj
    ,unsigned int& n_steps
    ,double start_t=0.0
    ,bool final  = false
    ,bool monitor = false
    );

  /*!\brief Returns a copy of best trajectory so far (m_steps_best values).
  * \param Vector SCTraj to be set with copy of current best trajectory.*/
  void copyCurrentTrajectory(SCTraj&);

  /// return reference of current sc trajectory
  SCTraj& getCurrentTrajectory(void);

  /// save approximation knot time values (clock angle)
  void saveApproximationClockTime(const std::string& filename) const;

  /// save approximation knot data values (clock angle)
  void saveApproximationClockValues(const std::string& filename) const;

  /// save approximation knot time values (cone angle)
  void saveApproximationConeTime(const std::string& filename) const;

  /// save approximation knot data values (cone angle)
  void saveApproximationConeValues(const std::string& filename) const;

  /// Tbd.
  void incp(void);

  /// Tbd.
  void decp(void);

  /// Tbd.
  void incg(void);

  /// Tbd.
  void decg(void);

  /// Tbd.
  unsigned int tellp(void);

  /// Tbd.
  unsigned int seekp(unsigned int i);

  /// Tbd.
  unsigned int tellg(void);

  /// Tbd.
  unsigned int seekg(unsigned int i);

  /// Sets reference trajectory with given trajectory.
  void setReferenceTrajectory(const SCTraj& t);
  const SCTraj& getRefTraj(void); ///< Returns current reference trajectory.

  /*!\brief Sets reference trajectory with given trajectory (data swapping).
  * \param t New reference trajectory data. After function call it'll contain
  * the former reference trajectory*/
  void swapReferenceTrajectory(SCTraj& t);

  void setReferenceTrajectoryTime(const std::vector<double>& data); ///< Tbd.
  void setReferenceTrajectoryTime(const SCTraj& traj); ///< Tbd.
  const std::vector<double>& getReferenceTrajectoryTime(void) const; ///< Tbd.
  void setReferenceTrajectoryClockAngleData(const SCTraj& traj); /// < Tbd.

  /// Tbd.
  void setReferenceTrajectoryConeAngleData(const std::vector<double>& data);
  void setReferenceTrajectoryConeAngleData(const SCTraj& traj); ///< Tbd.
  /// Tbd.
  const std::vector<double>& getReferenceTrajectoryConeAngleData(void) const;
  const std::vector<double>& getReferenceTrajectoryThrustData(void) const;

  /// Returns current reference trajectory.
  const SCTraj& getReferenceTrajectory(void);

  /*! \brief Copy of last saved s/c trajectory.
  * \param traj Reference of trajectory object to be filled with data from
  * last saved s/c trajectory.*/
  void copyLastSavedTrajectory(SCTraj& traj);

  /// Tbd.
  void calcLastSavedTrajectoryApprox(void);

  /*!\brief Make a copy of last s/c trajectory trapproximation.
  * \param traj Data structure to hold s/c trajectory*/
  void copyLastSavedTrajectoryApprox(SCTraj& traj);


  /// Calculates approximation of given input data
  void calcApproximation(
    const CDataArray& data
    ,const CDataArray& time
    ,const unsigned int degree
    ,const unsigned int seg_size
    ,const double resolution_factor
    ,      CDataArray& approx_time   
    ,      CDataArray& approx_data
    ,      CDataArray& aprrox_data_derivative
    );

  /// Sets new clock angle data of last saved s/c trajectory.
  void setLastTrajectoryTimeSteps(const valarray<double>& data);

  /// Clock angle data set of actual s/c trajectory.
  DataArray getClockAngleData(void);

  /*!
  * \brief Calculates new approximation of clock angle and stores it in object.
  *
  * \param data
  * \param time
  * \param new_time
  * \param new_data
  * \param new_data_d
  * \param degree
  * \param n_seg_points
  * \return
  */
  void calcClockAngleApproximation(
    const CDataArray& data
    ,const CDataArray& time
    ,      CDataArray& new_time
    ,      CDataArray& new_data
    ,      CDataArray& new_data_d
    ,const int         degree
    ,const int         n_seg_points);


  /// Returns a polynomial approximation of clock angle history.
  DataArray getClockAngleApproximation(void);

  /// Sets new clock angle data of last saved s/c trajectory.
  template<class T>
  void setLastTrajectoryClockAngleData(const T& data);

  /// Sets new clock angle approximation data of last 
  /// saved s/c trajectory.
  template<class T>
  void setLastTrajectoryClockAngleApproxData(const T& data);

  /// Sets new derivative of clock angle approximation data of 
  /// last saved s/c trajectory.
  template<class T>
  void setLastTrajectoryClockAngleApproxDataDerivative(const T& data_d);

  /*!\brief set data and first derivative of last approximated clock angles
  *
  * \param d data values structure
  * \param d_d first derivative of structure's data values
  */
  template<class T, class U>
  void setLastTrajectoryClockAngleApproxData(const T& d, const U& d_d);

  /// Cone angle data set of actual s/c trajectory.
  DataArray getConeAngleData(void);

  /// Sets new clock angle data of last saved s/c trajectory.
  template<class T>
  void setLastTrajectoryConeAngleData(const T& data);

  /// Sets new derivative of cone angle approximation data of 
  /// last saved s/c trajectory.
  template<class T>
  void setLastTrajectoryConeAngleApproxDataDerivative(const T& data_d);

  /*!\brief set data and first derivative of last approximated cone angles
  * \param d data values structure
  * \param d_d first derivative of structure's data values*/
  template<class T, class U>
  void setLastTrajectoryConeAngleApproxData(const T& d, const U& d_d);

  ///  Calculates new approximation of cone angle and stores it in object.
  void calcConeAngleApproximation(
    const CDataArray& data         ///<
    ,const CDataArray& time         ///<
    ,      CDataArray& new_time     ///<
    ,      CDataArray& new_data     ///<
    ,      CDataArray& new_data_d   ///<
    ,const int         degree       ///<
    ,const int         n_seg_points ///<
    ,const double      upper_bound  ///<
    );

  /// Returns a copy ofpolynomial approximation of one angle history.
  DataArray getConeAngleApproximation(void);

public:

  /// activate clock angle history data extraction
  void activateClockAngleHistory(void);
  /// deactivate clock angle history data extraction
  void deactivateClockAngleHistory(void);
  /// activate clock angle history approximation data extraction
  void activateClockAngleApproximation(void);
  /// deactivate clock angle history approximation data extraction
  void deactivateClockAngleApproximation(void);
  /// set clock angle approximation parameters
  void setClockAngleApproximationParameters(const int segSize
    , const int polyDegree);
  /// activate cone angle history data extraction
  void activateConeAngleHistory(void);
  /// deactivate cone angle history data extraction
  void deactivateConeAngleHistory(void);
  /// activate cone angle history approximation
  void activateConeAngleApproximation(void);
  /// deactivate cone angle history approximation
  void deactivateConeAngleApproximation(void);
  /// set cone angle approximation parameters
  void setConeAngleApproximationParameters(const int segSize
    , const int polyDegree);
  /// activate thrust history data extraction
  void activateThrustHistory(void);
  /// deactivate thrust history data extraction
  void deactivateThrustHistory(void);
  /// Tbd.
  DataArray getLastTrajectoryClockAngleApproximation(void);
  /// Returns current approximated trajectory expressed in Mathlib::Vectors
  Trajectory getLastTrajectoryApproxData(void);
  /// Sets new clock angle approximation data of last saved s/c trajectory.
  /// Returns cone angle approximation of last s/c trajectory
  DataArray getLastTrajectoryConeAngleApproximation(void);
  template<class T>
  void setLastTrajectoryConeAngleApproxData(const T& data);
  Trajectory getReferenceTrajectoryData(void); ///< Tbd.
  Trajectory getReferenceTrajectoryTargetData(void); ///< Tbd.
  DataArray getReferenceTrajectoryClockAngleData(void); ///< Tbd.
  DataArray getReferenceTrajectoryConeAngleData(void); ///< Tbd.
  DataArray getReferenceTrajectoryThrustData(void); ///< Tbd.
  /*!\brief Sets last computed trajectory (only m_steps_best values) as reference
  *trajectory.*/
  void setReferenceTrajectory(void);
  void extractDataFromReferenceTrajectory(void); ///< Tbd.
  const std::vector<double>& getRefTrajClockAngleArray(void) const; ///< Tbd.
  /// Tbd.
  template<class T>
  void postprocessAngleData(T& a);
  /// Tbd.
  void setReferenceTrajectoryClockAngleData(const std::vector<double>& data);
  /// wrapper funtion for saving last best approximation of clock angles
  void saveClockApproxParams(const std::string& filename) const;
  /// wrapper funtion for saving last best approximation of cone angles
  void CSimulation::saveConeApproxParams(const std::string& filename) const;

  // end of steering angles approximation part

}; // class CSimulation


//------------------------------------------------------------------------------
// implementation of inline methods
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
inline bool CSimulation::getDynamicIntegrationIntervalReduction(void) const
{
  return m_data.getDynamicIntegrationIntervalReduction();
} // bool CSimulation::getDynamicIntegrationIntervalReduction(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setDynamicIntegrationIntervalReduction(const bool flag)
{
  m_data.setDynamicIntegrationIntervalReduction(flag);
} // void CSimulation::setDynamicIntegrationIntervalReduction(const bool)


//------------------------------------------------------------------------------
inline double CSimulation::getStartTime(void) const
{
  return m_time_start;
} // double CSimulation::getStartTime(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setStartTime(const double start_time)
{
  m_time_start = start_time;
} // void CSimulation::setStartTime(const double)


//------------------------------------------------------------------------------
inline double CSimulation::getMaxFlightTime(void) const
{
  return getIntegrationInterval();
} // double CSimulation::getMaxFlightTime(void) const


//------------------------------------------------------------------------------
inline double CSimulation::getFlightTime(void) const
{
  return (m_t-getStartTime());
} // double CSimulation::getFlightTime(void) const


//------------------------------------------------------------------------------
inline double CSimulation::getBestFlightTime(void) const
{
  return m_time_best;
} // double getBestFlightTime(void) const


//------------------------------------------------------------------------------
inline const CProximity& CSimulation::getInitProximity(void) const
{
  return m_init_proximity;
} // const CProximity& CSimulation::getCurrProximity(void) const


//------------------------------------------------------------------------------
inline const CProximity& CSimulation::getBestProximity(void) const
{
  return m_best_proximity;
} // const CProximity& CSimulation::getCurrProximity(void) const


//------------------------------------------------------------------------------
inline CProximity& CSimulation::getCurrProximity(void)
{
  return m_curr_proximity;
} // const CProximity& CSimulation::getCurrProximity(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setBestProximity(CProximity& prox)
{
  m_best_proximity = prox;
} // void setBestProximity(CProximity& prox)


//------------------------------------------------------------------------------
inline void CSimulation::setInitProximity(CProximity& prox)
{
  m_init_proximity = prox;
} // void setBestProximity(CProximity& prox)


//------------------------------------------------------------------------------
inline double CSimulation::getRemainingFlightTime(void) const
{
  return ( getIntegrationInterval() - getFlightTime() );
} // double CSimulation::getRemainingFlightTime(void) const


//------------------------------------------------------------------------------
inline double CSimulation::getBestRemainingFlightTime(void) const
{
  return ( getIntegrationInterval() - getBestFlightTime());
} // double getBestRemainingFlightTime(void) const


//------------------------------------------------------------------------------
inline void CSimulation::resetCurrentStepCounter(void)
{
  m_steps = 0;
} // void CSimulation::resetCurrentStepCounter(void)


//------------------------------------------------------------------------------
inline int CSimulation::getCurrentStepCount(void) const
{
  return m_steps;
} // int CSimulation::getCurrentStepCount(void) const


//------------------------------------------------------------------------------
inline void CSimulation::incCurrentStepCounter(void)
{
  ++m_steps;
} // void CSimulation::incStepCounter(void)


//------------------------------------------------------------------------------
inline int CSimulation::getMaxStepsCount(void) const
{
  return m_data.getMaxStepsCount();
} // int CSimulation::getMaxStepsCount(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setMaxStepsCount(const int count)
{
  m_data.setMaxStepsCount(count);
} // void CSimulation::setMaxStepsCount(const int)


//------------------------------------------------------------------------------
inline StateType CSimulation::getInitialStateType(void) const
{
  return m_data.getInitialStateType();
} // StateType CSimulation::getInitialStateType(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setInitialStateType(const StateType state_type)
{
  m_data.setInitialStateType(state_type);
} // void CSimulation::setInitialStateType(const StateType state_type)


//------------------------------------------------------------------------------
inline BodyID CSimulation::getInitialBody(void) const
{
  return m_data.getInitialBody();
} // BodyID CSimulation::getInitialBody(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setInitialBody(const BodyID body_id)
{
  m_data.setInitialBody(body_id);
} // void CSimulation::setInitialBody(const BodyID body_id)


//------------------------------------------------------------------------------
inline int CSimulation::getMinOutputPoints(void) const
{
  return m_data.getMinOutputPoints();
} // int CSimulation::getMinOutputPoints(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setMinOutputPoints(const int n_min_points)
{
  m_data.setMinOutputPoints(n_min_points);
} // void CSimulation::setMinOutputPoints(const int n_min_points)


//------------------------------------------------------------------------------
inline double CSimulation::getMinLaunchVelocity(void) const
{
  return m_data.getMinLaunchVelocity();
} // double getMinStartVelocity(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setMinLaunchVelocity(const double vel)
{
  m_data.setMinLaunchVelocity(vel);
} // CSimulation::setMinLaunchVelocity(const double vel)


//------------------------------------------------------------------------------
inline double CSimulation::getMaxLaunchVelocity(void) const
{
  return m_data.getMaxLaunchVelocity();
} // double getMaxStartVelocity(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setMaxLaunchVelocity(const double vel)
{
  m_data.setMaxLaunchVelocity(vel);
} // CSimulation::setMaxLaunchVelocity(const double vel)


//------------------------------------------------------------------------------
inline double CSimulation::getMinLaunchAzimuthAngle(void) const
{
  return m_data.getMinLaunchAzimuthAngle();
} // double CSimulation::getMinLaunchAzimuthAngle(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setMinLaunchAzimuthAngle(const double angle)
{
  m_data.setMinLaunchAzimuthAngle(angle);
} // void CSimulation::setMinLaunchAzimuthAngle(const double angle)


//------------------------------------------------------------------------------
inline double CSimulation::getMaxLaunchAzimuthAngle(void) const
{
  return m_data.getMaxLaunchAzimuthAngle();
} // double CSimulation::getMaxLaunchAzimuthAngle(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setMaxLaunchAzimuthAngle(const double angle)
{
  m_data.setMaxLaunchAzimuthAngle(angle);
} // void CSimulation::setMaxLaunchAzimuthAngle(const double angle)


//------------------------------------------------------------------------------
inline double CSimulation::getMinLaunchElevationAngle(void) const
{
  return m_data.getMinLaunchElevationAngle();
} // double CSimulation::getMinLaunchElevationAngle(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setMinLaunchElevationAngle(const double angle)
{
  m_data.setMinLaunchElevationAngle(angle);
} // void CSimulation::setMinLaunchElevationAngle(const double angle)


//------------------------------------------------------------------------------
inline double CSimulation::getMaxLaunchElevationAngle(void) const
{
  return m_data.getMaxLaunchElevationAngle();
} // double CSimulation::getMaxLaunchElevationAngle(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setMaxLaunchElevationAngle(const double angle)
{
  m_data.setMaxLaunchElevationAngle(angle);
} // void CSimulation::setMaxLaunchElevationAngle(const double angle)


//------------------------------------------------------------------------------
inline BodyID CSimulation::getTargetBodyId(void) const
{
  return m_data.getTargetBodyId();
} // BodyID CSimulation::TargetBodyId(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setTargetBodyId(const BodyID target_id)
{
  m_data.setTargetBodyId(target_id);
}


//------------------------------------------------------------------------------
inline void CSimulation::printMessage(const std::ostringstream& m)
{
  statusMessage(m);
} // void CSimulation::printMessage(const std::ostringstream& m)


//------------------------------------------------------------------------------
inline void CSimulation::writeMessageToReportFile(const std::ostringstream& mess)
{
  if (m_rep.size() != 0)
  {
    ofstream report_file(m_rep.c_str(), ios::app);
    if (report_file.is_open() == true) report_file << mess.str();
  }
} // void CSimulation::writeMessageToReportFile(const std::ostringstream& mess)


//------------------------------------------------------------------------------
inline void CSimulation::writeMessageToReportFile(const std::string& mess)
{
  if (m_rep.size() != 0)
  {
    ofstream report_file(m_rep.c_str(), ios::app);
    if (report_file.is_open() == true) report_file << mess;
  }
} // void CSimulation::writeMessageToReportFile(const std::string& mess)


//------------------------------------------------------------------------------
inline void CSimulation::fillStatusRec(void) 
{
  fillStatusRec(m_status);
} // void CSimulation::fillStatusRec(void)


//------------------------------------------------------------------------------
inline void CSimulation::setStatusObject(SimStatusRec* p) 
{
  m_status = p;
} // void CSimulation::setStatusObject(SimStatusRec* p)


//------------------------------------------------------------------------------
inline SimStatusRec* CSimulation::getStatusObject(void)
{
  return m_status;
} // SimStatusRec* CSimulation::getStatusObject(void)


//------------------------------------------------------------------------------
inline void CSimulation::setReportFileName(const std::string& name)
{
  m_rep = name;
} // void CSimulation::setReportFileName(const std::string& name)


//------------------------------------------------------------------------------
inline const std::string& CSimulation::getReportFileName(void) const 
{ 
  return m_rep;
}


//------------------------------------------------------------------------------
inline void CSimulation::setAbortBoolPointer(bool* abort_pointer) 
{
  m_should_abort = abort_pointer;
} // void CSimulation::setAbortBoolPointer(bool* abort_pointer)


//------------------------------------------------------------------------------
inline void CSimulation::abort()
{
  if (m_should_abort!=0) (*m_should_abort)=true;
} // inline void CSimulation::abort()


//------------------------------------------------------------------------------
inline bool CSimulation::shouldAbort(void) const
{
  return (m_should_abort!=0)?(*m_should_abort):false;
} // bool CSimulation::shouldAbort(void) const


//------------------------------------------------------------------------------
inline bool CSimulation::useDAC(void) const
{
  return m_data.useDynamicAccuracyControl();
} // bool CSimulation::useDAC(void) const 


//------------------------------------------------------------------------------
inline double CSimulation::getDACMaxFactor(void) const 
{
  return m_data.getDACMaxFactor();
} // double CSimulation::getDACMaxFactor(void) const 


//------------------------------------------------------------------------------
inline double CSimulation::getDACMinFactor(void) const 
{
  return double(1.0);
} // double CSimulation::getDACMinFactor(void) const


//------------------------------------------------------------------------------
inline double CSimulation::getDACMaxHcSize(void) const
{
  return m_data.getDACMaxHypercubeSize();
} // double CSimulation::getDACMaxHcSize(void) const 


//------------------------------------------------------------------------------
inline double CSimulation::getDACFinalHcSize(void) const
{
  return m_data.getDACMinHypercubeSize();
} // double CSimulation::getDACFinalHcSize(void) const


//------------------------------------------------------------------------------
inline double CSimulation::getAccuracyFactor(void) const
{
  return m_proposed_accuracy_factor;
} // double getAccuracyFactor(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setAccuracyFactor(double factor)
{
  m_proposed_accuracy_factor = factor;
} // void CSimulation::setAccuracyFactor(double factor)


//------------------------------------------------------------------------------
inline void CSimulation::setSimType(const SimulationType type) 
{
  m_sim_type = type;
} // void CSimulation::setSimType(const SimulationType type)


//------------------------------------------------------------------------------
inline SimulationType CSimulation::getSimType(void) const
{
  return m_sim_type;
} // SimulationType CSimulation::getSimType(void) const


//------------------------------------------------------------------------------
inline double CSimulation::getMaxRelDistance(void) const
{
  return m_max_dr_AU;
} // double getMaxRelDistance(void) const

//------------------------------------------------------------------------------
inline double CSimulation::getMaxRelVelocity(void) const
{
  return m_max_dv_AUd;
} // double getMaxRelVelocity(void) const


//------------------------------------------------------------------------------
inline bool CSimulation::getStatusDSSC(void) const
{
  return m_data.getStatusDSSC();
} // bool CSimulation::getStatusDSSC(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setStatusDSSC(const bool status)
{
  m_data.setStatusDSSC(status);
} // void CSimulation::setStatusDSSC(const bool status)


//------------------------------------------------------------------------------
inline void CSimulation::enableDSSC(void)
{
  m_data.enableDSSC();
} // void CSimulation::enableDSSC(void)


//------------------------------------------------------------------------------
inline void CSimulation::disableDSSC(void)
{
  m_data.disableDSSC();
} // void CSimulation::disableDSSC(void)


//------------------------------------------------------------------------------
inline bool CSimulation::getStatusDSSCStepAngleControl(void) const
{
  return m_data.getStatusDSSCStepAngleControl();
} // bool CSimulation::getStatusDSSCStepAngleControl(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setStatusDSSCStepAngleControl(const bool status)
{
  m_data.setStatusDSSCStepAngleControl(status);
} // void CSimulation::setStatusDSSCStepAngleControl(const double status) const


//------------------------------------------------------------------------------
inline void CSimulation::enableDSSCStepAngleControl(void)
{
  m_data.enableDSSCStepAngleControl();
} // void CSimulation::enableDSSCStepAngleControl(void)


//------------------------------------------------------------------------------
inline void CSimulation::disableDSSCStepAngleControl(void)
{
  m_data.disableDSSCStepAngleControl();
} // void CSimulation::disableDSSCStepAngleControl(void);


//------------------------------------------------------------------------------
inline void CSimulation::setDSSCMaxStepAngle(const double angle)
{
  m_data.setDSSCMaxStepAngle(angle);
} // void CSimulation::setDSSCMaxStepAngle(const double angle)


//------------------------------------------------------------------------------
inline double CSimulation::getDSSCMaxStepAngle(void) const
{
  return m_data.getDSSCMaxStepAngle();
} // double CSimulation::getDSSCMaxStepAngle(void) const


//------------------------------------------------------------------------------
inline bool CSimulation::getStatusDSSCStepDistanceControl(void) const
{
  return m_data.getStatusDSSCStepDistanceControl();
} // bool CSimulation::getStatusDSSCStepDistanceControl(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setStatusDSSCStepDistanceControl(const bool status)
{
  m_data.setStatusDSSCStepDistanceControl(status);
} // void CSimulation::setStatusDSSCStepDistanceControl(const bool status)


//------------------------------------------------------------------------------
inline void CSimulation::enableDSSCStepDistanceControl(void)
{
  m_data.enableDSSCStepDistanceControl();
} // void CSimulation::enableDSSCStepDistanceControl(void)


//------------------------------------------------------------------------------
inline void CSimulation::disableDSSCStepDistanceControl(void)
{
  m_data.disableDSSCStepDistanceControl();
} // void CSimulation::disableDSSCStepDistanceControl(void)


//------------------------------------------------------------------------------
inline void CSimulation::setDSSCMaxStepDistance(const double distance)
{
  m_data.setDSSCMaxStepDistance(distance);
} // void CSimulation::setDSSCMaxStepDistance(const double distance)


//------------------------------------------------------------------------------
inline double CSimulation::getDSSCMaxStepDistance(void) const
{
  return m_data.getDSSCMaxStepDistance();
} // double CSimulation::getDSSCMaxStepDistance(void) const


//------------------------------------------------------------------------------
inline bool CSimulation::getStatusDSSCApproachControl(void) const
{
  return m_data.getStatusDSSCApproachControl();
} // bool CSimulation::getStatusDSSCApproachControl(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setStatusDSSCApproachControl(const bool status)
{
  m_data.setStatusDSSCApproachControl(status);
} // void CSimulation::setStatusDSSCApproachControl(const bool status)


//------------------------------------------------------------------------------
inline void CSimulation::enableDSSCApproachControl(void)
{
  m_data.enableDSSCApproachControl();
} // void CSimulation::enableDSSCApproachControl(void)


//------------------------------------------------------------------------------
inline void CSimulation::disableDSSCApproachControl(void)
{
  m_data.disableDSSCApproachControl();
} // void CSimulation::disableDSSCApproachControl(void)


//------------------------------------------------------------------------------
inline void CSimulation::setDSSCApproachFactor(const double factor)
{
  m_data.setDSSCApproachFactor(factor);
} // void CSimulation::setDSSCApproachFactor(const double factor)


//------------------------------------------------------------------------------
inline double CSimulation::getDSSCApproachFactor(void) const
{
  return m_data.getDSSCApproachFactor();
} // double CSimulation::getDSSCApproachFactor(void) const


//------------------------------------------------------------------------------
inline double CSimulation::getDSSCMinStepSize(void) const
{
  return m_data.getDSSCMinStepSize();
} // double CSimulation::getDSSCMinStepSize(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setDSSCMinStepSize(const double step_size)
{
  m_data.setDSSCMinStepSize(step_size);
} // void CSimulation::setDSSCMinStepSize(const double step_size)


//------------------------------------------------------------------------------
inline double CSimulation::getDSSCMaxStepSize(void) const
{
  return m_data.getDSSCMaxStepSize();
} // double CSimulation::getDSSCMaxStepSize(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setDSSCMaxStepSize(const double step_size)
{
  m_data.setDSSCMaxStepSize(step_size);
} // void CSimulation::setDSSCMaxStepSize(const double step_size) 


//------------------------------------------------------------------------------
inline TrajectoryType CSimulation::getTrajectoryType(void) const 
{
  return m_data.getTrajectoryType();
} // TrajectoryType CSimulation::getTrajectoryType(void) const


//------------------------------------------------------------------------------
inline OptimizationGoalType CSimulation::getOptGoalType(void) const
{
  return m_data.getOptGoalType();
} // OptimizationGoalType CSimulation::getOptGoalType(void) const


//------------------------------------------------------------------------------
inline StateType CSimulation::getTargetStateType(void) const
{
  return m_data.getTargetStateType();
} // StateType CSimulation::getTargetStateType(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setTargetStateType(const StateType target_state)
{
  m_data.setTargetStateType(target_state);
} // void CSimulation::setTargetStateType(const StateType target_state)


//------------------------------------------------------------------------------
inline const std::valarray<double>& CSimulation::getTargetOrbitalElementProxWeights(void) const
{
  return m_data.getTargetOrbitalElementProxWeights();
} // inline std::valarray<double>& getTargetOrbitalElementProxWeights(void) const


//------------------------------------------------------------------------------
inline std::valarray<double>& CSimulation::getTargetOrbitalElementProxWeights(void)
{
  return m_data.getTargetOrbitalElementProxWeights();
} // std::valarray<double>& getTargetOrbitalElementProxWeights(void)


//------------------------------------------------------------------------------
inline bool CSimulation::useTargetOrbitalElement(
  const SpaceLib::KEPLERIAN_ELEMENT element) const
{
  return m_data.useTargetOrbitalElement(element);
} // useTargetOrbitalElement(const SpaceLib::KEPLERIAN_ELEMENT element) const


//------------------------------------------------------------------------------
/// get minimum allowed start/launch date as Modified Julian Date [MJD]
inline double CSimulation::getMinStartTime(void) const
{
  return m_data.getMinStartTime();
} // double CSimulation::getMinStartTime(void) const 


/// set minimum allowed start/launch date as Modified Julian Date [MJD]
//------------------------------------------------------------------------------
inline void CSimulation::setMinStartTime(const double val) 
{ 
  m_data.setMinStartTime(val);
} // void CSimulation::setMinStartTime(const double val)


/// get maximum allowed start/launch date as Modified Julian Date [MJD]
//------------------------------------------------------------------------------
inline double CSimulation::getMaxStartTime(void) const 
{
  return m_data.getMaxStartTime();
} // double CSimulation::getMaxStartTime(void) const 


//------------------------------------------------------------------------------
/// set maximum allowed start/launch date as Modified Julian Date [MJD]
inline void CSimulation::setMaxStartTime(const double val) 
{ 
  m_data.setMaxStartTime(val);
} // void CSimulation::setMaxStartTime(const double val)


//------------------------------------------------------------------------------
/// get maximum integration interval in days [d]
inline double CSimulation::getIntegrationInterval(void) const
{
  return m_data.getIntegrationInterval();
} // double CSimulation::getIntegrationInterval(void) const


//------------------------------------------------------------------------------
/// set maximum integration interval in days [d]
inline void CSimulation::setIntegrationInterval(const double val)
{
  m_data.setIntegrationInterval(val);
} // void CSimulation::setIntegrationInterval(const double val)


//------------------------------------------------------------------------------
inline DESSolverType CSimulation::getIntegratorType(void) const
{
  return m_data.getIntegratorType();
} // DESSolverType CSimulation::getIntegratorType(void) const

//------------------------------------------------------------------------------
inline void CSimulation::setIntegratorType(const DESSolverType integrator_type)
{
  m_data.setIntegratorType(integrator_type);
} // void CSimulation::setIntegratorType(const DESSolverType integrator_type)

//------------------------------------------------------------------------------
inline double CSimulation::getMaxRelIntegrationError(void) const
{
  return m_data.getMaxRelIntegrationError();
} // double CSimulation::getMaxRelIntegrationError(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setMaxRelIntegrationError(const double rel_error)
{
  m_data.setMaxRelIntegrationError(rel_error);
} // void CSimulation::setMaxRelIntegrationError(const double rel_error)


//------------------------------------------------------------------------------
inline double CSimulation::getMaxAbsIntegrationError(void) const
{
  return m_data.getMaxAbsIntegrationError();
} // double CSimulation::getMaxAbsIntegrationError(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setMaxAbsIntegrationError(const double abs_error)
{
  m_data.setMaxAbsIntegrationError(abs_error);
} // void CSimulation::setMaxAbsIntegrationError(const double rel_error)


//------------------------------------------------------------------------------
inline double CSimulation::getMinSolarDistance(void) const
{
  return m_data.getMinSolarDistance();
} // double CSimulation::getMinSolarDistance(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setMinSolarDistance(const double distance)
{
  m_data.setMinSolarDistance(distance);
} // void CSimulation::setMinSolarDistance(const double)


//------------------------------------------------------------------------------
inline CSimulation::~CSimulation(void) 
{
  freeEOM();
  SafeDeleteArray(m_ga_body);
  if (m_solsys!=0)
  {
    m_solsys->removeBody(getCurrentSC().getId());
  }
} // CSimulation::~CSimulation(void)


//------------------------------------------------------------------------------
inline void CSimulation::getBasicThrustData(const SCTraj& t, CDataArray& a
                                            , const uint n_max_points)
{
  a.resize(n_max_points>0?(n_max_points<t.size()?n_max_points:t.size()):t.size());
  for (unsigned int i=0;i<a.size();++i) a[i]=t[i].getThrust();
} // void CSimulation::getBasicThrustData(...)


//------------------------------------------------------------------------------
inline void CSimulation::setReferenceTrajectoryConeAngleData(
  const std::vector<double>& data)
{
  m_ref_cone_angle = data;
} // void CSimulation::setReferenceTrajectoryConeAngleData(...)

//------------------------------------------------------------------------------
inline const std::vector<double>& 
CSimulation::getReferenceTrajectoryConeAngleData(void) const
{
  return m_ref_cone_angle;
} 

//------------------------------------------------------------------------------
inline const std::vector<double>&
CSimulation::getReferenceTrajectoryThrustData(void) const
{
  return m_ref_thrust;
}

//------------------------------------------------------------------------------
inline const SCTraj& CSimulation::getReferenceTrajectory(void)
{
  return m_ref_sc_traj;
} // const SCTraj& CSimulation::getReferenceTrajectory(void)

//------------------------------------------------------------------------------
inline const std::vector<double>& 
CSimulation::getRefTrajClockAngleArray(void) const
{
  return m_ref_clock_angle;
} // const vector<double>& getRefTrajClockAngleArray(void) const

//------------------------------------------------------------------------------
inline void 
CSimulation::setReferenceTrajectoryClockAngleData(const std::vector<double>& d)
{
  m_ref_clock_angle = d;
} // void setReferenceTrajectoryClockAngleData(const std::vector<double>&)

//------------------------------------------------------------------------------
inline const std::vector<double>& 
CSimulation::getReferenceTrajectoryTime(void) const
{
  return m_ref_time;
} // const vector<double>& getReferenceTrajectoryTime(void) const

//------------------------------------------------------------------------------
inline void 
CSimulation::setReferenceTrajectoryTime(const std::vector<double>& d) 
{
  m_ref_time = d;
} // void setRefrenceTrajectoryTime(const vector<double>&)

//------------------------------------------------------------------------------
inline const SCTraj& CSimulation::getRefTraj(void)
{
  return m_ref_sc_traj;
} // const SCTraj& getRefTraj(void)

//------------------------------------------------------------------------------
inline void CSimulation::setReferenceTrajectory(void)
{
  copyCurrentTrajectory(m_ref_sc_traj);
} // void CSimulation::setReferenceTrajectory(void)

//------------------------------------------------------------------------------
inline SCTraj& CSimulation::getCurrentTrajectory(void)
{
  return m_spacecraft_trajectory;
} // SCTraj& CSimulation::getCurrentTrajectory(void)


//------------------------------------------------------------------------------
inline void CSimulation::setLastTrajectoryTimeSteps(const valarray<double>& data)
{
#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_time_steps_lock);
#endif //COMMONCPP
  setData(data, m_last_saved_time_steps);
} // void setLastTrajectoryTimeSteps(const valarray<double>& data)


//------------------------------------------------------------------------------
template<class T>
inline void CSimulation::setLastTrajectoryClockAngleData(const T& data)
{
#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_clock_angle_history_lock);
#endif //COMMONCPP
  setData(data, m_last_saved_clock_angle_history);
} // void CSimulation::setLastTrajectoryClockAngleData(const T& data)


//------------------------------------------------------------------------------
template<class T>
inline void CSimulation::setLastTrajectoryClockAngleApproxData(const T& data)
{
#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_clock_angle_approx_values_lock);
#endif //COMMONCPP
  setData(data,m_last_saved_clock_angle_approx_values[tellp()]);
} // void CSimulation::setLastTrajectoryClockAngleApproxData(const T& data)


//------------------------------------------------------------------------------
template<class T, class U>
inline void CSimulation::setData(const T& s, U& d)
{
  d.resize(s.size());
  for (unsigned int i = 0; i< d.size();++i) d[i] = s[i];
} // void setData(const T& s, U& d)


//------------------------------------------------------------------------------
template<class T>
inline void 
CSimulation::setLastTrajectoryClockAngleApproxDataDerivative(const T& data_d)
{
#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_clock_angle_approx_values_dot_lock);
#endif //COMMONCPP
  setData(data_d,m_last_saved_clock_angle_approx_values_dot[tellp()]);
} // CSimulation::setLastTrajectoryClockAngleApproxDataDerivative(const T&)


//------------------------------------------------------------------------------
template<class T, class U>
inline void CSimulation::setLastTrajectoryClockAngleApproxData(const T& d, const U& d_d)
{
  setLastTrajectoryClockAngleApproxData(d);
  setLastTrajectoryClockAngleApproxDataDerivative(d_d);
} // void setLastTrajectoryClockAngleApproxData(const T& d, const U& d_d)


//------------------------------------------------------------------------------
template<class T>
inline void CSimulation::setLastTrajectoryConeAngleData(const T& data)
{
#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_cone_angle_history_lock);
#endif //COMMONCPP
  setData(data, m_last_saved_cone_angle_history);
} // void setLastTrajectoryConeAngleData(const T& data)


//------------------------------------------------------------------------------
template<class T>
inline void 
CSimulation::setLastTrajectoryConeAngleApproxDataDerivative(const T& data_d)
{
#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_cone_angle_approx_values_dot_lock);
#endif //COMMONCPP
  setData(data_d,m_last_saved_cone_angle_approx_values_dot[tellp()]);
} // void CSimulation::setLastTrajectoryConeAngleApproxDataDerivative(const T& data_d)


//------------------------------------------------------------------------------
template<class T, class U>
inline void 
CSimulation::setLastTrajectoryConeAngleApproxData(const T& d, const U& d_d)
{
  setLastTrajectoryConeAngleApproxData(d);
  setLastTrajectoryConeAngleApproxDataDerivative(d_d);
} // void CSimulation::setLastTrajectoryConeAngleApproxData(const T& d, const U& d_d)


//------------------------------------------------------------------------------
template<class T>
inline void CSimulation::setLastTrajectoryConeAngleApproxData(const T& data)
{
#if defined COMMONCPP
  SimpleGuard guard(m_last_saved_cone_angle_approx_values_lock);
#endif //COMMONCPP
  setData(data, m_last_saved_cone_angle_approx_values[tellp()]);
} // inline void setLastTrajectoryConeAngleApproxData(const T& data)


//------------------------------------------------------------------------------
inline void CSimulation::saveClockApproxParams(const std::string& filename) const
{
  saveApproxParams(m_best_clock_knot_time,m_best_clock_knot_values,filename,
    m_best_clock_knot_time.size());
} // void CSimulation::saveClockApproxParams(const std::string& filename) const


//------------------------------------------------------------------------------
inline void CSimulation::setFlyByInitConditions(CIndividual& i, double time)
{
  setRendezvousInitConditions(i,time);
} // void CSimulation::setFlyByInitConditions(CIndividual& i, double time)


//------------------------------------------------------------------------------
inline void CSimulation::setRetrogradeImpactInitConditions(CIndividual& i
                                                           , double time)
{
  setRendezvousInitConditions(i,time);
} // void setRetrogradeImpactInitConditions(CIndividual& i, double time)


//------------------------------------------------------------------------------
inline void CSimulation::setActProblemInitConditions(CIndividual& i, double time)
{
  setRendezvousInitConditions(i,time);
  m_best_prox_sum = INIT_BEST_PROX_SUM;
} // void CSimulation::setActProblemInitConditions(CIndividual& i, double time)


//------------------------------------------------------------------------------
inline void CSimulation::setSpiProblemInitConditions(CIndividual& i, double t)
{
  setRendezvousInitConditions(i, t);
} // void CSimulation::setSpiProblemInitConditions(CIndividual& i, double t)


//------------------------------------------------------------------------------
inline void CSimulation::saveConeApproxParams(const std::string& filename) const
{
  saveApproxParams(m_best_cone_knot_time,m_best_cone_knot_values,filename
    ,m_best_cone_knot_time.size());
} // void CSimulation::saveConeApproxParams(const std::string& filename) const


//------------------------------------------------------------------------------
inline void CSimulation::incp(void)
{
  ++m_approx_put_index;
  if (m_approx_put_index>=m_last_saved_trajectory_approx.size())
    m_approx_put_index=0;
} // void CSimulation::incp(void)


//------------------------------------------------------------------------------
inline void CSimulation::decp(void)
{
  if (m_approx_put_index>=m_last_saved_trajectory_approx.size()
    || m_approx_put_index==0)
    m_approx_put_index=m_last_saved_trajectory_approx.size()-1;
  else --m_approx_put_index;
} // void CSimulation::decp(void)


//------------------------------------------------------------------------------
inline void CSimulation::incg(void)
{
  ++m_approx_get_index;
  if (m_approx_get_index>=m_last_saved_trajectory_approx.size())
    m_approx_get_index=0;
} // void CSimulation::incg(void)

//------------------------------------------------------------------------------
inline void CSimulation::decg(void)
{
  if (m_approx_get_index>=m_last_saved_trajectory_approx.size()
    || m_approx_get_index==0)
    m_approx_get_index=m_last_saved_trajectory_approx.size()-1;
  else --m_approx_get_index;
} // void CSimulation::decg(void)


//------------------------------------------------------------------------------
inline unsigned int CSimulation::tellp(void)
{
  return bounded(m_approx_put_index,uint(0),m_last_saved_trajectory_approx.size());
} // unsigned int CSimulation::tellp(void)


//------------------------------------------------------------------------------
inline unsigned int CSimulation::seekp(unsigned int i)
{
  return m_approx_put_index=i%m_last_saved_trajectory_approx.size();
} // unsigned int CSimulation::seekp(unsigned int i)


//------------------------------------------------------------------------------
inline unsigned int CSimulation::tellg(void)
{
  return bounded(m_approx_get_index,uint(0),m_last_saved_trajectory_approx.size());
} // unsigned int CSimulation::tellg(void)


//------------------------------------------------------------------------------
inline unsigned int CSimulation::seekg(unsigned int i) 
{
  return m_approx_get_index=i%m_last_saved_trajectory_approx.size();
} // unsigned int CSimulation::seekg(unsigned int i)


//------------------------------------------------------------------------------
inline void CSimulation::activateClockAngleHistory(void) 
{ 
  m_clock_angle_history = true;
}


//------------------------------------------------------------------------------
inline void CSimulation::deactivateClockAngleHistory(void) 
{
  m_clock_angle_history = false;
}


//------------------------------------------------------------------------------
inline void CSimulation::activateClockAngleApproximation(void) 
{
  m_clock_angle_history_approximation = true;
}


//------------------------------------------------------------------------------
inline void CSimulation::deactivateClockAngleApproximation(void) 
{
  m_clock_angle_history_approximation = false;
}


//------------------------------------------------------------------------------
inline void CSimulation::setClockAngleApproximationParameters(
  const int segSize, const int polyDegree)
{
  m_clock_angle_approx_seg_size    = segSize;
  m_clock_angle_approx_poly_degree = polyDegree;
}

//------------------------------------------------------------------------------
inline void CSimulation::activateConeAngleHistory(void) 
{
  m_cone_angle_history = true;
}


//------------------------------------------------------------------------------
inline void CSimulation::deactivateConeAngleHistory(void) 
{
  m_cone_angle_history = false;
}


//------------------------------------------------------------------------------
inline void CSimulation::setConeAngleApproximationParameters(const int segSize
                                                             , const int polyDegree)
{
  m_cone_angle_approx_seg_size    = segSize;
  m_cone_angle_approx_poly_degree = polyDegree;
}


//------------------------------------------------------------------------------
inline void CSimulation::activateConeAngleApproximation(void) 
{
  m_cone_angle_history_approximation = true;
}


//------------------------------------------------------------------------------
inline void CSimulation::deactivateConeAngleApproximation(void) 
{
  m_cone_angle_history_approximation = false;
}


//------------------------------------------------------------------------------
inline void CSimulation::activateThrustHistory(void) 
{ 
  m_thrust_history = true;
}


//------------------------------------------------------------------------------
inline void CSimulation::deactivateThrustHistory(void) 
{
  m_thrust_history = false;
}


//------------------------------------------------------------------------------
inline void CSimulation::getBasicClockAngleData(const SCTraj& t
                                                , CDataArray& a, const uint n_max_points)
{
  a.resize(n_max_points>0?(n_max_points<t.size()?n_max_points:t.size()):t.size());
  for (unsigned int i=0;i<a.size();++i) a[i]=t[i].getClockAngle();
} // void CSimulation::getBasicClockAngleData(...)


//------------------------------------------------------------------------------
template<class T>
inline void CSimulation::postprocessAngleData(T& a)
{
  double max_jump(0.0),t(0.0);
  do
  {
    max_jump = 0.0;
    uint i(0);
    for ( i = 1;i<a.size();++i)
    {
      if (max_jump<(t= std::max(max_jump, std::max(a[i ]-a[i-1], a[i-1]-a[i]))))
      {
        max_jump = t;
      }
      if (max_jump > m_max_clock_jump) break;
    }
    if (max_jump > m_max_clock_jump)
    {
      if (a[i]<a[i-1]) a[i]+= TWOPI;
      else a[i-1] += TWOPI;
    }
  } while (max_jump > m_max_clock_jump);
} // 


//------------------------------------------------------------------------------
inline void CSimulation::getBasicConeAngleData(
  const SCTraj& t, CDataArray& a, const uint n_max_points)
{
  a.resize(n_max_points>0?(n_max_points<t.size()?n_max_points:t.size()):t.size());
  for (unsigned int i=0;i<a.size();++i) a[i]=t[i].getConeAngle();
} // void CSimulation::getBasicConeAngleData(...)

//------------------------------------------------------------------------------
inline const CSC& CSimulation::getCurrentSC(void) const
{
  return m_spacecraft;
} // const CSC& CSimulation::getCurrentSpacecraft(void) const


//------------------------------------------------------------------------------
inline CSC& CSimulation::getCurrentSC(void)
{
  return m_spacecraft;
} // const CSC& CSimulation::getCurrentSpacecraft(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setIntegrator(MathLib::RK4* itgr)
{
  assert( itgr != 0 );
  m_SpacecraftOrbit_RK4 = itgr;
} // void CSimulation::setIntegrator(MathLib::RK4* itgr)


//------------------------------------------------------------------------------
inline MathLib::RK4* CSimulation::getIntegratorRK4(void) const
{
  return m_SpacecraftOrbit_RK4;
} // MathLib::RK4* CSimulation::getIntegratorRK4(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setIntegrator(MathLib::RKF54* itgr)
{
  assert( itgr != 0 );
  m_SpacecraftOrbit_RKF54 = itgr;
} // void CSimulation::setIntegrator(MathLib::RKF54* itgr)


//------------------------------------------------------------------------------
inline void CSimulation::setIntegrator(MathLib::DOPRI8* itgr)
{
  assert( itgr != 0 );
  m_SpacecraftOrbit_DOPRI8  = itgr;
} // void CSimulation::setIntegrator(MathLib::DOPRI8* itgr)


//------------------------------------------------------------------------------
inline void CSimulation::setIntegrator(MathLib::DE* itgr)
{
  assert( itgr != 0 );
  m_SpacecraftOrbit_DESG = itgr;
} // void CSimulation::setIntegrator(MathLib::DE* itgr)


//------------------------------------------------------------------------------
inline MathLib::RKF54* CSimulation::getIntegratorRKF54(void) const
{
  return m_SpacecraftOrbit_RKF54;
} // MathLib::RKF54* CSimulation::getIntegratorRKF54(void) const


//------------------------------------------------------------------------------
inline MathLib::DOPRI8* CSimulation::getIntegratorDOPRI8(void) const
{
  return m_SpacecraftOrbit_DOPRI8;
} // MathLib::DOPRI8* CSimulation::getIntegratorDOPRI8(void) const


//------------------------------------------------------------------------------
inline MathLib::DE* CSimulation::getIntegratorDESG(void) const
{
  return m_SpacecraftOrbit_DESG;
} // MathLib::DE* CSimulation::getIntegratorDESG(void) const


//------------------------------------------------------------------------------
inline void CSimulation::setIntegrationStepSize(const double step_size)
{
  m_h_dyn = step_size;
} // void CSimulation::setIntegrationStepSize(const double step_size)


//------------------------------------------------------------------------------
inline double CSimulation::getIntegrationStepSize(void)
{
  return m_h_dyn;
} // double CSimulation::getIntegrationStepSize(void)


//------------------------------------------------------------------------------
inline int CSimulation::getIntegrationVectorDimension(void)
{
  return m_itgr_dim;
} // int CSimulation::getIntegrationVectorDimension(void)


//------------------------------------------------------------------------------
inline void CSimulation::setIntegrationVectorDimension(const int new_dim)
{
  m_itgr_dim = new_dim;
} // void CSimulation::setIntegrationVectorDimension(const int dim)


//------------------------------------------------------------------------------
inline const std::valarray<double>& 
CSimulation::getTOEElementDeviation(void) const
{
  return m_delements;
}


//------------------------------------------------------------------------------
inline std::valarray<double>& CSimulation::getTOEElementDeviation(void)
{
  return m_delements;
}


//------------------------------------------------------------------------------
inline const std::valarray<double>& 
CSimulation::getTOEBestElementDeviation(void) const
{
  return m_best_delements;
}


//------------------------------------------------------------------------------
inline 
void CSimulation::setTOEBestElementDeviation(const std::valarray<double>& new_v)
{
#if defined development_version
  assert(m_best_delements.size() == new_v.size());
#endif
  m_best_delements = new_v;
}


//------------------------------------------------------------------------------
inline std::valarray<double>& CSimulation::getTOERemChangeTime(void)
{
  return m_element_change_rest_times;
}


//------------------------------------------------------------------------------
inline std::valarray<double>& CSimulation::getTOEBestRemChangeTime(void)
{
  return m_best_element_change_rest_times;
}


//------------------------------------------------------------------------------
inline
void CSimulation::setTOEBestRemChangeTime(const std::valarray<double>& new_v)
{
#if defined development_version
  assert(m_best_element_change_rest_times.size() == new_v.size());
#endif
  m_best_element_change_rest_times = new_v;
}

//------------------------------------------------------------------------------
inline std::valarray<double>& CSimulation::getTOEInitialElementDeviation(void)
{
  return m_initial_delements;
} // std::valarray<double>& getTOEInitialDeviation(void)

//------------------------------------------------------------------------------
inline double CSimulation::getCurrMaxDSSCOrbAngle(void) const
{
  return m_curr_DSSC_max_angle;
} // inline double getCurrMaxDSSCOrbAngle(void) const

//------------------------------------------------------------------------------
inline void CSimulation::setCurrMaxDSSCOrbAngle(const double angle)
{
  m_curr_DSSC_max_angle = angle;
} // inline void setCurrMaxDSSCOrbAngle(const double angle)

//------------------------------------------------------------------------------
inline bool CSimulation::isHeliocentric(void) const
{
  return m_data.isHeliocentric();
} // bool isHeliocentric(void) const

//------------------------------------------------------------------------------
inline bool CSimulation::isBodyCentric(void) const
{
  return m_data.isBodycentric();
} // bool isBodyCentric(void) const

//------------------------------------------------------------------------------
inline SpaceLib::BodyID CSimulation::getInitCentralBody(void) const
{
  return m_data.getInitCentralBody();
} // SpaceLib::BodyID getInitCentralBody(void) const

#endif // CSIMULATION_H_INCLUDED
