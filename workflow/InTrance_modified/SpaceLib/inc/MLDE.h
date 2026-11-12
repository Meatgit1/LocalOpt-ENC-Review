/*!
\file MLDE.h
\version 1.7
\brief Classes for the numerical integration ordinary DES (RK4, RKF45, DOPRI8, DE).
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\date 11.11.2006 (AO)

\b Purpose:
   This module provides implemenations of the 4th-order Runge-Kutta method
   and the variable order variable stepsize multistep method of Shampine &
   Gordon.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.
   This code is based on the code from:
   Montenbruck, Gill: Satellite Orbits. Springer, 2000

\b Copyright:
   (c) 1999-2004 Oliver Montenbruck, Eberhard Gill
*/

#ifndef DE_H
#define DE_H

#include <MLMatrix.h>
#include <MLVector.h>

namespace MathLib
{

  class Vector;

  /// Information what kind of integrators can be used.
  enum DESSolverType {
    _NOT_SPECIFIED,             ///< no special solver
    _RK4,                       ///< Runge-Kutta 4(th) order integrator
    _RKF54,                     ///< Runge-Kutta 4(5) order predictor
                                ///< integrator (stepsize control)
    _DOPRI8,                    ///< Dormand-Prince 7(8) order predictor
                                ///< integrator (step size control +
                                ///< abort function)
    _DESG                       ///< DES integrator after Shampine and Gordon
  };

  /// Status code return values of Runge-Kutta integration step methods.
  enum RK_STATE
  {
    RK_DONE = 1,                ///< Successful step
    RK_ERROR = 2,               ///< Test for local error could not be performed
    RK_BADACC = 3,              ///< Accuracy requirement could not be achieved
    RK_NUMSTEPS = 4,            ///< Permitted number of steps exceeded
    RK_STOP_FUNCTION = 5,        ///< Numerical integration process has been cancelled
      ///< by user provided stop function
    RK_POSTCHECK_ERROR          ///< check of state vector showed error
  };


//------------------------------------------------------------------------------
//
// RK4 class (specification)
//
//------------------------------------------------------------------------------

  /// Function prototype for first order differential equations
  /// void f (double x, const Vector& y, Vector& yp )

    /**
       This declaration of the DES - function satisfies the need of the user
       provided nature of the problem to integrate. The function that returns
       the first derivation of the input vector \em y with respect to the
       independent variable \em x and possible auxiliary data that is accessible
       through the pointer \em pAux is implemented in the calling program. With
       the passing of a reference to this integration class it can be used to
       integrate the DES numerically.
       \param x Independent integration variable.
       \param y vector at time \em t
       \param pAux a pointer to data that could be needed by the DES-function to
              work properly and therefor should be passed by every call of the
              DES-function (optional)

       \returns \em yp first derivation of vector \em y at the time \em t

       \sa StopFunction
     */
  typedef void (*DESfunct) (double x, const Vector & y, Vector & yp, void *pAux);

  /// Declaration of a stop function.

    /**
       Currently only the DOPRI8 class is able to use a stop function.
       After every internal integration step this function is called. Depending
       on the returned value the integration continues or stops and returns
       the value of the vector \em y computed until the time \em t the stop signal
       occured. Because of the nature of such a stop function it is only possible
       to use it with variable step size integrators such as DOPRI8.

       \param y    The result vector of the last integration step.
       \param pAux A pointer to auxiliary data that might be needed by the stop function
       \return \em true if the integration should be stopped and \em false otherwise
       \sa RKfunct
    */
  typedef bool(*StopFunction) (double t, Vector& y, void *pAux);


/////////////////////////////////////////////////////////////////////////////

  /// An integrator class using the classical Runge-Kutta scheme.

    /**
       This is the most simpel type of an Runge-Kutta integration scheme.
       It does no accuracy checks and no step size control. It just
       computes a solution of the function \em f at the time \em t+delta_t
       with an accuracy depending on the step size set by the parameter
       \em h. In order to use this integrator for the solving of your DES
       you've to create an instance of it and supply it with the DES function.
       The integration can be started by calling the \em Step - method with
       the initial vector \em y and a step size \em h. The integration can be
       continued by subsequent calls of \em Step until the final value of
       the independent variable \em t is reached.
       \sa RKF54, DOPRI8
     */


  class RK4 {

  private:

    int n_eqn;         // number of elements in DES vector
    DESfunct f;        // the function that incorporates the DES
    void* pAux;        // pointer to a data structure if needed by function f
    Vector k_1;        // the k-multiplicators of the RUNGE-KUTTA scheme
    Vector k_2;
    Vector k_3;
    Vector k_4;

  public:

    /// Constructor of integrator class.
    /**
     \param f_     Differential equation system function(DES).
     \param n_eqn_ Dimension of DES.
     \param pAux_   Pointer to auxiliary data possibly used by \em f_.
    */
    RK4(DESfunct f_, int n_eqn_, void *pAux_)
      : n_eqn(n_eqn_),
        f(f_),
        pAux(pAux_),
        k_1(n_eqn_),
        k_2(n_eqn_),
        k_3(n_eqn_),
        k_4(n_eqn_)  {};

    /// Elementary integration step.
    void Step(double &t,  ///< Value of the independent variable; updated by \em t+h
              Vector & y, ///< Value of \em y(t); updated by \em y(t+h)
              double h    ///< Step size (unchanged)
      );
  };





  /// Runge-Kutta-Fehlberg 4(5) predictor integrator with step size control.
  /**
     The class uses a Runge-Kutta-Fehlberg 4(5) integration scheme for the
     integration of differential equation system (DES). The integration
     is bound to a maximum number of calls to the function f that
     is in fact the DES. It is provided by the user of this class.
     \sa RK4, DOPRI8
   */
  class RKF54
  {
    int n_eqn;
    DESfunct f;
    StopFunction sf;        // the stop function for earlier stop of integration
    void *pAux;
    void *pAux_stop_func;                 ///< data object to give stop function
    Vector k_1, k_2, k_3, k_4, k_5, k_6;
    Vector y_4th_order, y_5th_order, r;
    long maxcalls;
    double rel,abs;

  public:

    /**
    * \brief Constructor of integrator class.
    *
    * \param f_     Differential equation system function(DES).
    * \param n_eqn_ Dimension of DES.
    * \param pAux_   Pointer to auxiliary data possibly used by \em f_.
    */
    RKF54(DESfunct f_, int n_eqn_, void *pAux_);


    /// Initialization of the integrator class.
    void init(long _maxcalls   ///< maximum number of calls of des-function
      ,double _rel              ///< relative accuracy requirement
      ,double _abs               ///< absolute accuracy requirement
      ,MathLib::StopFunction asf = 0
      );

    ///< set new stop function
    void setStopFunc(MathLib::StopFunction new_sf, void* sf_data);

    void resetStopFunc(void);       // reset stop function related member fields

    /**
       The method tries to integrate DES \em f from \em t to \em t+delta_t.
       After the integration step, a check whether the accuracy conditions
       are obeyed is done. Depending on this check the integration step
       will be repeated with a smaller step size (limited to a lower boundary
       step size depending on the machine accuracy) or the result returned
       as \em y. There's also a maximum number of allowed calls to the
       DES-function. When the counter reaches this maximum and more calls
       are needed in order to perform the integration the method stops
       with the corresponding return value.

       \param t Value of the independent variable; updated by \em t+h
       \param y Value of DES vector \em y(t); updated by \em y(t+h)
       \param h Step Size (unchanged)
       \return A status indicating the result of the last integration step.
     */
    RK_STATE Step(double &t, Vector & y, double h);
  }; // RKF54

  // implementation

  //----------------------------------------------------------------------------
  inline void RKF54::setStopFunc(MathLib::StopFunction new_sf, void* sf_data)
   throw()
  {
    sf = new_sf;
    pAux_stop_func = sf_data;
  } // void setStopFunc(void* sf, void* sf_data)

  //----------------------------------------------------------------------------
  inline void RKF54::resetStopFunc(void) throw()
  {
    sf = 0;
    pAux_stop_func = 0;
  } // void resetStopFunc(void) throw()



/////////////////////////////////////////////////////////////////////////////

  /// Dormand-Prince 7(8) predictor integrator with step size control and
  /// user provided stop function.
  /**
     The class uses a Runge-Kutta 7(8) integration scheme developed by
     Dormand & Prince for the integration of differential equation system
     (DES). The integration is bound to a maximum number of calls of
     function \em f that is in fact the DES. It is provided by the user
     of this class. The integration process can also be stopped by a
     stop function which can (but does not have) to be provided by the
     user of this class. This is simply a function that returns \em true
     or \em false depending on the computed DES-vector \em y.
     \sa RKF54, RK4
   */
  class DOPRI8 {

  private:

    int n_eqn;                  // the dimension of the differential equation system (DES)
    DESfunct f;                 // the differential equation to solve
    StopFunction sf;        // the stop function for earlier stop of integration
    void *pAux;                 // pointer to auxiliary data
    Vector k_1, k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9, k_10, k_11, k_12, k_13;
    Vector y_7th_order, y_8th_order, r;
    long maxcalls;              // maximum number of calls of function f
    double rel;                 // relative accuracy requirement
    double abs;                 // absolute accuracy requirement

  public:

    /**
    * \brief Constructor of integrator class.
    *
    *  \param f_     Differential equation system function(DES).
    *  \param n_eqn_ Dimension of DES.
    *  \param pAux_   Pointer to auxiliary data possibly used by \em f_.
    */
    DOPRI8(DESfunct f_, int n_eqn_, void *pAux_);

    /// Initialization of the integrator class.

      /**
         Does the initialization of the integrtor in order to perform
         the numerical integration according to the requirements of
         the user.
         \param _maxcalls maximum number of calls of DES-function
         \param _rel relativ accuracy requirements
         \param _abs absolute accuracy requirements
         \param a_sF a function that can stop the integration depending on the
                     value of \em y
       */
    void Init(long _maxcalls, double _rel, double _abs, StopFunction a_sF = 0);


    /// One integrating step.

      /**
         The method tries to integrate DES \em f from \em t to \em t+delta_t.
         After the integration step, a check whether the accuracy conditions
         are obeyed is done. Depending on this check the integration step
         will be repeated with a smaller step size (limited to a lower boundary
         step size depending on the machine accuracy) or the result returned
         as \em y. There's also a maximum number of allowed calls to the
         DES-function. When the counter reaches this maximum and more calls
         are needed in order to perform the integration the method stops
         with the corresponding return value.
         There's another criteria that can stop the integration - the user
         provided stop function (see \em Init - method). After every
         successfull integration step (that means the accuracy requirements
         have been archieved) this function called. The return value
         (\em true or \em false) decides, whether the integration should be
         continued or not.

         \param t Value of the independent variable; updated by \em t+h
         \param y Value of DES vector \em y(t); updated by \em y(t+h)
         \param h Step Size (unchanged)
         \return A status indicating the result of the last integration step.
      */
    RK_STATE Step(double &t, Vector & y, double h);
  };



/////////////////////////////////////////////////////////////////////////////

  /// State codes that will be store in \em State
  /**
     After calling the 'Integ'-method the status code of the integrqation
     will be stored in the data field 'State'. The possible values of this
     data field are defined in this enum.
   */
  enum DE_STATE {
    DE_INIT = 1,                /// Restart integration
    DE_DONE = 2,                /// Successful step
    DE_BADACC = 3,              /// Accuracy requirement could not be achieved
    DE_NUMSTEPS = 4,            /// Permitted number of steps exceeded
    DE_STIFF = 5,               /// Stiff problem suspected
    DE_INVPARAM = 6             /// Invalid input parameters
  };


  /// 4th order Runge-Kutta multistep solver; variable order/ variable step size

  /**
     This module provides implemenations of the 4th-order Runge-Kutta method
     and the variable order variable stepsize multistep method of Shampine &
     Gordon.

     \b reference:

     Shampine, Gordon: "Computer solution of Ordinary Differential Equations",
     Freeman and Comp., San Francisco (1975).
   */
  class DE {
  public:


    /// Desired relative accuracy of the solution
    double relerr;

    /// Desired absolute accuracy of the solution
    double abserr;

    /// State code (default = DE_INIT)
    DE_STATE State;

    /// Flag for integrating past tout (default = true)
    bool PermitTOUT;

    /// Value of independent variable
    double t;

    /// default constructor

      /**
         If this constructor is used the integration can only be done
         if needed informations are provided with the 'Define' - method.
         \sa Define
       */
      DE(void);

      /**
      * \brief Init constructor.
      *
         To use the DE class for solving a DES it is also needed to
         provide some initial conditions like the beginning of
         integration. This can be done with the 'Init'-method.

         \param f_ The differential equation system (DES) function.
         \param n_eqn_ Dimension of DES.
         \param pAux_ Pointer to auxiliary data.

         \sa Init
       */
      DE(DESfunct f_, int n_eqn_, void *pAux_);


    /// Definition of the DES and related information.

      /**
         \param f_ The differential equation system (DES) function.
         \param n_eqn_ Dimension of DES.
         \param pAux_ Pointer to auxiliary data.
       */
    void Define(DESfunct f_, int n_eqn_, void *pAux_);


    /// Initialize the integrator with starting conditions and
    /// accuracy bounds

      /**
         \param t0  Initial value of the independent variable
         \param rel Relative accuracy requirement
         \param abs Absolute accuracy requirement
       */
    void Init(double t0, double rel, double abs);


    /// Integration (skips over warnings, aborts in case of error)

      /**
         \param tout Desired output point.
         \param y    The solution vector
         \return The computed solution in vector y.
       */
    void Integ(double tout, Vector & y);


    /// Interpolation of soulution \em y to a desired output point \em tout

      /**
         \param tout Desired output point.
         \param y reference of solution vector y
       */
    void Intrp(double tout, Vector & y);

  private:

    // Elements
    int n_eqn;
    DESfunct f;
    void *pAux;
    Vector yy, wt, p, yp, ypout;
    Matrix phi;
    double alpha[13], beta[13], v[13], w[13], psi[13];
    double sig[14], g[14];
    double x, h, hold, told, delsgn;
    int ns, k, kold;
    bool OldPermit, phase1, start, nornd;
    bool init;

    // Elementary integration step
    void Step(double &x, Vector & y, double &eps, bool & crash);

    // Interpolation
    void Intrp(double xout, Vector & yout, Vector & ypout);

    // Integration (with full control of warnings and error status codes)
    void Integ_(double &t,      // Initial value of the independent variable
      double tout,              // Desired output point
      Vector & y                // Solution Vector
      );
  };


  //----------------------------------------------------------------------------
  inline RKF54::RKF54(DESfunct f_, int n_eqn_, void *pAux_)
    :n_eqn(n_eqn_), f(f_), pAux(pAux_), k_1(n_eqn_),
    k_2(n_eqn_), k_3(n_eqn_), k_4(n_eqn_), k_5(n_eqn_),
    k_6(n_eqn_), y_4th_order(n_eqn_), y_5th_order(n_eqn_), r(n_eqn_) 
  {
  } // RKF54::RKF54

  //----------------------------------------------------------------------------
  inline DOPRI8::DOPRI8(DESfunct f_, int n_eqn_, void *pAux_) 
    :n_eqn(n_eqn_), f(f_), pAux(pAux_), k_1(n_eqn_),
     k_2(n_eqn_), k_3(n_eqn_), k_4(n_eqn_), k_5(n_eqn_),
     k_6(n_eqn_), k_7(n_eqn_), k_8(n_eqn_), k_9(n_eqn_),
     k_10(n_eqn_), k_11(n_eqn_), k_12(n_eqn_), k_13(n_eqn_),
     y_7th_order(n_eqn_), y_8th_order(n_eqn_), r(n_eqn_)
  {
  } // DOPRI8::DOPRI8

}                               // of namespace MathLib

#endif // of DE_H

// end of file

