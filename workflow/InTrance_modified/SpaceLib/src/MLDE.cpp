/*!
\file MLDE.cpp
\version 1.7
\brief Classes for the numerical integration ordinary DES (RK4, RKF45, DOPRI8, DE).
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\author Bernd Dachwald <bernd.dachwald@gmx.de>   (BD)
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

//SpaceLib
#include <configSpaceLib.h>

// C standard library
#include <math.h>

// STL
#include <limits>
#include <algorithm>

//SpaceLib
#include <Global.h>
#include <MLGlobal.h>
#include <MLDE.h>

namespace MathLib
{
  using namespace std;
  namespace
  {
    // Constants

    const int MAXNUM = 500;     ///< Maximum number of steps to take

    const double UMACH = std::numeric_limits < double >::epsilon(); ///< tbd
    const double TWOU = 2.0 * UMACH;                                ///< tbd
    const double FOURU = 4.0 * UMACH;                               ///< tbd

    /** \brief Namespace for auxiliary functions. */
    namespace helper
    {
      /** \brief Returns absolute value of a with sign of b. */
      double Sign(double a, double b)
      {
        return (b >= 0.0) ? fabs(a) : -fabs(a);
      }
    } // end of namespace helper

  } // end of (unnamed) namespace


  //------------------------------------------------------------------------------
  //
  // RK4 class (implementation)
  //
  //------------------------------------------------------------------------------

  // Step
  void RK4::Step(double &t, Vector & y, double h)
  {

    // Elementary RK4 step

    f(t, y, k_1, pAux);
    f(t + h / 2.0, y + (h / 2.0) * k_1, k_2, pAux);
    f(t + h / 2.0, y + (h / 2.0) * k_2, k_3, pAux);
    f(t + h, y + h * k_3, k_4, pAux);

    y = y + (h / 6.0) * (k_1 + 2.0 * k_2 + 2.0 * k_3 + k_4);

    // Update independent variable

    t = t + h;

  };


  //------------------------------------------------------------------------------
  //
  // RKF54 class (implementation)
  //
  // currently without step size control
  //
  //----------------------------------------------------------------------------

  // Initialization
  void RKF54::init(long _maxcalls, double _rel, double _abs
    , MathLib::StopFunction asF)
  {

    maxcalls = _maxcalls;
    rel = _rel;
    abs = _abs;
    sf = asF;
  };


  //----------------------------------------------------------------------------
  RK_STATE RKF54::Step(double &t, Vector & y, double h)
  {

    const int n = y.getDimension();     // DES dimension
    const double tend = t + h;  // end time
    const double hmax = h;      // max step size

    double hfinal = 0.0;        // optimum step size for the last integration step
    long fcalls = 0;            // number of function calls

    bool repeat = false;
    bool end = false;

    double et, tr, quot;

    while (true) 
    {

      if (!repeat) 
      {
        h = std::min(h, hmax); // max step size is current interval

        if (tend - t <= 1.25 * h)
        {
          hfinal = h;           // optimum step size for the last integration step
          end = true;
          h = tend - t;         // real step size for the last integration step
        }

        f(t, y, k_1, pAux);
        ++fcalls;             // update call counter
      } // end if (!repeat)

      f(t + 1.0 / 4.0 * h, y + h * (1.0 / 4.0 * k_1), k_2, pAux);
      f(t + 3.0 / 8.0 * h, y + h * (3.0 / 32.0 * k_1 + 9.0 / 32.0 * k_2), k_3, pAux);
      f(t + 12.0 / 13.0 * h, y + h * (1932.0 / 2197.0 * k_1
        - 7200.0 / 2197.0 * k_2 + 7296.0 / 2197.0 * k_3), k_4, pAux);
      f(t + h, y + h * (439.0 / 216.0 * k_1 - 8.0 * k_2 + 3680.0 / 513.0 * k_3 - 845.0 / 4104.0 * k_4), k_5, pAux);
      f(t + 1.0 / 2.0 * h, y + h * (-8.0 / 27.0 * k_1
        + 2.0 * k_2 - 3544.0 / 2565.0 * k_3 + 1859.0 / 4104.0 * k_4 - 11.0 / 40.0 * k_5), k_6, pAux);

      fcalls += 5;

      y_4th_order = y + h * (25.0 / 216.0 * k_1 // 4th order
        + 1408.0 / 2565.0 * k_3 + 2197.0 / 4104.0 * k_4 - 1.0 / 5.0 * k_5);

      y_5th_order = y + h * (16.0 / 135.0 * k_1 // 5th order
        + 6656.0 / 12825.0 * k_3 + 28561.0 / 56430.0 * k_4 - 9.0 / 50.0 * k_5 + 2.0 / 55.0 * k_6);

      r = y_5th_order - y_4th_order;

      quot = 0.0;

      for (int i = 0; i < n; ++i) 
      {
        et = rel * fabs(y_5th_order(i)) + abs;

        if (et == 0.0) 
        {
          cerr << "WARNING: RKF54::Step(): test for local error could not be performed" << endl;
          return RK_ERROR;
        } // if (et == 0.0)

        tr = fabs(r(i)) / et;
        quot = std::max(quot, tr);
      } // for (int i = 0; i < n; ++i)

      if (quot <= 1.0) 
      {
        y = y_5th_order;
        t += h;
        if (sf != 0)
        {
          if (sf(t, y, pAux_stop_func) == true)    // evaluate the stop function
            return RK_STOP_FUNCTION;            // and return if it signals true
        } // if (sf != 0)

        if (end) 
        {
          h = hfinal;
          return RK_DONE;
        }
        quot = std::max(quot, 0.00065336); // increase step size not more than 5 times
      }

      quot = std::min(quot, 4096.0); // decrease step size not more than 10 times

      h = 0.8 * h / sqrt(sqrt(quot)); // change step size

      if (h <= 1e3 * UMACH) 
      {
        cerr << "WARNING: RKF54::Step(): accuracy requirement could not be achieved" << endl;
        return RK_BADACC;
      }

      if (fcalls > maxcalls) 
      {
        cerr << "WARNING: RKF54::Step(): permitted number of steps exceeded" << endl;
        return RK_NUMSTEPS;
      }

      if (quot > 1.0) 
      {
        repeat = true;
        end = false;
      }
      else repeat = false;
    }  // while (true)
  };


  //------------------------------------------------------------------------------
  //
  // DOPRI8 class (implementation)
  //
  // with step size control and user supplied stop function
  //
  //------------------------------------------------------------------------------

  void DOPRI8::Init(long _maxcalls, double _rel, double _abs, StopFunction asF)
  {
    maxcalls = _maxcalls;
    rel = _rel;
    abs = _abs;
    sf = asF;
  };


  RK_STATE DOPRI8::Step(double &t, Vector & y, double h)
  {

    const int n = y.getDimension();     // DES dimension
    const double tend = t + h;  // end time
    const double hmax = h;      // max step size

    double hfinal = 0.0;        // optimum step size for the last integration step
    long fcalls = 0;            // number of function calls

    bool repeat = false;        // flag if computation should be repeated with modified step size
    bool end = false;

    while (true) {
      if (!repeat) {            // repeat with new step size == desired accuracy not met yet

        h = std::min(h, hmax); // max step size is current interval

        if (tend - t <= 1.25 * h) {
          hfinal = h;           // optimum step size for the last integration step
          end = true;
          h = tend - t;         // real step size for the last integration step
        }

        f(t, y, k_1, pAux);     // calculation of k1
        fcalls++;

      }                         // end if (!repeat)

      // calculation of the other k-values with user provided differential equation
      f(t + 1.0 / 18.0 * h, y + h * (1.0 / 18.0 * k_1), k_2, pAux);
      f(t + 1.0 / 12.0 * h, y + h * (1.0 / 48.0 * k_1 + 1.0 / 16.0 * k_2), k_3, pAux);
      f(t + 1.0 / 8.0 * h, y + h * (1.0 / 32.0 * k_1 + 3.0 / 32.0 * k_3), k_4, pAux);
      f(t + 5.0 / 16.0 * h, y + h * (5.0 / 16.0 * k_1 - 75.0 / 64.0 * k_3 + 75.0 / 64.0 * k_4), k_5, pAux);
      f(t + 3.0 / 8.0 * h, y + h * (3.0 / 80.0 * k_1 + 3.0 / 16.0 * k_4 + 3.0 / 20.0 * k_5), k_6, pAux);
      f(t + 59.0 / 400.0 * h, y + h * (29443841.0 / 614563906.0 * k_1
        + 77736538.0 / 692538347.0 * k_4
        - 28693883.0 / 1125000000.0 * k_5 + 23124283.0 / 1800000000.0 * k_6), k_7, pAux);
      f(t + 93.0 / 200.0 * h,
        y + h * (16016141.0 / 946692911.0 * k_1 +
        61564180.0 / 158732637.0 * k_4 +
        22789713.0 / 633445777.0 * k_5 +
        545815736.0 / 2771057229.0 * k_6 - 180193667.0 / 1043307555.0 * k_7), k_8, pAux);
      f(t + 5490023248.0 / 9719169821.0 * h,
        y + h * (39632708.0 / 573591083.0 * k_1 -
        433636366.0 / 683701615.0 * k_4 -
        421739975.0 / 2616292301.0 * k_5 +
        100302831.0 / 723423059.0 * k_6 +
        790204164.0 / 839813087.0 * k_7 + 800635310.0 / 3783071287.0 * k_8), k_9, pAux);
      f(t + 13.0 / 20.0 * h,
        y + h * (246121993.0 / 1340847787.0 * k_1 -
        37695042795.0 / 15268766246.0 * k_4 -
        309121744.0 / 1061227803.0 * k_5 -
        12992083.0 / 490766935.0 * k_6 +
        6005943493.0 / 2108947869.0 * k_7 +
        393006217.0 / 1396673457.0 * k_8 + 123872331.0 / 1001029789.0 * k_9), k_10, pAux);
      f(t + 1201146811.0 / 1299019798.0 * h,
        y + h * (-1028468189.0 / 846180014.0 * k_1 +
        8478235783.0 / 508512852.0 * k_4 +
        1311729495.0 / 1432422823.0 * k_5 -
        10304129995.0 / 1701304382.0 * k_6 -
        48777925059.0 / 3047939560.0 * k_7 +
        15336726248.0 / 1032824649.0 * k_8 -
        45442868181.0 / 3398467696.0 * k_9 + 3065993473.0 / 597172653.0 * k_10), k_11, pAux);
      f(t + h,
        y + h * (185892177.0 / 718116043.0 * k_1 -
        3185094517.0 / 667107341.0 * k_4 -
        477755414.0 / 1098053517.0 * k_5 -
        703635378.0 / 230739211.0 * k_6 +
        5731566787.0 / 1027545527.0 * k_7 +
        5232866602.0 / 850066563.0 * k_8 -
        4093664535.0 / 808688257.0 * k_9 +
        3962137247.0 / 1805957418.0 * k_10 + 65686358.0 / 487910083.0 * k_11), k_12, pAux);
      f(t + h,
        y + h * (403863854.0 / 491063109.0 * k_1 -
        5068492393.0 / 434740067.0 * k_4 -
        411421997.0 / 543043805.0 * k_5 +
        652783627.0 / 914296604.0 * k_6 +
        11173962825.0 / 925320556.0 * k_7 -
        13158990841.0 / 6184727034.0 * k_8 +
        3936647629.0 / 1978049680.0 * k_9 -
        160528059.0 / 685178525.0 * k_10 + 248638103.0 / 1413531060.0 * k_11), k_13, pAux);

      fcalls += 12;             // diff. equ. has been evaluated 12 times

      //      Vector y_7th_order(n), y_8th_order(n), r(n);

      y_7th_order = y + h * (13451932.0 / 455176623.0 * k_1 // 7th order result vector
        - 808719846.0 / 976000145.0 * k_6
        + 1757004468.0 / 5645159321.0 * k_7
        + 656045339.0 / 265891186.0 * k_8
        - 3867574721.0 / 1518517206.0 * k_9
        + 465885868.0 / 322736535.0 * k_10 + 53011238.0 / 667516719.0 * k_11 + 2.0 / 45.0 * k_12);

      y_8th_order = y + h * (14005451.0 / 335480064.0 * k_1 // 8th order result vector
        - 59238493.0 / 1068277825.0 * k_6
        + 181606767.0 / 758867731.0 * k_7
        + 561292985.0 / 797845732.0 * k_8
        - 1041891430.0 / 1371343529.0 * k_9
        + 760417239.0 / 1151165299.0 * k_10
        + 118820643.0 / 751138087.0 * k_11 - 528747749.0 / 2220607170.0 * k_12 + 1.0 / 4.0 * k_13);

      r = y_8th_order - y_7th_order;  // residue vector

      double quot = 0.0;

      for (int i = 0; i < n; ++i) 
      {
        const double et = rel * fabs(y_8th_order(i)) + abs; // error approximation

        if (et == 0.0) {
          cerr << "WARNING: DOPRI8::Step(): test for local error could not be performed" << endl;
          return RK_ERROR;
        }

        const double tr = fabs(r(i)) / et;

        quot = std::max(quot, tr);
      }

      if (quot <= 1.0) {        // desired accuracy has been reached
        y = y_8th_order;        // y_8th_order will is taken as value of y(n+1)

        t += h;                 // step size is increased by h

        if (sf != 0) {          // evaluate the stop function

          if (sf(t, y, pAux) == true)
            return RK_STOP_FUNCTION;  // and return if it signals true

        }

        if (end) {
          h = hfinal;
          return RK_DONE;
        }

        quot =  std::max(quot, 0.00065336); // increase step size not more than 5 times
      }

      quot = std::min(quot, 4096.0); // decrease step size not more than 10 times

      h = 0.8 * h / sqrt(sqrt(quot)); // change step size

      if (h <= 1e3 * UMACH) {
        cerr << "WARNING: DOPRI8::Step(): accuracy requirement could not be achieved" << endl;
        return RK_BADACC;
      }

      if (fcalls > maxcalls) {
        cerr << "WARNING: DOPPRI8::Step(): permitted number of steps exceeded" << endl;
        return RK_NUMSTEPS;
      }

      if (quot > 1.0) {
        repeat = true;
        end = false;
      }
      else {
        repeat = false;
      }

    }                           // of while(true) - loop
  };                            // of method MathLibDOPRI8::Step


  //------------------------------------------------------------------------------
  //
  // DE class (implementation)
  //
  //------------------------------------------------------------------------------

  /// Standard constructor.
  DE::DE(void):n_eqn(0), f(0), pAux(0)
  {
    State = DE_INVPARAM;        // Status flag
    PermitTOUT = false;         // Allow integration past tout by default
    // I HAVE CHANGED THAT TO FALSE
    t = 0.0;
    relerr = 0.0;               // Accuracy requirements
    abserr = 0.0;
  };



  DE::DE(DESfunct f_, int n_eqn_, void *pAux_):n_eqn(n_eqn_), f(f_), pAux(pAux_)
  {
    yy = Vector(n_eqn);         // Allocate Vectors with proper dimension
    wt = Vector(n_eqn);
    p = Vector(n_eqn);
    yp = Vector(n_eqn);
    ypout = Vector(n_eqn);
    phi = Matrix(n_eqn, 17);
    State = DE_INVPARAM;        // Status flag
    PermitTOUT = false;         // Allow integration past tout by default
    // I HAVE CHANGED THAT TO FALSE
    t = 0.0;
    relerr = 0.0;               // Accuracy requirements
    abserr = 0.0;
  };



  void DE::Define(DESfunct f_, int n_eqn_, void *pAux_)
  {
    n_eqn = n_eqn_;
    f = f_;
    pAux = pAux_;
    yy = Vector(n_eqn);         // Allocate Vectors with proper dimension
    wt = Vector(n_eqn);
    p = Vector(n_eqn);
    yp = Vector(n_eqn);
    ypout = Vector(n_eqn);
    phi = Matrix(n_eqn, 17);
    State = DE_INVPARAM;        // Status flag
    PermitTOUT = true;          // Allow integration past tout by default
    t = 0.0;
    relerr = 0.0;               // Accuracy requirements
    abserr = 0.0;
  };


  //
  // Integration step
  //

  void DE::Step(double &x, Vector & y, double &eps, bool & crash)
  {

    // Constants

    // Powers of two (two(n)=2**n)
    static const double two[14] = { 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0,
      256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0
    };

    static const double gstr[14] = { 1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188,
      0.0143, 0.0114, 0.00936, 0.00789, 0.00679,
      0.00592, 0.00524, 0.00468
    };


    // Variables

    bool success;
    int i, ifail, im1, ip1, iq, j, km1, km2, knew, kp1, kp2;
    int l, limit1, limit2, nsm2, nsp1, nsp2;
    double absh, erk, erkm1, erkm2, erkp1, err, hnew;
    double p5eps, r, reali, realns, rho, round, sum, tau;
    double temp1, temp2, temp3, temp4, temp5, temp6, xold;


    //
    // Begin block 0
    //
    // Check if step size or error tolerance is too small for machine
    // precision.  If first step, initialize phi array and estimate a
    // starting step size. If step size is too small, determine an
    // acceptable one.
    //

    if (fabs(h) < FOURU * fabs(x)) {
                h =  helper::Sign(FOURU * fabs(x), h);
      crash = true;
      return;                   // Exit
    };

    p5eps = 0.5 * eps;
    crash = false;
    g[1] = 1.0;
    g[2] = 0.5;
    sig[1] = 1.0;

    ifail = 0;

    // If error tolerance is too small, increase it to an
    // acceptable value.

    round = 0.0;
    for (l = 0; l < n_eqn; ++l)
      round += (y(l) * y(l)) / (wt(l) * wt(l));
    round = TWOU * sqrt(round);
    if (p5eps < round) {
      eps = 2.0 * round * (1.0 + FOURU);
      crash = true;
      return;
    };


    if (start) {
      // Initialize. Compute appropriate step size for first step.
      f(x, y, yp, pAux);
      sum = 0.0;
      for (l = 0; l < n_eqn; ++l)
      {
        phi(l, 1) = yp(l);
        phi(l, 2) = 0.0;
        sum += (yp(l) * yp(l)) / (wt(l) * wt(l));
      }
      sum = sqrt(sum);
      absh = fabs(h);
      if (eps < 16.0 * sum * h * h)
        absh = 0.25 * sqrt(eps / sum);
	  h = helper::Sign(std::max(absh, FOURU * fabs(x)), h);
      hold = 0.0;
      hnew = 0.0;
      k = 1;
      kold = 0;
      start = false;
      phase1 = true;
      nornd = true;
      if (p5eps <= 100.0 * round) 
      {
        nornd = false;
        for (l = 0; l < n_eqn; ++l)
          phi(l, 15) = 0.0;
      };
    };

    //
    // End block 0
    //


    //
    // Repeat blocks 1, 2 (and 3) until step is successful
    //
    do {

      //
      // Begin block 1
      //
      // Compute coefficients of formulas for this step. Avoid computing
      // those quantities not changed when step size is not changed.
      //

      kp1 = k + 1;
      kp2 = k + 2;
      km1 = k - 1;
      km2 = k - 2;

      // ns is the number of steps taken with size h, including the
      // current one. When k<ns, no coefficients change.

      if (h != hold)
        ns = 0;
      if (ns <= kold)
        ns = ns + 1;
      nsp1 = ns + 1;

      if (k >= ns) {

        // Compute those components of alpha[*],beta[*],psi[*],sig[*]
        // which are changed
        beta[ns] = 1.0;
        realns = ns;
        alpha[ns] = 1.0 / realns;
        temp1 = h * realns;
        sig[nsp1] = 1.0;
        if (k >= nsp1) {
          for (i = nsp1; i <= k; ++i)
          {
            im1 = i - 1;
            temp2 = psi[im1];
            psi[im1] = temp1;
            beta[i] = beta[im1] * psi[im1] / temp2;
            temp1 = temp2 + h;
            alpha[i] = h / temp1;
            reali = i;
            sig[i + 1] = reali * alpha[i] * sig[i];
          };
        };
        psi[k] = temp1;

        // Compute coefficients g[*]; initialize v[*] and set w[*].
        if (ns > 1) {
          // If order was raised, update diagonal part of v[*]
          if (k > kold) {
            temp4 = k * kp1;
            v[k] = 1.0 / temp4;
            nsm2 = ns - 2;
            for (j = 1; j <= nsm2; ++j) 
            {
              i = k - j;
              v[i] = v[i] - alpha[j + 1] * v[i + 1];
            };
          };
          // Update V[*] and set W[*]
          limit1 = kp1 - ns;
          temp5 = alpha[ns];
          for (iq = 1; iq <= limit1; ++iq) 
          {
            v[iq] = v[iq] - temp5 * v[iq + 1];
            w[iq] = v[iq];
          };

          g[nsp1] = w[1];
        }
        else {
          for (iq = 1; iq <= k; ++iq) 
          {
            temp3 = iq * (iq + 1);
            v[iq] = 1.0 / temp3;
            w[iq] = v[iq];
          };
        };

        // Compute the g[*] in the work Vector w[*]
        nsp2 = ns + 2;
        if (kp1 >= nsp2) {
          for (i = nsp2; i <= kp1; ++i) 
          {
            limit2 = kp2 - i;
            temp6 = alpha[i - 1];
            for (iq = 1; iq <= limit2; ++iq)
              w[iq] = w[iq] - temp6 * w[iq + 1];
            g[i] = w[1];
          };
        };

      };                        // if K>=NS

//
// End block 1
//


      //
      // Begin block 2
      //
      // Predict a solution p[*], evaluate derivatives using predicted
      // solution, estimate local error at order k and errors at orders
      // k, k-1, k-2 as if constant step size were used.
      //

      // Change phi to phi star
      if (k >= nsp1) {
        for (i = nsp1; i <= k; ++i)
        {
          temp1 = beta[i];
          for (l = 0; l < n_eqn; ++l)
            phi(l, i) = temp1 * phi(l, i);
        }
      };

      // Predict solution and differences
      for (l = 0; l < n_eqn; ++l) 
      {
        phi(l, kp2) = phi(l, kp1);
        phi(l, kp1) = 0.0;
        p(l) = 0.0;
      };
      for (j = 1; j <= k; ++j) 
      {
        i = kp1 - j;
        ip1 = i + 1;
        temp2 = g[i];
        for (l = 0; l < n_eqn; ++l) 
        {
          p(l) = p(l) + temp2 * phi(l, i);
          phi(l, i) = phi(l, i) + phi(l, ip1);
        };
      };
      if (nornd) {
        p = y + h * p;
      }
      else {
        for (l = 0; l < n_eqn; ++l) 
        {
          tau = h * p(l) - phi(l, 15);
          p(l) = y(l) + tau;
          phi(l, 16) = (p(l) - y(l)) - tau;
        };
      };
      xold = x;
      x = x + h;
      absh = fabs(h);
      f(x, p, yp, pAux);

      // Estimate errors at orders k, k-1, k-2
      erkm2 = 0.0;
      erkm1 = 0.0;
      erk = 0.0;

      for (l = 0; l < n_eqn; ++l) 
      {
        temp3 = 1.0 / wt(l);
        temp4 = yp(l) - phi(l, 1);
        if (km2 > 0)
          erkm2 = erkm2 + ((phi(l, km1) + temp4) * temp3)
            * ((phi(l, km1) + temp4) * temp3);
        if (km2 >= 0)
          erkm1 = erkm1 + ((phi(l, k) + temp4) * temp3)
            * ((phi(l, k) + temp4) * temp3);
        erk = erk + (temp4 * temp3) * (temp4 * temp3);
      };

      if (km2 > 0)
        erkm2 = absh * sig[km1] * gstr[km2] * sqrt(erkm2);
      if (km2 >= 0)
        erkm1 = absh * sig[k] * gstr[km1] * sqrt(erkm1);

      temp5 = absh * sqrt(erk);
      err = temp5 * (g[k] - g[kp1]);
      erk = temp5 * sig[kp1] * gstr[k];
      knew = k;

      // Test if order should be lowered
      if (km2 > 0)
		  if (std::max(erkm1, erkm2) <= erk)
          knew = km1;
      if (km2 == 0)
        if (erkm1 <= 0.5 * erk)
          knew = km1;

      //
      // End block 2
      //


      //
      // If step is successful continue with block 4, otherwise repeat
      // blocks 1 and 2 after executing block 3
      //

      success = (err <= eps);

      if (!success) {

        //
        // Begin block 3
        //

        // The step is unsuccessful. Restore x, phi[*,*], psi[*]. If
        // 3rd consecutive failure, set order to 1. If step fails more
        // than 3 times, consider an optimal step size. Double error
        // tolerance and return if estimated step size is too small
        // for machine precision.
        //

        // Restore x, phi[*,*] and psi[*]
        phase1 = false;
        x = xold;
        for (i = 1; i <= k; ++i) 
        {
          temp1 = 1.0 / beta[i];
          ip1 = i + 1;
          for (l = 0; l < n_eqn; ++l)
            phi(l, i) = temp1 * (phi(l, i) - phi(l, ip1));
        };

        if (k >= 2)
          for (i = 2; i <= k; ++i)
            psi[i - 1] = psi[i] - h;

        // On third failure, set order to one.
        // Thereafter, use optimal step size
        ifail++;
        temp2 = 0.5;
        if (ifail > 3)
          if (p5eps < 0.25 * erk)
            temp2 = sqrt(p5eps / erk);
        if (ifail >= 3)
          knew = 1;
        h = temp2 * h;
        k = knew;
        if (fabs(h) < FOURU * fabs(x)) {
          crash = true;
          h = helper::Sign(FOURU * fabs(x), h);
          eps *= 2.0;
          return;               // Exit
        };

        //
        // End block 3, return to start of block 1
        //

      };                        // end if(success)

    }
    while (!success);


    //
    // Begin block 4
    //
    // The step is successful. Correct the predicted solution, evaluate
    // the derivatives using the corrected solution and update the
    // differences. Determine best order and step size for next step.
    //

    kold = k;
    hold = h;


    // Correct and evaluate
    temp1 = h * g[kp1];
    if (nornd)
      for (l = 0; l < n_eqn; ++l)
        y(l) = p(l) + temp1 * (yp(l) - phi(l, 1));
    else
      for (l = 0; l < n_eqn; ++l) {
        rho = temp1 * (yp(l) - phi(l, 1)) - phi(l, 16);
        y(l) = p(l) + rho;
        phi(l, 15) = (y(l) - p(l)) - rho;
      };

    f(x, y, yp, pAux);


    // Update differences for next step
    for (l = 0; l < n_eqn; ++l) 
    {
      phi(l, kp1) = yp(l) - phi(l, 1);
      phi(l, kp2) = phi(l, kp1) - phi(l, kp2);
    };
    for (i = 1; i <= k; ++i)
      for (l = 0; l < n_eqn; ++l)
        phi(l, i) = phi(l, i) + phi(l, kp1);


    // Estimate error at order k+1 unless
    // - in first phase when always raise order,
    // - already decided to lower order,
    // - step size not constant so estimate unreliable
    erkp1 = 0.0;
    if ((knew == km1) || (k == 12))
      phase1 = false;

    if (phase1) {
      k = kp1;
      erk = erkp1;
    }
    else {
      if (knew == km1) {
        // lower order
        k = km1;
        erk = erkm1;
      }
      else {

        if (kp1 <= ns) {
          for (l = 0; l < n_eqn; l++)
            erkp1 = erkp1 + (phi(l, kp2) / wt(l)) * (phi(l, kp2) / wt(l));
          erkp1 = absh * gstr[kp1] * sqrt(erkp1);

          // Using estimated error at order k+1, determine
          // appropriate order for next step
          if (k > 1) {
			  if (erkm1 <= std::min(erk, erkp1)) {
              // lower order
              k = km1;
              erk = erkm1;
            }
            else {
              if ((erkp1 < erk) && (k != 12)) {
                // raise order
                k = kp1;
                erk = erkp1;
              };
            };
          }
          else {
            if (erkp1 < 0.5 * erk) {
              // raise order
              // Here erkp1 < erk < max(erkm1,ermk2) else
              // order would have been lowered in block 2.
              // Thus order is to be raised
              k = kp1;
              erk = erkp1;
            };
          };

        };                      // end if kp1<=ns

      };                        // end if knew!=km1

    };                          // end if !phase1


    // With new order determine appropriate step size for next step
    if (phase1 || (p5eps >= erk * two[k + 1]))
      hnew = 2.0 * h;
    else {
      if (p5eps < erk) {
        temp2 = k + 1;
        r = pow(p5eps / erk, 1.0 / temp2);
		hnew = absh * std::max(0.5, std::min(0.9, r));
		hnew = helper::Sign(std::max(hnew, FOURU * fabs(x)), h);
      }
      else
        hnew = h;
    };

    h = hnew;

    //
    // End block 4
    //
  };


  //
  // Interpolation
  //

  void DE::Intrp(double xout, Vector & yout, Vector & ypout)
  {

    // Variables

    int i, j, ki;
    double eta, gamma, hi, psijm1;
    double temp1, term;
    double g[14], rho[14], w[14];


    g[1] = 1.0;
    rho[1] = 1.0;

    hi = xout - x;
    ki = kold + 1;

    // Initialize w[*] for computing g[*]
    for (i = 1; i <= ki; ++i) 
    {
      temp1 = i;
      w[i] = 1.0 / temp1;
    }

    // Compute g[*]
    term = 0.0;
    for (j = 2; j <= ki; ++j)
    {
      psijm1 = psi[j - 1];
      gamma = (hi + term) / psijm1;
      eta = hi / psijm1;
      for (i = 1; i <= ki + 1 - j; ++i)
      {
        w[i] = gamma * w[i] - eta * w[i + 1];
      }
      g[j] = w[1];
      rho[j] = gamma * rho[j - 1];
      term = psijm1;
    };

    // Interpolate for the solution yout and for
    // the derivative of the solution ypout
    ypout = 0.0;
    yout = 0.0;
    for (j = 1; j <= ki; ++j) 
    {
      i = ki + 1 - j;
      yout = yout + g[i] * phi.Col(i);
      ypout = ypout + rho[i] * phi.Col(i);
    };
    yout = yy + hi * yout;
  };


  //
  // DE integration
  // (with full control of warnings and errros status codes)
  //

  void DE::Integ_(double &t,    // Value of independent variable
    double tout,                // Desired output point
    Vector & y                  // Solution Vector
    )
  {

    // Variables

    bool stiff, crash;          // Flags
    int nostep;                 // Step count
    int kle4 = 0;
    double releps, abseps, tend;
    double absdel, del, eps;


    // Return, if output time equals input time

    if (t == tout)
      return;                   // No integration


    // Test for improper parameters

	eps = std::max(relerr, abserr);

    if ((relerr < 0.0) ||       // Negative relative error bound
      (abserr < 0.0) ||         // Negative absolute error bound
      (eps <= 0.0) ||           // Both error bounds are non-positive
      (State > DE_INVPARAM) ||  // Invalid status flag
      ((State != DE_INIT) && (t != told))) {
      State = DE_INVPARAM;      // Set error code
      return;                   // Exit
    };


    // On each call set interval of integration and counter for
    // number of steps. Adjust input error tolerances to define
    // weight Vector for subroutine STEP.

    del = tout - t;
    absdel = fabs(del);

    tend = t + 100.0 * del;
    if (!PermitTOUT)
      tend = tout;

    nostep = 0;
    kle4 = 0;
    stiff = false;
    releps = relerr / eps;
    abseps = abserr / eps;

    if ((State == DE_INIT) || (!OldPermit) || (delsgn * del <= 0.0)) {
      // On start and restart also set the work variables x and yy(*),
      // store the direction of integration and initialize the step size
      start = true;
      x = t;
      yy = y;
      delsgn = helper::Sign(1.0, del);
	  h = helper::Sign(std::max(FOURU * fabs(x), fabs(tout - x)), tout - x);
    }

    while (true) {              // Start step loop

      // If already past output point, interpolate solution and return
      if (fabs(x - t) >= absdel) {
        Intrp(tout, y, ypout);
        State = DE_DONE;        // Set return code
        t = tout;               // Set independent variable
        told = t;               // Store independent variable
        OldPermit = PermitTOUT;
        return;                 // Normal exit
      };

      // If cannot go past output point and sufficiently close,
      // extrapolate and return
      if (!PermitTOUT && (fabs(tout - x) < FOURU * fabs(x))) {
        h = tout - x;
        f(x, yy, yp, pAux);     // Compute derivative yp(x)
        y = yy + h * yp;        // Extrapolate Vector from x to tout
        State = DE_DONE;        // Set return code
        t = tout;               // Set independent variable
        told = t;               // Store independent variable
        OldPermit = PermitTOUT;
        return;                 // Normal exit
      };

      // Test for too much work
      if (nostep >= MAXNUM) {
        State = DE_NUMSTEPS;    // Too many steps
        if (stiff)
          State = DE_STIFF;     // Stiffness suspected
        y = yy;                 // Copy last step
        t = x;
        told = t;
        OldPermit = true;
        return;                 // Weak failure exit
      };

      // Limit step size, set weight Vector and take a step
	  h = helper::Sign(std::min(fabs(h), fabs(tend - x)), h);
      for (int l = 0; l < n_eqn; l++)
        wt(l) = releps * fabs(yy(l)) + abseps;

      Step(x, yy, eps, crash);


      // Test for too small tolerances
      if (crash) {
        State = DE_BADACC;
        relerr = eps * releps;  // Modify relative and absolute
        abserr = eps * abseps;  // accuracy requirements
        y = yy;                 // Copy last step
        t = x;
        told = t;
        OldPermit = true;
        return;                 // Weak failure exit
      }

      nostep++;                 // Count total number of steps

      // Count number of consecutive steps taken with the order of
      // the method being less or equal to four and test for stiffness
      kle4++;
      if (kold > 4)
        kle4 = 0;
      if (kle4 >= 50)
        stiff = true;

    }                           // End step loop


  };


  //
  // Initialization
  //

  void DE::Init(double t0,      // Initial value of the independent variable
    double rel,                 // Relative accuracy requirement
    double abs                  // Absolute accuracy requirement
    )
  {
    t = t0;
    relerr = rel;
    abserr = abs;
    State = DE_INIT;
  };


  //
  // DE integration with simplified state code handling
  // (skips over warnings, aborts in case of error)
  //

  void DE::Integ(double tout,   // Desired output point
    Vector & y                  // Solution Vector
    )
  {
    do {
      Integ_(t, tout, y);
      if (State == DE_INVPARAM) {
        std::cerr << "ERROR: invalid parameters in DE::Integ" << std::endl;
        exit(1);
      }
      if (State == DE_BADACC) {
        std::cerr << "WARNING: Accuracy requirement not achieved in DE::Integ" << std::endl;
      }
      if (State == DE_STIFF) {
        std::cerr << "WARNING: Stiff problem suspected in DE::Integ" << std::endl;
      }
    }
    while (State > DE_DONE);
  };


  // Interpolation
  //----------------------------------------------------------------------------
  void DE::Intrp(double tout,   // Desired output point
    Vector & y                  // Solution Vector
    )
  {
    Intrp(tout, y, ypout); // Interpolate and discard interpolated deriv. ypout
  }

}                               // of namespace MathLib

// end of file

