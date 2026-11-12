/*!
\file Euler.cpp
\version 1.7
\brief A class to simplify the handling of rotation sequences/transformation matrices (Euler).
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 06.01.2005 (AO)

\b Purpose:
   Class Euler is intended to make handling of Euler angle rotation
   sequences/transformations easier. It provides functions to
   calculate transformation matrices and also to get the desired Euler angles
   back from a given transformation matrix. These operations are supported
   for every one of the 12 possible rotation sequences (three subsequent
   rotations about three different axis result in 12 possible rotation sequences).

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.
   
\b Copyright:
   (c) 2005 Andreas Ohndorf
*/

//SpaceLib
#include <configSpaceLib.h>

// standard C library
#include <math.h>

//SpaceLib
#include <Euler.h>
#include <MLConst.h>

using namespace MathLib;

namespace SpaceLib {

  // This structure stores information to get the right Attitude matrix
  // element. It is used in an array to store needed element information
  // (seven elements) for all 12 possible rotation sequences.
  // Those information is basically the row and column index in the matrix
  // and a factor to multiply with the matrix element's content. The factor
  // is whether 1.0, which won't change the element, or -1.0 to change the
  // element's sign.
  /** Tbd. */
  struct EulerLocation {
    double factor; ///< tbd
    int row;       ///< tbd 
    int col;       ///< tbd

    // Constructor to set the data fields. It's needed in the "const array"
    // below.
    /** Tbd. */
    inline EulerLocation(double _f, int _r, int _c):factor(_f), row(_r), col(_c) {};
  };



  /** This array of constant data structures of EulerLocation type holds
      information about the needed key elements of the Attitude matrix to
      determine correct values for the Euler angles. Correct means, that an
      Attitude matrix calculated from the resulting angles with the rotation
      sequence given before to determine these angles, is equal to the origin
      Attitude matrix. So the transformation is the same. */
  const struct EulerLocation ELocations[EULER_N_SEQUENCES][EULER_USED_ELEMENTS]={
    {  // Euler rotation sequence 1 (1-2-3)
      EulerLocation(1.0,2,0),EulerLocation(-1.0,2,1),EulerLocation(1.0,2,2),
      EulerLocation(-1.0,1,0),EulerLocation(1.0,0,0),EulerLocation(1.0,0,1),
      EulerLocation(1.0,1,2)
    },
    {  // Euler rotation sequence 2 (1-3-2)
      EulerLocation(-1.0,1,0),EulerLocation(1.0,1,2),EulerLocation(1.0,1,1),
      EulerLocation(1.0,2,0),EulerLocation(1.0,0,0),EulerLocation(1.0,0,1),
      EulerLocation(1.0,2,2)
    },
    { // Euler rotation sequence 3 (2-3-1)
      EulerLocation(1.0,0,1),EulerLocation(-1.0,0,2),EulerLocation(1.0,0,0),
      EulerLocation(-1.0,2,1),EulerLocation(1.0,1,1),EulerLocation(1.0,1,2),
      EulerLocation(1.0,2,0)
    },
    { // Euler rotation sequence 4 (2-1-3)
      EulerLocation(-1.0,2,1),EulerLocation(1.0,2,0),EulerLocation(1.0,2,2),
      EulerLocation(1.0,0,1),EulerLocation(1.0,1,1),EulerLocation(1.0,1,2),
      EulerLocation(1.0,0,0)
    },
    { // Euler rotation sequence 5 (3-1-2)
      EulerLocation(1.0,1,2),EulerLocation(-1.0,1,0),EulerLocation(1.0,1,1),
      EulerLocation(-1.0,0,2),EulerLocation(1.0,2,2),EulerLocation(1.0,2,0),
      EulerLocation(1.0,0,1)
    },
    { // Euler rotation sequence 6 (3-2-1)
      EulerLocation(-1.0,0,2),EulerLocation(1.0,0,1),EulerLocation(1.0,0,0),
      EulerLocation(1.0,1,2),EulerLocation(1.0,2,2),EulerLocation(1.0,2,0),
      EulerLocation(1.0,1,1)
    },
    { // Euler rotation sequence 7 (1-2-1)
      EulerLocation(1.0,0,0),EulerLocation(1.0,0,1),EulerLocation(-1.0,0,2),
      EulerLocation(1.0,1,0),EulerLocation(1.0,2,0),EulerLocation(1.0,1,1),
      EulerLocation(1.0,1,2)
    },
    { // Euler rotation sequence 8 (1-3-1)
      EulerLocation(1.0,0,0),EulerLocation(1.0,0,2),EulerLocation(1.0,0,1),
      EulerLocation(1.0,2,0),EulerLocation(-1.0,1,0),EulerLocation(1.0,2,2),
      EulerLocation(-1.0,2,1)
    },
    { // Euler rotation sequence 9 (2-1-2)
      EulerLocation(1.0,1,1),EulerLocation(1.0,1,0),EulerLocation(1.0,1,2),
      EulerLocation(1.0,0,1),EulerLocation(-1.0,2,1),EulerLocation(1.0,0,0),
      EulerLocation(-1.0,0,2)
    },
    { // Euler rotation sequence 10 (2-3-2)
      EulerLocation(1.0,1,1),EulerLocation(1.0,1,2),EulerLocation(-1.0,1,0),
      EulerLocation(1.0,2,1),EulerLocation(1.0,0,1),EulerLocation(1.0,2,2),
      EulerLocation(1.0,2,0),
    },
    { // Euler rotation sequence 11 (3-1-3)
      EulerLocation(1.0,2,2),EulerLocation(1.0,2,0),EulerLocation(-1.0,2,1),
      EulerLocation(1.0,0,2),EulerLocation(1.0,1,2),EulerLocation(1.0,0,0),
      EulerLocation(1.0,0,1),
    },
    { // Euler rotation sequence 12 (3-2-3)
      EulerLocation(1.0,2,2),EulerLocation(1.0,2,1),EulerLocation(1.0,2,0),
      EulerLocation(1.0,1,2),EulerLocation(-1.0,0,2),EulerLocation(1.0,1,1),
      EulerLocation(-1.0,1,0)
    }
  };


  //---------------------------------------------------------------------------
  Matrix Euler::getEulerTransformation(
    const EulerSequences seq,
    const double a1,
    const double a2,
    const double a3
  ) {

    const double EPS = 1E-10;

    double s1 = sin(a1);
    double c1 = cos(a1);
    double s2 = sin(a2);
    double c2 = cos(a2);
    double s3 = sin(a3);
    double c3 = cos(a3);

    Matrix T(3,3);

    switch(seq) {

      case E123:
        T(0,0) = c3*c2;   T(0,1) = c3*s2*s1+s3*c1;   T(0,2) = -c3*s2*c1+s3*s1;
        T(1,0) = -s3*c2;  T(1,1) = -s3*s2*s1+c3*c1;  T(1,2) = s3*s2*c1+c3*s1;
        T(2,0) = s2;      T(2,1) = -c2*s1;           T(2,2) = c2*c1;

      break;

      case E132:
        T(0,0) = c3*c2;  T(0,1) = c3*s2*c1+s3*s1;  T(0,2) = c3*s2*s1-s3*c1;
        T(1,0) = -s2;    T(1,1) = c2*c1;           T(1,2) = c2*s1;
        T(2,0) = s3*c2;  T(2,1) = s3*s2*c1-c3*s1;  T(2,2) = s3*s2*s1+c3*c1;
      break;

      case E231:
        T(0,0) = c2*c1;            T(0,1) = s2;      T(0,2) = -c2*s1;
        T(1,0) = -c3*s2*c1+s3*s1;  T(1,1) = c3*c2;   T(1,2) = c3*s2*s1+s3*c1;
        T(2,0) = s3*s2*c1+c3*s1;   T(2,1) = -s3*c2;  T(2,2) = -s3*s2*s1+c3*c1;
      break;

      case E213:
        T(0,0) = c3*c1+s3*s2*s1;   T(0,1) = s3*c2;  T(0,2) = -c3*s1+s3*s2*c1;
        T(1,0) = -s3*c1+c3*s2*s1;  T(1,1) = c3*c2;  T(1,2) = s3*s1+c3*s2*c1;
        T(2,0) = c2*s1;            T(2,1) = -s2;    T(2,2) = c2*c1;
      break;

      case E312:
        T(0,0) = c3*c1-s3*s2*s1;  T(0,1) = c3*s1+s3*s2*c1;  T(0,2) = -s3*c2;
        T(1,0) = -c2*s1;          T(1,1) = c2*c1;           T(1,2) = s2;
        T(2,0) = s3*c1+c3*s2*s1;  T(2,1) = s3*s1-c3*s2*c1;  T(2,2) = c3*c2;
      break;

      case E321:
        T(0,0) = c2*c1;            T(0,1) = c2*s1;           T(0,2) = -s2;
        T(1,0) = -c3*s1+s3*s2*c1;  T(1,1) = c3*c1+s3*s2*s1;  T(1,2) = s3*c2;
        T(2,0) = s3*s1+c3*s2*c1;   T(2,1) = -s3*c1+c3*s2*s1; T(2,2) = c3*c2;
      break;

      case E121:
        T(0,0) = c2;    T(0,1) = s2*s1;           T(0,2) = -s2*c1;
        T(1,0) = s3*s2; T(1,1) = c3*c1-s3*c2*s1;  T(1,2) = c3*s1+s3*c2*c1;
        T(2,0) = c3*s2; T(2,1) = -s3*c1-c3*c2*s1; T(2,2) = -s3*s1+c3*c2*c1;

        // The following lines are needed because due to rounding/truncating
        // errors, it's possible to get values below -1.0 or greater 1.0 for
        // certain elements of the Attitude matrix. This happens if angle 2 is
        // a multiple of Pi (0�/180�). Here, in case of E121, it is needed
        // to treat the elements T(1,1), T(1,2), T(2,1) and T(2,2) especially.
        // The same procedure has to be applied for all kind of Type 2 rotation
        // sequences.

        if ((fabs(T(1,0)) < EPS) && (fabs(T(2,0)) < EPS)) {
          if (c2>0.0) {
            T(1,1) = cos(a1+a3); T(1,2) = sin(a1+a3);
            T(2,1) = -T(1,2);    T(2,2) = T(1,1);
          } else {
            T(1,1) = cos(a1-a3); T(1,2) = sin(a1-a3);
            T(2,1) = T(1,2);     T(2,2) = -T(1,1);
          }
        }
      break;

      case E131:
        T(0,0) = c2;     T(0,1) = s2*c1;           T(0,2) = s2*s1;
        T(1,0) = -c3*s2; T(1,1) = c3*c2*c1-s3*s1;  T(1,2) = c3*c2*s1+s3*c1;
        T(2,0) = s3*s2;  T(2,1) = -s3*c2*c1-c3*s1; T(2,2) = -s3*c2*s1+c3*c1;

        if ((fabs(T(1,0)) < EPS) && (fabs(T(2,0)) < EPS)) {
          if (c2>0.0) {
            T(1,1) = cos(a1+a3); T(1,2) = sin(a1+a3);
            T(2,1) = -T(1,2);    T(2,2) = T(1,1);
          } else {
            T(1,1) = -cos(a1-a3); T(1,2) = -sin(a1-a3);
            T(2,1) = T(1,2);     T(2,2) = -T(1,1);
          }
        }
      break;

      case E212:
        T(0,0) = c3*c1-s3*c2*s1; T(0,1) = s3*s2;  T(0,2) = -c3*s1-s3*c2*c1;
        T(1,0) = s2*s1;          T(1,1) = c2;     T(1,2) = s2*c1;
        T(2,0) = s3*c1+c3*c2*s1; T(2,1) = -c3*s2; T(2,2) = -s3*s1+c3*c2*c1;

        if ((fabs(T(0,1)) < EPS) && (fabs(T(2,1)) < EPS)) {
          if (c2>0.0) {
            T(0,0) = cos(a1+a3); T(0,2) = -sin(a1+a3);
            T(2,0) = -T(0,2);    T(2,2) = T(0,0);
          } else {
            T(0,0) = cos(a1-a3); T(0,2) = -sin(a1-a3);
            T(2,0) = T(0,2);     T(2,2) = -T(0,0);
          }
        }
      break;

      case E232:
        T(0,0) = c3*c2*c1-s3*s1; T(0,1) = c3*s2; T(0,2) = -c3*c2*s1-s3*c1;
        T(1,0) = -s2*c1;         T(1,1) = c2;    T(1,2) = s2*s1;
        T(2,0) = s3*c2*c1+c3*s1; T(2,1) = s3*s2; T(2,2) = -s3*c2*s1+c3*c1;

        if ((fabs(T(0,1)) < EPS) && (fabs(T(2,1)) < EPS)) {
          if (c2>0.0) {
            T(0,0) = cos(a1+a3); T(0,2) = -sin(a1+a3);
            T(2,0) = -T(0,2);    T(2,2) = T(0,0);
          } else {
            T(0,0) = -cos(a1-a3); T(0,2) = sin(a1-a3);
            T(2,0) = T(0,2);     T(2,2) = -T(0,0);
          }
        }
      break;

      case E313:
        T(0,0) = c3*c1-s3*c2*s1;  T(0,1) = c3*s1+s3*c2*c1;  T(0,2) = s3*s2;
        T(1,0) = -s3*c1-c3*c2*s1; T(1,1) = -s3*s1+c3*c2*c1; T(1,2) = c3*s2;
        T(2,0) = s2*s1;           T(2,1) = -s2*c1;          T(2,2) = c2;

        if ((fabs(T(1,2)) < EPS) && (fabs(T(0,2)) < EPS)) {
          if (c2>0.0) {
            T(0,0) = cos(a1+a3); T(0,1) = sin(a1+a3);
            T(1,0) = -T(0,1);    T(1,1) = T(0,0);
          } else {
            T(0,0) = cos(a1-a3); T(0,1) = sin(a1-a3);
            T(1,0) = T(0,1);     T(1,1) = -T(0,0);
          }
        }

      break;

      case E323:
        T(0,0) = c3*c2*c1-s3*s1;  T(0,1) = c3*c2*s1+s3*c1;  T(0,2) = -c3*s2;
        T(1,0) = -s3*c2*c1-c3*s1; T(1,1) = -s3*c2*s1+c3*c1; T(1,2) = s3*s2;
        T(2,0) = s2*c1;           T(2,1) = s2*s1;           T(2,2) = c2;

        if ((fabs(T(1,2)) < EPS) && (fabs(T(0,2)) < EPS)) {
          if (c2>0.0) {
            T(0,0) = cos(a1+a3); T(0,1) = sin(a1+a3);
            T(1,0) = -T(0,1);    T(1,1) = T(0,0);
          } else {
            T(0,0) = -cos(a1-a3); T(0,1) = -sin(a1-a3);
            T(1,0) = T(0,1);      T(1,1) = -T(0,0);
          }
        }
      break;

    } // end of switch(seq)

    return T;

  }


  //---------------------------------------------------------------------------
  void Euler::getEulerAnglesType1(
    const Matrix& T,
    double& a1,
    double& a2,
    double& a3,
    const enum EulerSequences seq
  ) {

    const double EPS = 1E-6;         // the used numerical accuracy
    double Psi[2], Theta[2], Phi[2]; // temporary arrays to store possible angles

    // sin(a2)
    double vE1 = ELocations[seq][E1].factor *
      T(ELocations[seq][E1].row, ELocations[seq][E1].col);

    // cos(a2)sin(a1)
    double vE2 = ELocations[seq][E2].factor *
      T(ELocations[seq][E2].row, ELocations[seq][E2].col);

    // cos(a2)cos(a1)
    double vE3 = ELocations[seq][E3].factor *
      T(ELocations[seq][E3].row, ELocations[seq][E3].col);

    // cos(a2)sin(a3)
    double vE4 = ELocations[seq][E4].factor *
      T(ELocations[seq][E4].row, ELocations[seq][E4].col);

    // cos(a2)cos(a3)
    double vE5 = ELocations[seq][E5].factor *
      T(ELocations[seq][E5].row, ELocations[seq][E5].col);

    // E123,E231,E312: sin(a3)cos(a1)+cos(a3)sin(a2)sin(a1)
    // E321,E213,E132: sin(a3)sin(a1)+cos(a3)sin(a2)cos(a1)
//    double vE6 = ELocations[seq][E6].factor *
//      T(ELocations[seq][E6].row, ELocations[seq][E6].col);

    // E123,E231,E312: cos(a3)sin(a1)+sin(a3)sin(a2)cos(a1)
    // E321,E213,E132: cos(a3)cos(a1)+sin(a3)sin(a2)sin(a1)
//    double vE7 = ELocations[seq][E7].factor *
//      T(ELocations[seq][E7].row, ELocations[seq][E7].col);

    //from 231
    Theta[0] = asin(vE1);
    if (Theta[0] < 0.0) {
      Theta[1] = TWOPI + Theta[0];
      Theta[0] = PI - Theta[0];
    } else {
      Theta[1] = PI - Theta[0];
    }

    Psi[0] = atan2(vE2, vE3);
    if (Psi[0] < 0.0) {
      Psi[1] = PI + Psi[0];
      Psi[0] = TWOPI + Psi[0];
    } else {
      Psi[1] = PI - Psi[0];
    }

    Phi[0] = atan2(vE4, vE5);
    if (Phi[0] < 0.0) {
      Phi[1] = PI + Phi[0];
      Phi[0] = TWOPI + Phi[0];
    } else {
      Phi[1] = PI - Phi[0];
    }

    // The problem now is to estimate the right values of the three angles
    // The easiest way, but at the same time the slowest, is to calculate
    // the transformation matrix of the 313-Euler transformation and to
    // compare it with the given one. If they're equal, the angles have been
    // found. This way is not choosen due to performance needs.
    // Thus, the way it is done consists of three loops that test four
    // certain transformation matrix elements to get the right angles.

    int i_Psi, i_Theta, i_Phi;
    double s1,c1,c2,s3,c3;
    bool found = false;

    for (i_Psi = 0; i_Psi < 2; ++i_Psi) {
      s1 = sin(Psi[i_Psi]);
      c1 = cos(Psi[i_Psi]);

      for (i_Theta = 0; i_Theta < 2; ++i_Theta) {
        c2 = cos(Theta[i_Theta]);

        if ((fabs(c2*c1-vE3) < EPS) && (fabs(s1*c2-vE2) < EPS)){

          for (i_Phi = 0; i_Phi < 2; ++i_Phi) {

            c3 = cos(Phi[i_Phi]);
            s3 = sin(Phi[i_Phi]);

            if ((fabs(c3*c2-vE5) < EPS) && (fabs(s3*c2-vE4) < EPS)) {

              a1   = Psi[i_Psi];
              a2 = Theta[i_Theta];
              a3   = Phi[i_Phi];
              found = true;

            }
            if (found) break;

          } // Phi - loop

        }

        if (found) break;
      } // Theta - loop

      if (found) break;

    } // Psi - loop

  }


  //---------------------------------------------------------------------------
  void Euler::getEulerAnglesType2(
    const Matrix& T,
    double& a1,
    double& a2,
    double& a3,
    const enum EulerSequences seq
  ) {

    if (seq > E321) {

      const double EPS = 1E-6;         // the used numerical accuracy
      double Psi[2], Theta[2], Phi[2]; // temporary arrays to store possible angles
      int n_Psi, n_Theta, n_Phi;       // indices for the arrays

      // cos(a2)
      double vE1 = ELocations[seq][E1].factor *
        T(ELocations[seq][E1].row, ELocations[seq][E1].col);

      // sin(a2)sin(a1)
      double vE2 = ELocations[seq][E2].factor *
        T(ELocations[seq][E2].row, ELocations[seq][E2].col);

      // sin(a2)cos(a1)
      double vE3 = ELocations[seq][E3].factor *
        T(ELocations[seq][E3].row, ELocations[seq][E3].col);

      // sin(a2)sin(a3)
      double vE4 = ELocations[seq][E4].factor *
        T(ELocations[seq][E4].row, ELocations[seq][E4].col);

      // sin(a2)cos(a3)
      double vE5 = ELocations[seq][E5].factor *
        T(ELocations[seq][E5].row, ELocations[seq][E5].col);

      // cos(a3)cos(a1)-sin(a3)cos(a2)sin(a1)
      double vE6 = ELocations[seq][E6].factor *
        T(ELocations[seq][E6].row, ELocations[seq][E6].col);

      // cos(a3)sin(a1)+sin(a3)cos(a2)cos(a1)
      double vE7 = ELocations[seq][E7].factor *
        T(ELocations[seq][E7].row, ELocations[seq][E7].col);

      Theta[0] = acos(vE1); //cos(a2)

      // is the resulting angle greater than zero?
      if (Theta[0] > EPS) {
        // Angle Theta is not zero; This will be the case most of the time.

        // is Theta = 180 degree or expressed differently: is T33 = -1.0?
        if (fabs(vE1+1.0)<EPS) {   // Theta is 180 degree

          n_Theta = 1;  // there is only this solution

          // T11 now expresses the cosinus of the difference angle (Psi-Phi)

          if (fabs(vE6-1.0)<EPS) { // difference between Psi and Phi is zero

          // Because the difference angle is zero, an indefinite number
          // of angles of equal angles of Psi and Phi produce the same
          // result. So the trivial solution of Psi = Phi = 0.0 will do also.
          // This is the way it's done here.

            Psi[0] = 0.0;
            n_Psi  = 1;

            Phi[0] = 0.0;
            n_Phi  = 1;

          } else {

          // T11 now is a real angle difference (thus not zero). In this
          // constellation angle Psi will be set to zero and the resulting
          // angle difference assigned to Phi.

            Psi[0] = 0.0;
            n_Psi  = 1;

            Phi[0] = acos(vE6);
            Phi[1] = TWOPI-Phi[0];
            n_Phi = 2;

          }

        } else {

          // Angle Theta is not 180 degree

          Theta[1] = TWOPI - Theta[0];
          n_Theta = 2;

          Psi[0] = -HALFPI-atan2(vE3, vE2);
          while (Psi[0] < 0.0) Psi[0]+=PI;
          Psi[1] = Psi[0] + PI;
          n_Psi = 2;

          Phi[0] = -HALFPI-atan2(vE5, vE4);
          while (Phi[0] < 0.0) Phi[0]+=PI;
          Phi[1] = Phi[0] + PI;
          n_Phi = 2;
        }

      } else{

        // Angle Theta is zero. In this case angle Psi is not defined, because
        // the angles Psi and Phi are in the same plane of rotation. Thus
        // we are only able to compute the sum of Psi and Theta. Here, we set
        // Psi to zero and assign the resulting sum to Phi.

        Theta[0] = 0.0;           // Theta is set to zero
        n_Theta = 1;

        Psi[0] = 0.0;             // Psi is set to zero
        n_Psi = 1;

        Phi[0] = acos(vE6);    // and the sum of Psi and Phi stored
        Phi[1] = TWOPI-Phi[0];   // in Phi
        n_Phi = 2;
      }

      // The problem now is to estimate the right values of the three angles
      // The easiest way, but at the same time the slowest, is to calculate
      // the transformation matrix of the 313-Euler transformation and to
      // compare it with the given one. If they're equal, the angles have been
      // found. This way is has not been choosen doue to performance needs.
      // Thus, the way it is done is somehow complicated and consists of
      // certain checks following one after another.

      int i_Psi, i_Theta, i_Phi;
      double s1,c1,s2,s3,c3;
      bool found = false;

      // First, we check whether Theta is 0 or 180 degree
      // if it is, T33 (which is cos(Theta)) raised to the
      // power of two should be 1.
      if (fabs(vE1*vE1 - 1.0) < EPS) {

        // Yes, it is 0 or 180 degree.

        a1 = Psi[0];                // So save Psi and Theta.
        a2 = Theta[0];

        if (fabs(vE1+1.0) < EPS) {

          // Theta is 180 degree

          if ((fabs(cos(Phi[0])-vE6)<EPS) &&
              (fabs(-sin(Phi[0])-vE7)<EPS)) a3=Phi[0];
          else a3=Phi[1];

        } else {

          // Theta is 0 degree

          if ((fabs(cos(Phi[0])-vE6)<EPS) &&
              (fabs(sin(Phi[0])-vE7)<EPS)) a3=Phi[0];
          else a3=Phi[1];

        }

      } else {

        // No, Theta is neither 0 nor 180 degree.
        // So we have to probe the different combinations of angles and stop
        // when we've found the right ones. Instead of comparing the
        // whole matix, we just need to check the Elements T31, T32, T13 and T23.

        for (i_Psi = 0; i_Psi < n_Psi; ++i_Psi) {
          s1 = sin(Psi[i_Psi]);
          c1 = cos(Psi[i_Psi]);

          for (i_Theta = 0; i_Theta < n_Theta; ++i_Theta) {
            s2 = sin(Theta[i_Theta]);

            if ((fabs(s2*s1-vE2) < EPS) && (fabs(s2*c1-vE3) < EPS)){

              for (i_Phi = 0; i_Phi < n_Phi; ++i_Phi) {
                c3 = cos(Phi[i_Phi]);
                s3 = sin(Phi[i_Phi]);

                if ((fabs(s3*s2-vE4) < EPS) && (fabs(c3*s2-vE5) < EPS)) {

                  a1   = Psi[i_Psi];
                  a2 = Theta[i_Theta];
                  a3   = Phi[i_Phi];
                  found = true;

                }
                if (found) break;

              } // Phi - loop

            } // T31,T32 check

            if (found) break;

          } // Theta - loop

          if (found) break;

        } // Psi - loop

      } // end of "Theta is not 0 or 180 degree" - branch

    } else a1 = a2 = a3 = 0.0;

  }


  //---------------------------------------------------------------------------
  void Euler::getEulerAngles(
    const Matrix& T,
    double& a1,
    double& a2,
    double& a3,
    const enum EulerSequences seq
  ) {

    if (seq < E121) getEulerAnglesType1(T, a1, a2, a3, seq);
    else getEulerAnglesType2(T, a1, a2, a3, seq);

  }

} // of namespace SpaceLib
