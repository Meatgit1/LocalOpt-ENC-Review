/*!
\file DE405.cpp
\version 1.7
\brief A class to access the JPL ephemeris data files DE405 (DE405).
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\date 06.11.2004 (AO)

\b Purpose:
   The provided class is the C++ interface to the data file of JPL
   that holds the ephemeris data for the nine major planets, the
   sun and the moon. Currently it ranges from dates of 1600AC until
   2200AC.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.
   This code is based on the code by:
          David Hoffman/EG5
          NASA, Johnson Space Center
          Houston, TX 77058
          e-mail: david.a.hoffman1@jsc.nasa.gov

\b Copyright:
   (c) 2005 Andreas Ohndorf, Bernd Dachwald
*/

//SpaceLib
#include <configSpaceLib.h>

//STL
#include <ios>
#include <iomanip>
#include <sstream>

// standard C library
#include <math.h>

//SpaceLib
#include <DE405.h>

using namespace std;
using namespace MathLib;
using namespace SpaceLib;

namespace SpaceLib {


  // standard constructor
  //---------------------------------------------------------------------------
  DE405::DE405(void) : initialized(false) {
  }


  // data file constructor
  //---------------------------------------------------------------------------
  DE405::DE405(const std::string data_file_name)
   :_data_file_name(data_file_name), initialized(false) {
  }


  // destructor
  //---------------------------------------------------------------------------
  DE405::~DE405(void) {
    if (initialized) _data_file.close();
  }


  // set the data file name
  //---------------------------------------------------------------------------
  void DE405::setDataFile(const std::string data_file_name) {
    _data_file_name = data_file_name;
  }


  // initializes the object with the data file's content
  //---------------------------------------------------------------------------
  bool DE405::initialize(void) throw(FileErrorEXC) {

    if (initialized) {    // if true, it'll be a change of data file
      _data_file.close();  // close the file (maybe it is open)
      _data_file.clear();  // clear status flags to reset file stream
    }

    // now try to open data file
    if (_data_file_name.length() == 0) throw FileErrorEXC("Could not open ephemeris data file (no file name)!\n");
    _data_file.open(_data_file_name.c_str(), ios::binary | ios::in);

    if (!_data_file.is_open()) throw FileErrorEXC("Could not open ephemeris data file!\n");

    _data_file.read((char *)&H1, ARRAY_SIZE * sizeof(double));
    _data_file.read((char *)&H2, ARRAY_SIZE * sizeof(double));
    _data_file.read((char *)&Coeff_Array, ARRAY_SIZE * sizeof(double));

    // Store header data in global variables
    R1 = H1.data;
    R2 = H2.data;

    // Set current time variables
    T_beg = Coeff_Array[0];
    T_end = Coeff_Array[1];
    T_span = T_end - T_beg;

    // Debug Print (optional)
//    if (Debug) {
//    cout << "\n  In: Initialize_Ephemeris \n";
//    cout << "\n      ARRAY_SIZE = " << dec << setw(4) << ARRAY_SIZE;
//    cout << "\n      headerID   = " << dec << setw(3) << headerID;
//    cout << "\n      T_Beg      = " << setw(11) << setprecision(3) << T_beg;
//    cout << "\n      T_End      = " << setw(11) << setprecision(3) << T_end;
//    cout << "\n      T_Span     = " << setw(11) << setprecision(3) << T_span;
//    }

    stringstream err_msg;

    if ((int)R1.DENUM == EPHEMERIS) {
      return initialized = true;
    }
    else {
      err_msg << "\n Opened wrong file: " << _data_file_name
              << " for ephemeris: " << dec << EPHEMERIS << endl;
      throw FileErrorEXC(err_msg.str());
    }
  }


  // gets the position/velocity information out of the data file
  //---------------------------------------------------------------------------
  MathLib::Vector DE405::getState(const DE405StateIndex body,
                                  const SpaceLib::CTime& time)
   throw(DE405::IndexOutOfRangeEXC,
         DE405::TimeNotCoveredByFileEXC,
         DE405::NotInitializedEXC,
         FileErrorEXC) {

    if (!initialized)
     throw NotInitializedEXC("Object not initialized with data file!");

    double Time = time.JD();

    if (Time < H1.data.timeData[0] || Time > H1.data.timeData[1])
      throw TimeNotCoveredByFileEXC();

    double  A[50], B[50], Cp[50], P_Sum[3], V_Sum[3]
           ,Up[50], T_break, T_seg=0, T_sub=0, Tc=0;
    int i, j;
    long int C, G, N, offset = 0;
    stateType X;

    // This function doesn't "do" nutations or librations.
    // Don't call this member function to compute those data.
    if (body >= 11)
      throw IndexOutOfRangeEXC("Function called without proper body index!\n");

    // Initialize local coefficient array.
    for (i = 0; i < 50; i++) {
      A[i] = 0.0;
      B[i] = 0.0;
    }

    // Determine if a new record needs to be input
    // and read it, if needed.
    // Attention: This function may generate a FileErrorEXC - exception!
    if (Time < T_beg || Time > T_end) getCoefficients(Time);

    C = R1.coeffPtr[body][0] - 1; // Coeff array entry point
    N = R1.coeffPtr[body][1];     // Number of coeff's
    G = R1.coeffPtr[body][2];     // Granules in current record

    // debug print
/*
    if (Debug) {
      cout << "\n  In: Interpolate_State\n";
      cout << "\n  Target = " << setw(2) << Target;
      cout << "\n  C      = " << setw(4) << C << " (before)";
      cout << "\n  N      = " << setw(4) << N;
      cout << "\n  G      = " << setw(4) << G << endl;
    }
*/

    //  Compute the normalized time, then load the Tchebeyshev coefficients
    //  into array A[]. If T_span is covered by a single granule this is easy.
    //  If not, the granule that contains the interpolation time is found, and
    //  an offset from the array entry point for the ephemeris body is used to
    //  load the coefficients.

    if (G == 1) {
      Tc = 2.0 * (Time - T_beg) / T_span - 1.0;
      for (i = C; i < (C + 3 * N); i++) A[i - C] = Coeff_Array[i];
    } else if (G > 1) {
      T_sub = T_span / ((double)G); // Compute subgranule interval

      for (j = G; j > 0; j--) {
        T_break = T_beg + ((double)j - 1) * T_sub;
        // the following greater/equal comparison had to be changed from
        // just greater comparison, because there's been misbehaviour if
        // the given time matched interpolation interval start exactly
        if (Time >= T_break) {
          T_seg = T_break;
          offset = j - 1;
          break;
        }
      }

      Tc = 2.0 * (Time - T_seg) / T_sub - 1.0;
      C = C + 3 * offset * N;

      for (i = C; i < (C + 3 * N); i++) A[i - C] = Coeff_Array[i];
    }  else {            // Something has gone terribly wrong
      cout << "\n Number of granules must be >= 1: check header data.\n";
    }

    // Debug print (optional)
/*
    if (Debug) {
      cout << "\n  C      = " << setw(4) << dec << C << " (after)";
      cout << "\n  offset = " << setw(4) << dec << offset;
      cout << "\n  Time   = " << setw(12) << setprecision(7) << Time;
      cout << "\n  T_sub  = " << setw(12) << setprecision(7) << T_sub;
      cout << "\n  T_seg  = " << setw(12) << setprecision(7) << T_seg;
      cout << "\n  Tc     = " << setw(12) << setprecision(7) << Tc;
      cout << "\n  Array Coefficients:\n";
      for (i = 0; i < 3 * N; i++) {
        cout << "\n  A[" << setw(2) << i << "] = " << setw(22)
             << setprecision(15) << A[i];
      }
      cout << endl << endl;
    }
*/
    // Compute the interpolated position & velocity
    for (i = 0; i < 3; i++) {     // Compute interpolating polynomials
      Cp[0] = 1.0;
      Cp[1] = Tc;
      Cp[2] = 2.0 * Tc * Tc - 1.0;

      Up[0] = 0.0;
      Up[1] = 1.0;
      Up[2] = 4.0 * Tc;

      for (j = 3; j < N; j++) {
        Cp[j] = 2.0 * Tc * Cp[j - 1] - Cp[j - 2];
        Up[j] = 2.0 * Tc * Up[j - 1] + 2.0 * Cp[j - 1] - Up[j - 2];
      }

      P_Sum[i] = 0.0;             // Compute interpolated position & velocity
      V_Sum[i] = 0.0;

     for (j = N - 1; j > -1; j--) P_Sum[i] = P_Sum[i] + A[j + i * N] * Cp[j];
     for (j = N - 1; j > 0; j--)  V_Sum[i] = V_Sum[i] + A[j + i * N] * Up[j];

     X.Position[i] = P_Sum[i];
     X.Velocity[i] = V_Sum[i] * 2.0 * ((double)G) / (T_span * 86400.0);
    }

    return Vector( X.Position[0], X.Position[1], X.Position[2]
                  ,X.Velocity[0], X.Velocity[1], X.Velocity[2]);
  }


  // This function computes an array of libration angles from Chebys hev
  // coefficients read in from an ephemeris data file. The coefficients
  // are read in from a file by the function getCoefficients (when
  // necessary).
  //---------------------------------------------------------------------------
  MathLib::Vector DE405::getNutation(const DE405StateIndex body,
                                     const SpaceLib::CTime& time)
   throw(DE405::IndexOutOfRangeEXC,
         DE405::TimeNotCoveredByFileEXC,
         DE405::NotInitializedEXC,
         FileErrorEXC) {

    double A[50], Cp[50], sum[3], T_break, T_seg=0, T_sub=0, Tc=0;
    double Nutation[3];
    int i, j;
    long int C, G, N, offset = 0;

    // This function only computes nutations.
    if (body != DE405_SS_BARY)
     throw IndexOutOfRangeEXC("Nutation only for solar system barycenter!");

    if (!initialized)
     throw NotInitializedEXC("Object not initialized with data file!\n");

    double Time = time.JD();
    if (Time < H1.data.timeData[0] || Time > H1.data.timeData[1])
     throw TimeNotCoveredByFileEXC();

    // Initialize local coefficient array.
    for (i = 0; i < 50; i++) A[i] = 0.0;

    // Determine if a new record needs to be input
    // and read it, if needed.
    // Attention: This function may generate a FileErrorEXC - exception!
    if (Time < T_beg || Time > T_end) getCoefficients(Time);


    // Read the coefficients from the binary record.
    C = R1.coeffPtr[body][0] - 1; // Coefficient array entry point
    N = R1.coeffPtr[body][1];     // Number of coefficients
    G = R1.coeffPtr[body][2];     // Granules in current record

/*
    // Debug print (optional)
    if (Debug) {
      cout << "\n  In: Interpolate_Nutation\n";
      cout << "\n  Target = " << setw(2) << Target;
      cout << "\n  C      = " << setw(4) << C << " (before)";
      cout << "\n  N      = " << setw(4) << N;
      cout << "\n  G      = " << setw(4) << G << endl;
    }
*/


    // Compute the normalized time, then load the Tchebeyshev coefficients
    // into array A[]. If T_span is covered by a single granule this is easy.
    // If not, the granule that contains the interpolation time is found, and
    // an offset from the array entry point for the nutation angles is used
    // to load the coefficients.

    if (G == 1) {
      Tc = 2.0 * (Time - T_beg) / T_span - 1.0;
      for (i = C; i < (C + 3 * N); i++) A[i - C] = Coeff_Array[i];
    } else if (G > 1) {
      T_sub = T_span / ((double)G); // Compute subgranule interval

      for (j = G; j > 0; j--) {
        T_break = T_beg + ((double)j - 1) * T_sub;
        if (Time > T_break) {
          T_seg = T_break;
          offset = j - 1;
          break;
        }
      }

      Tc = 2.0 * (Time - T_seg) / T_sub - 1.0;
      C = C + 3 * offset * N;

      for (i = C; i < (C + 3 * N); i++) A[i - C] = Coeff_Array[i];
    } else {                        // Something has gone terribly wrong
      cout << "\n Number of granules must be >= 1: check header data.\n";
    }
/*
    // Debug print (optional)
    if (Debug) {
      cout << "\n  C      = " << setw(4) << dec << C << " (after)";
      cout << "\n  offset = " << setw(4) << dec << offset;
      cout << "\n  Time   = " << setw(12) << setprecision(7) << Time;
      cout << "\n  T_sub  = " << setw(12) << setprecision(7) << T_sub;
      cout << "\n  T_seg  = " << setw(12) << setprecision(7) << T_seg;
      cout << "\n  Tc     = " << setw(12) << setprecision(7) << Tc;
      cout << "\n  Array Coefficients:\n";
      for (i = 0; i < 3 * N; i++) {
        cout << "\n  A[" << setw(2) << i << "] = " << setw(22) << setprecision(15)
             << A[i];
      }
      cout << endl << endl;
    }
*/

    // Compute interpolated the nutation.
    for (i = 0; i < 2; i++) {
      Cp[0] = 1.0;                // Begin polynomial sum
      Cp[1] = Tc;
      sum[i] = A[i * N] + A[1 + i * N] * Tc;

      for (j = 2; j < N; j++) {   // Finish it
        Cp[j] = 2.0 * Tc * Cp[j - 1] - Cp[j - 2];
        sum[i] = sum[i] + A[j + i * N] * Cp[j];
      }
      Nutation[i] = sum[i];
    }

    return MathLib::Vector(Nutation[0], Nutation[1], Nutation[2]);
  }



  //---------------------------------------------------------------------------
  MathLib::Vector DE405::getLibration(const DE405StateIndex body,
                                      const SpaceLib::CTime& time)
   throw(DE405::IndexOutOfRangeEXC,
         DE405::TimeNotCoveredByFileEXC,
         DE405::NotInitializedEXC,
         FileErrorEXC) {


    double A[50], Cp[50], sum[3], T_break, T_seg=0, T_sub=0, Tc=0;
    double Libration[3];
    int i, j;
    long int C, G, N, offset = 0;

    // This function only computes librations.
    if (body != DE405_EM_BARY)
     throw IndexOutOfRangeEXC("Libration only for earth -moon barycenter!\n");

    if (!initialized)
     throw NotInitializedEXC("Object not initialized with data file!\n");

    double Time = time.JD();
    if (Time < H1.data.timeData[0] || Time > H1.data.timeData[1])
     throw TimeNotCoveredByFileEXC();

    for (i = 0; i < 50; i++) A[i] = 0.0;

    // Determine if a new record needs to be input (if so, get it).
    if (Time < T_beg || Time > T_end) getCoefficients(Time);

    // Read the coefficients from the binary record.
    C = R1.libratPtr[0] - 1;      // Coefficient array entry point
    N = R1.libratPtr[1];          // Number of coefficients
    G = R1.libratPtr[2];          // Granules in current record

/*
    // Debug print (optional)
    if (Debug) {
     cout << "\n  In: Interpolate_Libration\n";
     cout << "\n  Target = " << setw(2) << Target;
     cout << "\n  C      = " << setw(4) << C;
     cout << "\n  N      = " << setw(4) << N;
     cout << "\n  G      = " << setw(4) << G;
    }
*/


    //  Compute the normalized time, then load the Tchebeyshev coefficients
    //  into array A[]. If T_span is covered by a single granule this is easy.
    //  If not, the granule that contains the interpolation time is found, and
    //  an offset from the array entry point for the libration angles is used
    //  to load the coefficients.

    if (G == 1) {
     Tc = 2.0 * (Time - T_beg) / T_span - 1.0;
     for (i = C; i < (C + 3 * N); i++)
     A[i - C] = Coeff_Array[i];
    } else if (G > 1) {
      T_sub = T_span / ((double)G); // Compute subgranule interval

      for (j = G; j > 0; j--) {
        T_break = T_beg + ((double)j - 1) * T_sub;
        if (Time > T_break) {
          T_seg = T_break;
          offset = j - 1;
          break;
        }
      }

      Tc = 2.0 * (Time - T_seg) / T_sub - 1.0;
      C = C + 3 * offset * N;

      for (i = C; i < (C + 3 * N); i++) A[i - C] = Coeff_Array[i];
    } else {                        // Something has gone terribly wrong
      cout << "\n Number of granules must be >= 1: check header data.\n";
    }

/*
    // Debug print (optional)
    if (Debug) {
      cout << "\n  C      = " << setw(4) << dec << C << " (after)";
      cout << "\n  offset = " << setw(4) << dec << offset;
      cout << "\n  Time   = " << setw(12) << setprecision(7) << Time;
      cout << "\n  T_sub  = " << setw(12) << setprecision(7) << T_sub;
      cout << "\n  T_seg  = " << setw(12) << setprecision(7) << T_seg;
      cout << "\n  Tc     = " << setw(12) << setprecision(7) << Tc;
      cout << "\n  Array Coefficients:\n";
      for (i = 0; i < 3 * N; i++) {
        cout << "\n  A[" << setw(2) << i << "] = " << setw(22) << setprecision(15)
             << A[i];
      }
      cout << "\n\n";
    }
*/

    // Compute interpolated the libration.
    for (i = 0; i < 3; i++) {
      Cp[0] = 1.0;                // Begin polynomial sum
      Cp[1] = Tc;
      sum[i] = A[i * N] + A[1 + i * N] * Tc;
      for (j = 2; j < N; j++) {   // Finish it
        Cp[j] = 2.0 * Tc * Cp[j - 1] - Cp[j - 2];
        sum[i] = sum[i] + A[j + i * N] * Cp[j];
      }
      Libration[i] = sum[i];
    }
    return MathLib::Vector(Libration[0], Libration[1], Libration[2]);
  }


  // returns a constant stored in the data file
  //---------------------------------------------------------------------------
  double DE405::getConstant(const char* Name)
    throw(DE405::NotInitializedEXC) {

    if (!initialized)
      throw NotInitializedEXC("Object not initialized with data file!\n");

    int i;
    for (i = 0; i < R1.numConst; i++) {
      if (memcmp(R1.constName[i], Name, 6) == 0)
      break;
    }

    if (i != R1.numConst) return R2.constValue[i];
    else return -1.0;
  }


  // This function is used by the member functions to read an array of
  // Tchebeychev coefficients from a binary ephemeris data file.
  //---------------------------------------------------------------------------
  void DE405::getCoefficients(const double time) throw(FileErrorEXC) {

    double T_delta = 0.0;
    int Offset = 0;

    if (time < T_beg) {           // Compute backwards location offset
      T_delta = T_beg - time;
      Offset = (int)-ceil(T_delta / T_span);
    }

    if (time > T_end) {           // Compute forewards location offset
      T_delta = time - T_end;
      Offset = (int)ceil(T_delta / T_span);
    }

    //  Retrieve ephemeris data from new record.
    _data_file.seekg((static_cast<std::streamoff>(Offset) - 1) * static_cast<std::streamoff>(ARRAY_SIZE) * static_cast<std::streamoff>(sizeof(double)), ios_base::cur);
    if (_data_file.fail())
      throw FileErrorEXC("getCoefficients: Error during seekg-operation!\n");

    _data_file.read((char *)&Coeff_Array, ARRAY_SIZE * sizeof(double));
    if (_data_file.fail())
      throw FileErrorEXC("getCoefficients: Error during read-operation!\n");

    T_beg = Coeff_Array[0];
    T_end = Coeff_Array[1];
    T_span = T_end - T_beg;

  }


  // returns the gravitational parameter of the desired body
  //---------------------------------------------------------------------------
  double DE405::getGM(const DE405StateIndex body)
   throw(DE405::IndexOutOfRangeEXC, DE405::NotInitializedEXC) {

    double value = 0.0;
    double temp;

    switch (body) {
      case DE405_MERCURY:
        value = getConstant("GM1   ");
        break;
      case DE405_VENUS:
        value = getConstant("GM2   ");
        break;
      case DE405_EARTH:
        if ((temp = getConstant("EMRAT ")) == -1.0) break;
        if ((value = getConstant("GMB   ")) == -1.0) break;
        value = value - value / (1.0 + temp);
        break;
      case DE405_MARS:
        value = getConstant("GM4   ");
        break;
      case DE405_JUPITER:
        value = getConstant("GM5   ");
        break;
      case DE405_SATURN:
        value = getConstant("GM6   ");
        break;
      case DE405_URANUS:
        value = getConstant("GM7   ");
        break;
      case DE405_NEPTUNE:
        value = getConstant("GM8   ");
        break;
     case DE405_PLUTO:
       value = getConstant("GM9   ");
       break;
     case DE405_MOON:
       if ((temp = getConstant("EMRAT ")) == -1.0) break;
       if ((value = getConstant("GMB   ")) == -1.0) break;
       value = value / (1.0 + temp);
       break;
     case DE405_SUN:
       value = getConstant("GMS   ");
       break;
     case DE405_SS_BARY:
       return 0.0;
     case DE405_EM_BARY:
       value = getConstant("GMB   ");
       break;
     case DE405_NUTATIONS:
     case DE405_LIBRATIONS:
     case DE405_NO_BODY:
       value = 0.0;
       break;
  }

  if (value != -1.0)
    return value;
  else
    return 0.0;

  }


}
