/*!
\file DE405.h
\version 1.7
\brief A class to access the JPL ephemeris data files DE405 (DE405).
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\date 15.11.2004 (AO)

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
   (c) 2004 Andreas Ohndorf, Bernd Dachwald
*/

#ifndef DE405_H
#define DE405_H

#include <SLTime.h>
#include <MLVector.h>
#include <SpaceConst.h>

namespace SpaceLib {

  const int EPHEMERIS  = 405;   ///< tbd
  const int ARRAY_SIZE = 1018;  ///< the data file's record size

  /** Tbd. */
  enum DE405StateIndex {
    DE405_MERCURY    = 0,  ///< Mercury
    DE405_VENUS      = 1,  ///< Venus
    DE405_EARTH      = 2,  ///< Earth-Moon Barycenter
    DE405_MARS       = 3,  ///< Mars
    DE405_JUPITER    = 4,  ///< Jupiter
    DE405_SATURN     = 5,  ///< Saturn
    DE405_URANUS     = 6,  ///< Uranus
    DE405_NEPTUNE    = 7,  ///< Neptune
    DE405_PLUTO      = 8,  ///< Pluto
    DE405_MOON       = 9,  ///< Moon relative to geocenter
    DE405_SUN        = 10, ///< tbd
    DE405_SS_BARY    = 11, ///< tbd
    DE405_EM_BARY    = 12, ///< tbd
    DE405_NUTATIONS  = 13, ///< tbd
    DE405_LIBRATIONS = 14, ///< tbd
    DE405_NO_BODY    = 15  ///< tbd
  };

  /** Tbd. */
  struct stateData
  {
    double Position[3]; ///< tbd
    double Velocity[3]; ///< tbd
  };

  typedef struct stateData stateType; ///< tbd

  /** Define the content of binary header records. */
  struct DE405RecOneData {
    char     label[3][84];      ///< tbd
    char     constName[400][6]; ///< tbd
    double   timeData[3];       ///< tbd
    long int numConst;          ///< tbd
    double   AU;                ///< tbd
    double   EMRAT;             ///< tbd
    long int coeffPtr[12][3];   ///< tbd
    long int DENUM;             ///< tbd
    long int libratPtr[3];      ///< tbd
  };

  /** Tbd. */
  struct DE405RecTwoData
  {
    double constValue[400]; ///< tbd
  };

  typedef struct DE405RecOneData DE405RecOneType; ///< tbd
  typedef struct DE405RecTwoData DE405RecTwoType; ///< tbd

  /** Define binary header record formats. */
  struct DE405HeaderOne
  {
    DE405RecOneType data;                                            ///< tbd
    char pad[ARRAY_SIZE * sizeof(double) - sizeof(DE405RecOneType)]; ///< tbd
  };

  /** Tbd. */
  struct DE405HeaderTwo
  {
    DE405RecTwoType data;                                            ///< tbd
    char pad[ARRAY_SIZE * sizeof(double) - sizeof(DE405RecTwoType)]; ///< tbd
  };

  typedef struct DE405HeaderOne DE405HeadOneType; ///< tbd
  typedef struct DE405HeaderTwo DE405HeadTwoType; ///< tbd


  /** Exceptions for error handling in C++ style. */
  class FileErrorEXC : public SLEXC {
   public:
   /** Tbd. */
            FileErrorEXC(void) : SLEXC() {};
            FileErrorEXC(const std::string& em):SLEXC(em){};
    virtual ~FileErrorEXC(void) throw() {};
  };

  /** Tbd. */
  class DE405 {

   public:

    /** Tbd. */
    class IndexOutOfRangeEXC : public SLEXC {
     public:
       /** Tbd. */
               IndexOutOfRangeEXC(const std::string& em):SLEXC(em) {};
       virtual ~IndexOutOfRangeEXC(void) throw() {};
    };

    /** Tbd. */
    class TimeNotCoveredByFileEXC : public SLEXC {
      public:
              TimeNotCoveredByFileEXC(void) {};
              TimeNotCoveredByFileEXC(const std::string& em):SLEXC(em) {};
      virtual ~TimeNotCoveredByFileEXC(void) throw() {};
    };

    /** Tbd. */
    class NotInitializedEXC : public SLEXC {
      public:
        /** Tbd. */
                NotInitializedEXC(const std::string& em):SLEXC(em){};
        virtual ~NotInitializedEXC(void) throw() {};
    };

    // constructors
    DE405(void);
    /** Tbd. */
    DE405(const std::string data_file_name);

    // destructor
    ~DE405(void);

    /** Tbd. */
    void setDataFile(const std::string data_file_name);

    /** Tbd. */
    bool initialize(void) throw(FileErrorEXC);

    /** Tbd. */
    MathLib::Vector getState(const DE405StateIndex body,
                             const SpaceLib::CTime& time)
     throw(IndexOutOfRangeEXC,TimeNotCoveredByFileEXC, NotInitializedEXC, FileErrorEXC);

    /** Tbd. */
    MathLib::Vector getNutation(const DE405StateIndex body,
                                const SpaceLib::CTime& time)
     throw(IndexOutOfRangeEXC,TimeNotCoveredByFileEXC, NotInitializedEXC, FileErrorEXC);

    /** Tbd. */
    MathLib::Vector getLibration(const DE405StateIndex body,
                                 const SpaceLib::CTime& time)
     throw(IndexOutOfRangeEXC,TimeNotCoveredByFileEXC, NotInitializedEXC, FileErrorEXC);

    /** Tbd. */
    double getConstant(const char* Name) throw(NotInitializedEXC);

    /// Tbd.
    double getGM(const DE405StateIndex body)
     throw(IndexOutOfRangeEXC, NotInitializedEXC);

   private:
    std::ifstream _data_file;      ///< tbd
    std::string   _data_file_name; ///< tbd

    DE405HeadOneType H1;            ///< tbd
    DE405HeadTwoType H2;            ///< tbd
    DE405RecOneType R1;             ///< tbd
    DE405RecTwoType R2;             ///< tbd
    double Coeff_Array[ARRAY_SIZE]; ///< tbd
    double T_beg;                   ///< tbd
    double T_end;                   ///< tbd
    double T_span;                  ///< tbd

    bool initialized; ///< tbd

    /** Tbd. */
    void getCoefficients(const double time) throw(FileErrorEXC);

  }; // of class DE405

} // of namespace SpaceLib

#endif
