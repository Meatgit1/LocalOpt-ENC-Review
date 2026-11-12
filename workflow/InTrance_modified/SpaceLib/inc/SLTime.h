/*!
\file SLTime.h
\version 1.7
\brief A class for time calculations in astrodynamical problems (CTime).
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\date 01.11.2004 (AO)

\b Purpose:
   The class CTime simplifies the time calculations in
   astrodynamical problems. Especially the conversion from a given
   date and time to a julian or modified julian date is capsuled
   in the class for convienent use. The possibility of using
   different output formats and units for calculating operations
   completes the class.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.

\b References:
   (1) Montenbruck: Grundlagen der Ephemeridenrechnung (6. Auflage). Sterne und Weltraum, 2001
   (2) Fränz, Harper: Heliospheric coordinate systems. Planetary and Space Science 50 (2002), pp. 217-233, 2001

\b Copyright:
   (c) 2003-2004 Andreas Ohndorf, Bernd Dachwald
*/

#ifndef SLTIME_H
#define SLTIME_H

#include <Global.h>

namespace SpaceLib {

  // important julian dates
  const double JD_J1900 = 2415020.0;            ///< julian date of J1900 (01.01.1900 12:00 TDB)
  const double JD_J1950 = 2433282.5;            ///< julian date of J1950 (01.01.1950 00:00 TDB)
  const double JD_J2000 = 2451545.0;            ///< julian date of J2000 (01.01.2000 12:00 TDB)
  const double JD_B1950 = 2433282.42345905;     ///< julian date of Bessel epoch 1950
  const double JD_MJD_0 = 2400000.5;            ///< MJD 0.0 (julian date of 17.11.1858, 0h UT)
  const double MJD_J2000 = JD_J2000 - JD_MJD_0; ///< MJD of J2000 (01.01.2000 12:00 TDB)

  // conversion factors
  const double JCY2JYR = 100.0;              ///< conversion factor: julian centuries => julian years (obvious)
  const double JYR2DAY = 365.25;             ///< conversion factor: julian years => days
  const double DAY2SEC = 24.0 * 60.0 * 60.0; ///< conversion factor: days => seconds
  const double SEC2DAY = 1.0 / DAY2SEC;      ///< conversion factor: seconds => days
  const double DAY2JYR = 1.0 / JYR2DAY;      ///< conversion factor: julian days => years
  const double JYR2JCY = 1.0 / JCY2JYR;      ///< conversion factor: julian year => julian century
  const double JCY2DAY = JCY2JYR * JYR2DAY;  ///< conversion factor: julian centuries => julian days
  const double DAY2JCY = 1.0 / JCY2DAY;      ///< conversion factor: julian days => julian centuries
  const double JYR2SEC = JYR2DAY * DAY2SEC;  ///< conversion factor: julian years => seconds
  const double SEC2JYR = 1.0 / JYR2SEC;      ///< conversion factor: julian years => seconds
  const double SEC2JCY = SEC2JYR * JYR2JCY;  ///< conversion factor: seconds => julian centuries

  /** Calculates the gregorian date out of the julian date.
      Algorithm from: Montenbruck: Grundlagen der Ephemeridenrechnung (6. Auflage), 2001 */
  void JD2date(double jd, int &year, int &month, int &day, double &hourUT);

  /** Calculates the gregorian date out of the modified julian date.
      Algorithm from: Montenbruck: Grundlagen der Ephemeridenrechnung (6. Auflage), 2001 */
  void MJD2date(double mjd, int &year, int &month, int &day, double &hourUT);

  /** Calculates the julian date out of the gregorian date.
    Algorithm from: Montenbruck: Grundlagen der Ephemeridenrechnung (6. Auflage), 2001 */
  double date2JD(int year, int month, int day, double hourUT);

  /** Calculates the modified julian date out of the gregorian date. */
  double date2MJD(int year, int month, int day, double hourUT);

  /** Calculates the approximation for TT-UT time difference.
    Algorithm from: Montenbruck: Grundlagen der Ephemeridenrechnung (6. Auflage), 2001 */
  double deltaT(double jd);

  /** Estimates the number of leap seconds between TAI and UTC */
  int get_leapseconds(int y, int m, int d);

  /** \brief Output format of CTime class
      This is the enumeration in which the output format information
      is stored. */
  enum TimeFormat
  {
    TimeFormat_JD,       ///< julian date
    TimeFormat_MJD,      ///< modified julian date
    TimeFormat_HDMY_UT,  ///< hourUT.day.month.year (universal time)
    TimeFormat_HDMY_TAI, ///< hourTAI.day.month.year (atomic time)
    TimeFormat_HDMY_UTC, ///< hourUTC.day.month.year (corrected universal time)
    TimeFormat_HDMY_TT   ///< hourTT.day.month.year (terrestrial time)
  };

  /** \brief Time unit for incremental/decremental operators.
      The enumeration describes the time unit, by which a CTime
      object should be changed when the incremental (++) or
      decremental operator (--) is called. */
  enum TimeUnit
  {
    TimeUnit_second,  ///< one second
    TimeUnit_minute,  ///< one minute
    TimeUnit_hour,    ///< one hour
    TimeUnit_day,     ///< one day
    TimeUnit_Jmonth,  ///< one julian month = julian year / 12
    TimeUnit_Jyear,   ///< one julian year
    TimeUnit_Jcentury ///< one julian century
  };

  /** \brief Class CTime
  *
      CTime can be used to handle different time systems, e.g. julian date,
      modified julian date, normal date, julian centuries since a given date
      or simple mathematical operations like addition and subtraction. */
  class CTime
  {
  private:

    double m_value;  ///< time as a julian date
    TimeFormat   m_format; ///< time out put format
    TimeUnit   m_unit;   ///< time unit for incremental and decremental operators

  public:

    /// Class for CTime exceptions.
    class TimeEXC : public SLEXC {
    public:
      TimeEXC(void) : SLEXC() {};
      /** Tbd. */
      TimeEXC(const char* s) : SLEXC(s) {};
      virtual ~TimeEXC(void) throw() {};
    };

    /** \brief Default constructor
        Initializes the CTime object with the julian date of J2000. */
    CTime(void);

    /** \brief Copy constructor
        Initializes the CTime object with another CTime object.
        \param t CTime object template */
    CTime(const CTime& t);

    /** \brief Julian Date - constructor
        Initializes the CTime object with a julian date.
        \param _JD julian date */
    CTime(double _JD);

    /** \brief Normal date-constructor
        Initializes the CTime object with a gregorian date (year, month, day, hour) (coordinated universal time, UTC).
        The hour value has to be a double value, so that eventually remaining minutes
        and seconds have to be expressed as the fraction of an hour.
        \param year year
        \param month month
        \param day day
        \param hourUTC hourUTC = hours + (minutes/60.0) + (seconds/3600.0) */
    CTime(int year, int month, int day, double hourUTC);

    /** \brief Destructor
        The CTime destructor does nothing special. */
    ~CTime();

    /** \brief Sets the CTime object to a julian date.
        Sets the CTime object to the given point in time.
        \param _JD julian date
        \exception TimeEXC julian date < 0.0 and therefore not defined */
    void setTime(double _JD) throw(TimeEXC);

    /** \brief Sets the CTime object with a gregorian date (universal time, UT).
        Sets the CTime object to the given point in time.
        The hour value has to be a double value, so that eventually remaining minutes
        and seconds have to be expressed as the fraction of an hour.
        \param year year
        \param month month
        \param day day
        \param hourUT hourUT = hours + (minutes/60.0) + (seconds/3600.0)
        \exception TimeEXC at least one of the input parameters is invalid */
    void setTimeUT(int year, int month, int day, double hourUT) throw(TimeEXC); /// \todo why UT and UTC above?

    /** \brief Year/Month/Day/Hour representation of date (universal time, UT)
        The parameters given by reference will be updated according to
        the information stored in the CTime object. The hourUT value is a decimal
        one so the minutes and seconds are expressed as the fraction of one.
        \return year year
        \return month The month of date.
        \return day The day of date.
        \return hourUT hour in UT as decimal value */
    void dateUT(int& year, int& month, int& day, double& hourUT) const;

    /** \brief Year/Month/Day/Hour representation of date (terrestrial time, TT)
        The parameters given by reference will be updated according to
        the information stored in the CTime object. The hourTT value is a decimal
        one so the minutes and seconds are expressed as the fraction of one.
        \return year year
        \return month The month of date.
        \return day The day of date.
        \return hourTT hour in TT as decimal value */
    void dateTT(int& year, int& month, int& day, double& hourTT) const;

    /** \brief Year/Month/Day/Hour representation of date (international atomic time, TAI)
        The parameters given by reference will be updated according to
        the information stored in the CTime object. The hourTAI value is a decimal
        one so the minutes and seconds are expressed as the fraction of one.
        \return year year
        \return month The month of date.
        \return day The day of date.
        \return hourTAI hour in TAI as decimal value */
    void dateTAI(int& year, int& month, int& day, double& hourTAI) const;

    /** \brief Year/Month/Day/Hour representation of date (international atomic time, UTC)
        The parameters given by reference will be updated according to
        the information stored in the CTime object. The hourUTC value is a decimal
        one so the minutes and seconds are expressed as the fraction of one.
        \return year year
        \return month The month of date.
        \return day The day of date.
        \return hourUTC hour in UTC as decimal value */
    void dateUTC(int& year, int& month, int& day, double& hourUTC) const;

    /** \brief Returns the time stored in the CTime object as julian date.
        Returns the time stored in the CTime object as a julian date.
        \return julian date stored in CTime object.
        \sa MJD */
    double JD(void) const;

    /** \brief Returns the time stored in the CTime object as modified julian date.
        Returns the time stored in the CTime object as a modified julian date.
        \return modified julian date stored in the CTime object.
        \sa JD */
    double MJD(void) const;

    /** \brief Julian centuries since given date.
        Calculates the decimal value of the Julian Centuries which have passed
        since the given date. This particular date has to be expressed as an
        CTime object. The expression 'since' isn't meant literally because
        the function also allows values that are in the future.
        Therefore, the returned value might also be negativ.
        \return julian centuries that have passed since (decimal value) */
    double JC(const CTime& since) const;

    /** \brief Assign operator.
        The whole content of the given object will be copied (the data and the
        output format). The ouput of a CTime object after an assignment
        operation will therefore depend on the output settings of the object that
        has been given (as input parameter).
        \return reference to the CTime object */
    CTime& operator=(const CTime&);

    /** \brief Addition operator.
        This operator takes the Julian Date from the given CTime object (as
        input parameter) and adds this value to its own value. The result
        is returned as an CTime object.
        \param right CTime object to be added
        \return sum (as CTime object)
    */
    CTime operator+(const CTime& right);

    /** \brief Subtraction operator.
        This operator takes the Julian Date from the given CTime object (as
        input parameter) and subtracts this value from its own value. The result
        is returned as an CTime object.
        \param right CTime object to be subtracted
        \return difference (as CTime object)
    */
    CTime operator-(const CTime& right);

    /** \brief Addition/Assign operator.
        This operator takes the Julian Date from the given CTime object (as
        input parameter) and adds this value to its own value. Then the CTime
        object returns a reference to itself.
        \param right CTime object to be added
        \return reference to the modified CTime object */
    CTime& operator+=(const CTime& right);

    /** \brief Subtraction/Assign operator.
        This operator takes the Julian Date from the given CTime object (as
        input parameter) and subtracts this value from its own value. Then the CTime
        object returns a reference to itself.
        \param right CTime object to be subtracted
        \return reference to the modified CTime object */
    CTime& operator-=(const CTime& right);

    /** \brief Gets the time unit for inc/dec operations.
        This member function returns the time unit used for the decrement
        and increment operations.
        Example: If the return value of this member function is TimeUnit_hour and
                 the ++ operator is called, the content of the CTime object
                 will be increased by 1/24 julian days and therefore by one hour.
        \return time unit for inc/dec operations
        \sa setOpUnit */
    TimeUnit getOpUnit(void) const;

    /** \brief Sets the time unit for inc/dec operations.
        This member function sets the time unit used for the decrement
        and increment operations.
        \param _u unit for inc/dec operations
        \sa setOpUnit */
    void setOpUnit(TimeUnit _u);

    /** \brief The increment operator (prefix). */
    CTime& operator++();

    /** \brief The decrement operator (prefix). */
    CTime& operator--();

    /** \brief The increment operator (postfix). */
    CTime operator++(int);

    /** \brief The decrement operator (postfix). */
    CTime operator--(int);

    /** \brief Get the actual ouput settings of the CTime object.
        \return the object's output settings
        \sa setOutputFormat */
    TimeFormat getOutputFormat(void) const;

    /** \brief Get the ouput settings of the CTime object.
        The given parameter will be stored and used for further output operations.
        \param _format new output settings
        \sa getOutputFormat  */
    void setOutputFormat(TimeFormat _format);

    friend std::ostream& operator<<(std::ostream&, const CTime&);

  };  // end of class CTime

} // end of namespace SpaceLib

#endif // of SLTIME_H
