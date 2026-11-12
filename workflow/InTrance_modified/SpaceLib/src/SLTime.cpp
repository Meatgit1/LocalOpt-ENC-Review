/*!
\file SLTime.cpp
\version 1.7
\brief A class for time calculations in astrodynamical problems (CTime).
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\date 06.01.2005 (AO)

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
   (c) 2003-2005 Andreas Ohndorf, Bernd Dachwald

\todo Improvement/Enhancement of the implementation of DT=TT-UT1 according to
      the information of the IERS files, polynomial approximation formulas
      and the stored (dlr/literature/topic/time) website
      http://www.phys.uu.nl/~vgent/astro/deltatime.htm
\todo finish documentation on the "date" member functions
*/

//SpaceLib
#include <configSpaceLib.h>

#ifdef DEBUG
// #ifndef SLTIME_DEBUG
//  #define SLTIME_DEBUG 1
// #endif //SLTIME_DEBUG
#endif // DEBUG

//STL
#include <iomanip>

//standard C library
#include <math.h>

//SpaceLib
#include <SLTime.h>

namespace SpaceLib
{
// ----------------------------------------------------------------------------
// Table with the number of leap seconds between TAI and UTC and the date from
//  when they are to use. Data from: http://maia.usno.navy.mil
//  valid for dates after 01.01.1972 (JD 2441317.5, MJD 41317.0)

  /** Tbd. */
  struct LeapSecond
  {
    int year;     ///< tbd
    int month;    ///< tbd
    int day;      ///< tbd
    int leapsecs; ///< tbd

    /** Tbd. */
    LeapSecond(int y, int m, int d, int lp) : year(y), month(m), day(d), leapsecs(lp) {}
  };

  /** Tbd. */
  const LeapSecond ARRAY_LS[] =
  {
    LeapSecond(1972, 1, 1, 10),
    LeapSecond(1972, 7, 1, 11),
    LeapSecond(1973, 1, 1, 12),
    LeapSecond(1974, 1, 1, 13),
    LeapSecond(1975, 1, 1, 14),
    LeapSecond(1976, 1, 1, 15),
    LeapSecond(1977, 1, 1, 16),
    LeapSecond(1978, 1, 1, 17),
    LeapSecond(1979, 1, 1, 18),
    LeapSecond(1980, 1, 1, 19),
    LeapSecond(1981, 7, 1, 20),
    LeapSecond(1982, 7, 1, 21),
    LeapSecond(1983, 7, 1, 22),
    LeapSecond(1985, 7, 1, 23),
    LeapSecond(1988, 1, 1, 24),
    LeapSecond(1990, 1, 1, 25),
    LeapSecond(1991, 1, 1, 26),
    LeapSecond(1992, 7, 1, 27),
    LeapSecond(1993, 7, 1, 28),
    LeapSecond(1994, 7, 1, 29),
    LeapSecond(1996, 1, 1, 30),
    LeapSecond(1997, 7, 1, 31),
    LeapSecond(1999, 1, 1, 32)
  };

// ----------------------------------------------------------------------------
// Coefficent table for calculation of TT-UT difference.
// Data from: Montenbruck: Grundlagen der Ephemeridenrechnung (6. Auflage), 2001

  /** Tbd. */
  const double ARRAY_TT2UT_COEFF[][6] =
  {
    {1825.0,  10.4, - 80.8,  413.9, - 572.3, 1.75},
    {1850.0,   6.6,   46.3, -358.4,    18.8, 1.50},
    {1875.0, - 3.9, - 10.8, -166.2,   867.4, 1.25},
    {1900.0, - 2.6,  114.1,  327.5, -1467.4, 1.00},
    {1925.0,  24.2,  - 6.3, -  8.2,   483.4, 0.75},
    {1950.0,  29.3,   32.5, -  3.8,   550.7, 0.50},
    {1975.0,  45.3,  130.5, -570.5,  1516.7, 0.25},
    {2000.0,   0.0,    0.0,    0.0,     0.0,  0.0}
  };

// ======================= auxiliary functions ================================

  void JD2date(double jd, int &year, int &month, int &day, double &hourUT)
  {
    double b, c;

    const double a = floor(jd + 0.5);
    if (a < 2299161.0)
      c = a + 1524.0;
    else
    {
      b = floor((a - 1867216.25) / 36524.25);
      c = a + b - floor(b/4.0) + 1525.0;
    }
    const double d = floor((c - 122.1) / 365.25);
    const double e = floor(365.25 * d);
    const double f = floor((c - e) / 30.6001);

    day   = int (c - e - floor(30.6001 * f) + (jd + 0.5 - a));
    month = int (f - 1.0 - 12.0 * floor(f/14.0));
    year  = int (d - 4715.0 - floor((7.0 + month) / 10.0));

    hourUT = 24.0 * (jd - 0.5 - floor(jd - 0.5));
  }

// ----------------------------------------------------------------------------

  void MJD2date(double mjd, int &year, int &month, int &day, double &hourUT)
  {
    JD2date(mjd + JD_MJD_0, year, month, day, hourUT);
  }

// ----------------------------------------------------------------------------

  double date2JD(int year, int month, int day, double hourUT)
  {
    double b, jd;

    if (month <= 2) {
      month += 12;
      --year;
    }

    if ((10000L * year + 100L * month + day) <= 15821004L) // julian calendar
      b = -2;
    else // gregorian calendar
      b = floor(year/400.0) - floor(year/100.0);

    if (year > 0)
      jd = floor(365.25 * year) + floor(30.6001 * (month + 1)) + b + 1720996.5 + day + hourUT / 24.0;
    else
      jd = int(365.25 * year - 0.75) + floor(30.6001 * (month + 1)) + b + 1720996.5 + day + hourUT / 24.0;

    return jd;
  }

// ----------------------------------------------------------------------------

  double date2MJD(int year, int month, int day, double hourUT)
  {
    return date2JD(year, month, day, hourUT) - JD_MJD_0;
  }

// ----------------------------------------------------------------------------

  double deltaT(double jd)
  {
    int          y, m, d;
    double       h;
    unsigned int i;
    double       t, T;

    JD2date(jd, y, m, d, h); // get calender date with hour in UT
    // find index of time span to use
    for (i = 0; i < sizeof ARRAY_TT2UT_COEFF / sizeof *ARRAY_TT2UT_COEFF; i++)
      if (ARRAY_TT2UT_COEFF[i][0] > y)
        break;

    // is the given date covered by the coefficient table?
    if (i == 0 || ARRAY_TT2UT_COEFF[i][0] < y)
      throw CTime::TimeEXC("deltaT::TT2UTArrayFault: time not covered by coefficient table");

    i -= 1; // got to correct entry in table

    // julian centuries since 2000/01/01 12:00 UT
    T = (jd - JD_J2000) * DAY2JCY;

    t = T + ARRAY_TT2UT_COEFF[i][5]; // polynomial variable

    return (  ARRAY_TT2UT_COEFF[i][1]
            + ARRAY_TT2UT_COEFF[i][2] * t
            + ARRAY_TT2UT_COEFF[i][3] * t*t
            + ARRAY_TT2UT_COEFF[i][4] * t*t*t);
  }

// ----------------------------------------------------------------------------

  int get_leapseconds(int y, int m, int d)
  {
    unsigned int i;
    int          v1 = y * 365 + m * 30 + d;
    int          v2 = 0;

    for (i = 0; i < sizeof ARRAY_LS / sizeof *ARRAY_LS; i++)
    {
      v2 = ARRAY_LS[i].year * 365 + ARRAY_LS[i].month * 30 + ARRAY_LS[i].day;
      if (v2 >= v1)
        break;
    }
    if (v2 != v1 && i == 0)
      throw CTime::TimeEXC("LeapSecondsArrayFault: time not covered by leap seconds table");

    return ARRAY_LS[v2 == v1 ? i : i - 1].leapsecs;
  }

// =========================== member functions ==============================
  CTime::CTime(void)
  {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::CTime(void)=>";
    std::cout << JD_J2000 << std::endl;
#endif

    m_value = JD_J2000;
    m_format = TimeFormat_JD;
    m_unit = TimeUnit_day;
  }


// ----------------------------------------------------------------------------
  CTime::CTime(const CTime & t) {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::CTime(CTime &)=>";
    std::cout << "&.JD()=" << m_value << ",outputFormat=" << m_format << std::endl;
#endif

    m_value = t.JD();
    m_format = t.getOutputFormat();
    m_unit = t.getOpUnit();

  }


// ----------------------------------------------------------------------------
  CTime::CTime(double _JD) {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::CTime(double)=>";
    std::cout << _JD << std::endl;
#endif

    try {
      setTime(_JD);
    }
    catch(TimeEXC & _exc) {
      std::cout << _exc.what() << std::endl;
    }

    m_value = _JD;
    m_format = TimeFormat_JD;
    m_unit = TimeUnit_day;
  }


// ----------------------------------------------------------------------------
  CTime::CTime(int year, int month, int day, double hourUT) {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::CTime(int,int,int,double)=>"
      << hourUT << " UT;" << day << "." << month << "." << year << std::endl;
#endif
    try {
      setTimeUT(year, month, day, hourUT);
    }
    catch(TimeEXC & _exc) {
    std::cout << _exc.what() << std::endl;
    }
    m_format = TimeFormat_JD;
    m_unit = TimeUnit_day;
  }


// ----------------------------------------------------------------------------
  CTime::~CTime(void) {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::~CTime()=>" << std::endl;
#endif
  }


// ----------------------------------------------------------------------------
  void CTime::setTime(double _JD) throw(CTime::TimeEXC) {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::setTime(double)=>" << _JD << std::endl;
#endif

    if (_JD < 0.0)
      throw TimeEXC("CTime::setTime(double): julian date below zero!");

    m_value = _JD;
  }


// ----------------------------------------------------------------------------
  void CTime::setTimeUT(int year, int month, int day, double hourUT)
    throw(CTime::TimeEXC) {
   // algorithm from: Montenbruck: Grundlagen der Ephemeridenrechnung (6. Auflage), 2001

    double y;
    double m, B;

    // year 0 ?
    if (year == 0)
      throw TimeEXC("CTime::setTime(int, int, int, double): invalid date(year)");

    if (year < -4713.0)
      throw TimeEXC("CTime::setTime(int,int,int,double): year < -4713)");

    // 0<month<13 ?
    if (month < 1 || month > 12)
      throw TimeEXC("CTime::setTime(int,int,int,double):invalid date(month)");

    // 0<day<31 ?
    if (day < 1 || day > 31)
      throw TimeEXC("CTime::setTime(int,int,int,double):invalid date(day)");

    // 29th of february ?
    if (month == 2 && day > 29)
      throw TimeEXC("CTime::setTime(int,int,int,double):invalid date: day of february>29?");

    // hour<0 or >=24
    if (hourUT < 0.0 || hourUT >= 24.00)
      throw TimeEXC("CTime::setTime(int,int,int,double):invalid date(hour)");

    if (month <= 2) {
      y = year - 1;
      m = month + 12;
    }
    else {
      y = year;
      m = month;
    }

    if (year * 365 + month * 30 + day <= 577734)
      B = -2.0;
    else
      B = floor(year / 400.0) - floor(year / 100.0);

      m_value = floor(365.25 * y) +
                floor(30.6001 * (m + 1)) + B + 1720996.5 + day + hourUT / 24.0;

#ifdef SLTIME_DEBUG
    std::cout << "CTime::setTime(int year, int month,";
    std::cout << " int day, double hourUT)=>" << m_value << std::endl;
#endif

  }


// ----------------------------------------------------------------------------
  double CTime::JD(void) const {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::JD(void)=>" << m_value << std::endl;
#endif

    return m_value;
  }


// ----------------------------------------------------------------------------
  double CTime::MJD(void) const {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::MJD(void)=>" << (m_value - JD_MJD_0) << std::endl;
#endif

    return m_value - JD_MJD_0;
  }


// ----------------------------------------------------------------------------
  double CTime::JC(const CTime & since) const {

#ifdef SLTIME_DEBUG
    std::cout << "CTime::JC(const CTime &)=>"
              << (m_value - since.JD() * DAY2JYR) << std::endl;
#endif

    return (m_value - since.JD()) * DAY2JYR * JYR2JCY;
  }


// ----------------------------------------------------------------------------
  void CTime::dateUT(int &year, int &month, int &day, double &hourUT) const {

    JD2date(m_value, year, month, day, hourUT);

#ifdef SLTIME_DEBUG
    std::cout << "CTime::date(int&, int&, int&,double &) =>"
              << m_value << std::endl;
#endif
  }


// ----------------------------------------------------------------------------
  void CTime::dateTT(int &year, int &month, int &day, double &hourTT) const {

    double dT;

    JD2date(m_value, year, month, day, hourTT); // get UT from julian date

    try {
      dT = deltaT(m_value);     // try to calculate dT=TT-UT
    }
    catch(TimeEXC _exc) {           // catch exception, which indicates
      std::cout << _exc.what() << std::endl;             // whether dT could have been
      return;                   // calculated or not
    }

    hourTT += dT / 3600.0;      // transform UT -> TT

    if (hourTT >= 24.00) {      // overflow; hour >=24.00?
      double temp = m_value + dT * SEC2DAY; // temporary julian date with dT

      JD2date(temp, year, month, day, hourTT);  // calc new calendar date
    }
  }


// ----------------------------------------------------------------------------
  void CTime::dateTAI(int &year, int &month, int &day, double &hourTAI) const {

    double dT, temp;

    JD2date(m_value, year, month, day, hourTAI);  // get UT from julian date

    try {
      dT = deltaT(m_value);     // try to calculate dT=TT-UT
    }
    catch(TimeEXC _exc) {           // catch exception, which indicates
      std::cout << _exc.what() << std::endl;   // whether dT could be
      return;                   // calculated or not
    }

    hourTAI += ((dT - 32.184) / 3600.0);  // transform UT -> TAI

    // overflow; hour >=24.00 or hour < 00.00?
    if (hourTAI >= 24.00 || hourTAI < 0.0) {
      temp = m_value + dT * SEC2DAY - 32.84 * SEC2DAY;
      JD2date(temp, year, month, day, hourTAI); // calc new calendar date
    }
  }


// ----------------------------------------------------------------------------
  void CTime::dateUTC(int &year, int &month, int &day, double &hourUTC) const {

    double dT;              // TT-UT
    double dT2;             // TAI-UTC
    double temp;

    JD2date(m_value, year, month, day, hourUTC);  // get UT from julian date

    try {
      dT = deltaT(m_value);     // try to calculate dT=TT-UT
    }
    catch(TimeEXC & _exc) {         // catch exception, which indicates
      std::cout << _exc.what() << std::endl; // whether there dT could have been
      return;                   // calculated or not
    }

    try {
      dT2 = get_leapseconds(year, month, day) * 1.0;  // get leap seconds from table
    }
    catch(TimeEXC & _exc) {         // catch exception, whether the
      std::cout << _exc.what() << std::endl; // desired date is covered by our
      return;                   // table
    }

    hourUTC += ((dT - 32.184 - dT2) / 3600.0);  // transform UT -> TAI

    if (hourUTC >= 24.00 || hourUTC < 0.0) {  // overflow ?
      temp = m_value + (dT - 32.84 - dT2) * SEC2DAY;
      JD2date(temp, year, month, day, hourUTC); // calc new calendar date
    }
  }


// ----------------------------------------------------------------------------
  CTime & CTime::operator=(const CTime & right) {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::operator= =>" << right << std::endl;
#endif

    m_value = right.JD();
    m_format = right.getOutputFormat();
    return *this;
  }


// ----------------------------------------------------------------------------
  CTime CTime::operator+(const CTime & right) {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::operator+ =>" << right << std::endl;
#endif

    return CTime(JD() + right.JD());
  }


// ----------------------------------------------------------------------------
  CTime CTime::operator-(const CTime & right) {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::operator- =>" << right << std::endl;
#endif

    return CTime(JD() - right.JD());
  }


// ----------------------------------------------------------------------------
  CTime & CTime::operator+=(const CTime & right) {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::operator+= =>" << right << std::endl;
#endif

    m_value += right.JD();
    return *this;
  }


// ----------------------------------------------------------------------------
  CTime & CTime::operator-=(const CTime & right) {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::operator-= =>" << right << std::endl;
#endif

    m_value -= right.JD();
    return *this;
  }


// ----------------------------------------------------------------------------
  TimeUnit CTime::getOpUnit(void) const {
    return m_unit;
  }


// ----------------------------------------------------------------------------
  void CTime::setOpUnit(TimeUnit _u) {
    m_unit = _u;
  }


// ----------------------------------------------------------------------------
  CTime & CTime::operator++() {
    switch (m_unit) {
         case TimeUnit_second:
           m_value += SEC2DAY;
           break;
         case TimeUnit_minute:
           m_value += SEC2DAY * 60.0;
           break;
         case TimeUnit_hour:
           m_value += SEC2DAY * 60.0 * 60.0;
           break;
         case TimeUnit_day:
           m_value += 1.0;
           break;
         case TimeUnit_Jmonth:
           m_value += JYR2DAY / 12.0;
           break;
         case TimeUnit_Jyear:
           m_value += JYR2DAY;
           break;
         case TimeUnit_Jcentury:
           m_value += JCY2DAY;
           break;
         default:
           break;
    }
    return *this;
  }


// ----------------------------------------------------------------------------
  CTime & CTime::operator--() {
    switch (m_unit) {
         case TimeUnit_second:
           m_value -= SEC2DAY;
           break;
         case TimeUnit_minute:
           m_value -= SEC2DAY * 60.0;
           break;
         case TimeUnit_hour:
           m_value -= SEC2DAY * 60.0 * 60.0;
           break;
         case TimeUnit_day:
           m_value -= 1.0;
           break;
         case TimeUnit_Jmonth:
           m_value -= JYR2DAY / 12.0;
           break;
         case TimeUnit_Jyear:
           m_value -= JYR2DAY;
           break;
         case TimeUnit_Jcentury:
           m_value -= JCY2DAY;
           break;
         default:
           break;
    }
    return *this;
  }


// ----------------------------------------------------------------------------
  CTime CTime::operator++(int) {
    CTime t;

    t = *this;
    ++(*this);
    return t;
  }

  CTime CTime::operator--(int) {
    CTime t = *this;

    --(*this);
    return t;
  }


// ----------------------------------------------------------------------------
  TimeFormat CTime::getOutputFormat(void) const {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::getOutputFormat() =>" << std::endl;
#endif

    return m_format;
  }


// ----------------------------------------------------------------------------
  void CTime::setOutputFormat(TimeFormat _format) {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::setOutputFormat() =>" << _format << std::endl;
#endif

    m_format = _format;
  }


// ----------------------------------------------------------------------------
/** \brief The stream operator. */
  std::ostream & operator<<(std::ostream & h, const CTime & t) {
#ifdef SLTIME_DEBUG
    std::cout << "CTime::operator<< =>" << std::endl;
#endif

    int ly, lm, ld;
    double lh;
    double temp = 0.0;

    int _format;

    _format = t.getOutputFormat();

    int w = h.width();
    int p = h.precision();

    switch (_format) {
         case TimeFormat_JD:

           h << std::setw(14) << std::setprecision(6);
           h << std::setiosflags(std::ios::fixed);
           h << t.JD();
           h << std::resetiosflags(std::ios::fixed);

           break;
         case TimeFormat_MJD:

           h.width(15);
           h << std::setprecision(4);
           h << (t.MJD());

           break;
         case TimeFormat_HDMY_UT:
           t.dateUT(ly, lm, ld, lh);
           temp = lh;
           break;
         case TimeFormat_HDMY_TAI:
           t.dateTAI(ly, lm, ld, lh);
           temp = lh;
           break;
         case TimeFormat_HDMY_UTC:
           t.dateUTC(ly, lm, ld, lh);
           temp = lh;
           break;
         case TimeFormat_HDMY_TT:
           t.dateTT(ly, lm, ld, lh);
           temp = lh;
           break;
    }

    if (temp != 0.0) {
      h << std::setiosflags(std::ios::fixed);

      // date output
      h << std::setw(5) << std::setprecision(0) << ly << "/";
      h << std::setw(2) << std::setprecision(0) << lm << "/";
      h << std::setw(2) << std::setprecision(0) << ld << ", ";

      // calc hour, minute, second
      double hh = floor(temp);
      double mm = temp = (temp - floor(temp)) * 60.0;
      double ss = temp = (temp - floor(temp)) * 60.0;

      // time output
      h << std::setw(2) << std::setprecision(0) << hh << ":";
      h << std::setw(2) << std::setprecision(0) << mm << ":";
      h << std::setw(6) << std::setprecision(3) << ss;

      h << std::resetiosflags(std::ios::fixed);
    }

    h.width(w);
    h.precision(p);

    return h;
  }  // end of operator<<

}  // end of namespace SpaceLib

// end of file

