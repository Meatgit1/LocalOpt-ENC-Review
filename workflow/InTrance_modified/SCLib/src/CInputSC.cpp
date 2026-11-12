/*!
\file CInputSC.cpp
\brief An auxiliarary class that contains the values of the spacecraft input file.

\version

  1.7

\author

  Bernd Dachwald (BD)\n
  German Aerospace Center (DLR)\n
  Institute of Space Simulation\n
  Linder Hoehe\n
  51147 Cologne\n
  Germany\n
  bernd.dachwald@dlr.de


\author

  Andreas Ohndorf (AO)\n
  German Aerospace Center (DLR)\n
  Institute of Space Simulation\n
  Linder Hoehe\n
  51147 Cologne\n
  Germany\n
  andreas.ohndorf@gmx.de

\date

  29.10.2006 (AO)

\b Compiler:

  GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++

\b Copyright:

  (c) 2001-2006 by Bernd Dachwald, Andras Ohndorf

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/

//SCLib
#include <configSC.h>

// STL
#include <algorithm>
#include <sstream>

// SpaceLib
#include <SpaceLib.h>

// SCLib
#include <SCConst.h>
#include <CInputSC.h>

using namespace std;

//------------------------------------------------------------------------------
void CInputSC::init(void)
{
  psystype = SOLARSAIL;
  name = DFLT_SC_NAME;
  moments_of_inertia.Set(Vector(3));
  max_steering_moments.Set(Vector(3));
  sailcraft_performance = 0.0;
  ideal_or_real_sail = true;
  max_temperature = 0.0;
  reflection_coefficient_new = 0.0;
  specular_reflection_coefficient_new = 0.0;
  emission_coefficient_front_new = 0.0;
  emission_coefficient_back_new = 0.0;
  reflection_coefficient_aged = 0.0;
  specular_reflection_coefficient_aged = 0.0;
  emission_coefficient_front_aged = 0.0;
  emission_coefficient_back_aged = 0.0;
  reflection_coefficient_halflife = 0.0;
  specular_reflection_coefficient_halflife = 0.0;
  emission_coefficient_front_halflife = 0.0;
  emission_coefficient_back_halflife = 0.0;
  non_lambertian_coefficient_front = 0.0;
  non_lambertian_coefficient_back = 0.0;
  sc_mass = 0.0;
  min_mass_initial_propellant = 0.0;
  max_mass_initial_propellant = 0.0;
  throttletype = UNDEFINED_TT;
  thrust_max = 0.0;
  isp = 0.0;
  char_power = 0.0;
  p_exponent = 0.0;
  thrusters = 0.0;
} // CInputSC::init(void)

//------------------------------------------------------------------------------
CInputSC::CInputSC(ParamList& parameters, const string& trailing_chars)
{
  init();

  ParamList::iterator i;
  ParamList::iterator first(parameters.begin());
  ParamList::iterator last(parameters.end());

  i = find_if(first, last, name_equal<Param>("FILE_TYPE"));
  if (i!=last)
  {
    removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
    removeChars((*i).second, string(" \""));          // remove possible quotes
    upper((*i).second);                               // upcase whole string
    if ((*i).second.compare("INTRANCE1.7SPACECRAFTPARAMETERFILE")!=0)
      throw invalid_argument(string("CInputSC::CInputSC(): not an InTrance 1.7 spacecraft parameter file\n"));
  }
  else throw invalid_argument(string("CInputSC::CInputSC(): not an InTrance 1.7 spacecraft parameter file\n"));
  parameters.erase(i);first=parameters.begin();last=parameters.end();

  string sc_type;
  i = find_if(first, last, name_equal<Param>("SC_TYPE"));
  if (i!=last) 
  {
    removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
    removeChars((*i).second, string(" \""));          // remove possible quotes
    upper((*i).second);                               // upcase whole string
    sc_type=(*i).second;
  } else throw invalid_argument("CInputSC::CInputSC(): missing parameter (SC_TYPE)\n");
  parameters.erase(i);first=parameters.begin();last=parameters.end();

  if (sc_type == "SOLARSAIL") psystype = SOLARSAIL;
  else if (sc_type == "NEP")       psystype = NEP;
  else if (sc_type == "NSTAR")     psystype = NSTAR;
  else if (sc_type == "T6")        psystype = T6;
  else if (sc_type == "RIT")       psystype = RIT;
  else throw invalid_argument(string("CInputSC::CInputSC(): unknown spacecraft type\n"));

  if (psystype == SOLARSAIL)
  {
    i = find_if(first, last, name_equal<Param>("LIGHTNESS_NUMBER"));
    if (i!=last)
    {
      removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
      removeChars((*i).second, string(" \""));          // remove possible quotes
      istringstream((*i).second) >> sailcraft_performance;
      parameters.erase(i);first=parameters.begin();last=parameters.end();
    }
    else
    {
      i = find_if(first, last, name_equal<Param>("CHARACTERISTIC_ACCELERATION"));
      if (i!=last)
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> sailcraft_performance;
        sailcraft_performance/=5.93010834;
        parameters.erase(i);first=parameters.begin();last=parameters.end();
      }
      else throw invalid_argument(string("CInputSC::CInputSC(): missing parameter (LIGHTNESS_NUMBER or CHARACTERISTIC_ACCELERATION)\n"));
    }

    string ideal_or_real;
    i = find_if(first, last, name_equal<Param>("OPTICAL_SAIL_MODEL"));
    if (i!=last) 
    {
      removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
      removeChars((*i).second, string(" \""));          // remove possible quotes
      upper((*i).second);                               // upcase whole string
      ideal_or_real=(*i).second;
    }
    else throw invalid_argument(string("CInputSC::CInputSC(): missing parameter (OPTICAL_SAIL_MODEL)\n"));
    parameters.erase(i);first=parameters.begin();last=parameters.end();


    if (ideal_or_real == "IDEALSAIL")
    {
      ideal_or_real_sail = true;
    }

    else if (ideal_or_real == "REALSAIL")
    {
      ideal_or_real_sail = false;

      i = find_if(first, last, name_equal<Param>("MAX_SAIL_TEMP"));
      if (i!=last) {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> max_temperature;
      }
      else throw invalid_argument(string(" CInputSC::InputSC(): missing parameter (MAX_SAIL_TEMP)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("REFLECTION_COEFF_NEW"));
      if (i!=last) 
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> reflection_coefficient_new;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (REFLECTION_COEFF_NEW)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("SPECULAR_REFLECTION_COEFF_NEW"));
      if (i!=last)
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> specular_reflection_coefficient_new;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (SPECULAR_REFLECTION_COEFF_NEW)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("FRONT_EMISSION_COEFF_NEW"));
      if (i!=last) 
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> emission_coefficient_front_new;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (FRONT_EMISSION_COEFF_NEW)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("BACK_EMISSION_COEFF_NEW"));
      if (i!=last)
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> emission_coefficient_back_new;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (BACK_EMISSION_COEFF_NEW)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("REFLECTION_COEFF_AGED"));
      if (i!=last) 
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> reflection_coefficient_aged;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (REFLECTION_COEFF_AGED)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("SPECULAR_REFLECTION_COEFF_AGED"));
      if (i!=last) 
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> specular_reflection_coefficient_aged;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (SPECULAR_REFLECTION_AGED)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("FRONT_EMISSION_COEFF_AGED"));
      if (i!=last)
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> emission_coefficient_front_aged;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (FRONT_EMISSION_COEFF_AGED)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("BACK_EMISSION_COEFF_AGED"));
      if (i!=last) 
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> emission_coefficient_back_aged;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (BACK_EMISSION_COEFF_AGED)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("REFLECTION_COEFF_HALFLIFE"));
      if (i!=last) 
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> reflection_coefficient_halflife;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (REFLECTION_COEFF_HALFLIFE)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("SPECULAR_REFLECTION_COEFF_HALFLIFE"));
      if (i!=last)
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> specular_reflection_coefficient_halflife;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (SPECULAR_REFLECTION_COEFF_HALFLIFE)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("FRONT_EMISSION_COEFF_HALFLIFE"));
      if (i!=last) 
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> emission_coefficient_front_halflife;
      } 
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (FRONT_EMISSION_COEFF_HALFLIFE)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("BACK_EMISSION_COEFF_HALFLIFE"));
      if (i!=last) 
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> emission_coefficient_back_halflife;
      } 
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (BACK_EMISSION_COEFF_HALFLIFE)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("FRONT_NON_LAMBERT_COEFF"));
      if (i!=last)
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> non_lambertian_coefficient_front;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (FRONT_NON_LAMBERT_COEFF)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("BACK_NON_LAMBERT_COEFF"));
      if (i!=last) 
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> non_lambertian_coefficient_back;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (BACK_NON_LAMBERT_COEFF)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();

    } // end if (ideal_or_real == "REALSAIL")

    else throw invalid_argument(string("CInputSC::InputSC(): CInputSC::CInputSC(): unknown optical sail model\n"));

  } // end if (psystype == SOLARSAIL)

  else if (psystype == NEP || psystype == NSTAR || psystype == T6 || psystype == RIT)
  {
    i = find_if(first, last, name_equal<Param>("SC_DRYMASS"));
    if (i!=last)
    {
      removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
      removeChars((*i).second, string(" \""));          // remove possible quotes
      istringstream((*i).second) >> sc_mass;
    }
    else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (SC_DRYMASS)\n"));
    parameters.erase(i);first=parameters.begin();last=parameters.end();


    i = find_if(first, last, name_equal<Param>("MIN_PROP_MASS"));
    if (i!=last) 
    {
      removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
      removeChars((*i).second, string(" \""));          // remove possible quotes
      istringstream((*i).second) >> min_mass_initial_propellant;
    }
    else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (MIN_PROP_MASS)\n"));
    parameters.erase(i);first=parameters.begin();last=parameters.end();


    i = find_if(first, last, name_equal<Param>("MAX_PROP_MASS"));
    if (i!=last)
    {
      removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
      removeChars((*i).second, string(" \""));          // remove possible quotes
      istringstream((*i).second) >> max_mass_initial_propellant;
    }
    else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (MAX_PROP_MASS)\n"));
    parameters.erase(i);first=parameters.begin();last=parameters.end();


    string thr_type;
    i = find_if(first, last, name_equal<Param>("THROTTLE_TYPE"));
    if (i!=last)
    {
      removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
      removeChars((*i).second, string(" \""));          // remove possible quotes
      upper((*i).second);
      thr_type=(*i).second;
    }
    else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (THROTTLE_TYPE)\n"));
    parameters.erase(i);first=parameters.begin();last=parameters.end();


    if (thr_type == "BANG-BANG")
      throttletype = BANGBANG;
    else if (thr_type == "VARIABLE")
      throttletype = VARIABLE;
    else throw invalid_argument(string("CInputSC::InputSC():  unknown throttle type\n"));


    if (psystype == NEP)
    {
      i = find_if(first, last, name_equal<Param>("MAX_THRUST"));
      if (i!=last) 
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> thrust_max;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (MAX_THRUST)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("SPECIFIC_IMPULSE"));
      if (i!=last) 
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> isp;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (SPECIFIC_IMPULSE)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();
    }


    else // if (psystype == NSTAR || psystype == T6 || psystype == RIT)
    {
      i = find_if(first, last, name_equal<Param>("SOLAR_ARRAY_CHAR_POWER"));
      if (i!=last) 
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> char_power;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (SOLAR_ARRAY_CHAR_POWER)\n"));
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("POWER_VARIATION_EXPONENT"));
      if (i!=last) 
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> p_exponent;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (POWER_VARIATION_EXPONENT)\n")); 
      parameters.erase(i);first=parameters.begin();last=parameters.end();


      i = find_if(first, last, name_equal<Param>("N_THRUSTERS"));
      if (i!=last)
      {
        removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
        removeChars((*i).second, string(" \""));          // remove possible quotes
        istringstream((*i).second) >> thrusters;
      }
      else throw invalid_argument(string("CInputSC::InputSC(): missing parameter (N_THRUSTERS)\n")); 
      parameters.erase(i);first=parameters.begin();last=parameters.end();

    } // end if (psystype == NEP)

  } // end if (psystype == NEP || psystype == NSTAR || psystype == T6 || psystype == RIT)

  // parse possibly given moments of inertia
  try
  {
    InputParameter("MOMENT_OF_INERTIA_1", moments_of_inertia[0], parameters,"", trailing_chars);
    InputParameter("MOMENT_OF_INERTIA_2", moments_of_inertia[1], parameters,"", trailing_chars);
    InputParameter("MOMENT_OF_INERTIA_3", moments_of_inertia[2], parameters,"", trailing_chars);
  }
  catch (invalid_argument&) {};

  // parse possibly given maximum steering moments
  try
  {
    InputParameter("MAX_STEERING_MOMENT_1", max_steering_moments[0], parameters,"", trailing_chars);
    InputParameter("MAX_STEERING_MOMENT_2", max_steering_moments[1], parameters,"", trailing_chars);
    InputParameter("MAX_STEERING_MOMENT_3", max_steering_moments[2], parameters,"", trailing_chars);
  } catch (std::invalid_argument&) {}

  try 
  {
    string s_temp("");
    InputStringParameter("SC_NAME", s_temp, parameters, "", trailing_chars);
    if (!s_temp.empty())
    {
      name = s_temp;
      upper(name);
    } // if (!s_temp.empty())
  } // try
  catch (...) {}
};

//------------------------------------------------------------------------------
void CInputSC::printSettings(void) const
{
  const unsigned int w(37);   // output field width
  ostringstream ss;

  // name
  ss << left << setw(w) << "SC: name" << ": ";
  ss << name << endl;Output(ss.str());

  // characteristic power
  ss.str("");
  ss << left << setw(w) << "SC: char. power"<<": ";
  ss << char_power << endl;Output(ss.str());

  // specific impulse
  ss.str("");
  ss << left << setw(w) << "SC: specific impulse [s]"<<": ";
  ss << isp << endl;Output(ss.str());

  // sc mass
  ss.str("");
  ss << left << setw(w) << "SC: spacecraft mass [kg]"<<": ";
  ss << sc_mass << endl;Output(ss.str());

  // max. initial propellant
  ss.str("");
  ss << left << setw(w) << "SC: max. mass initial prop. [kg]"<<": ";
  ss << max_mass_initial_propellant << endl;Output(ss.str());  

  // min. initial propellant
  ss.str("");
  ss << left << setw(w) << "SC: min. mass initial prop. [kg]"<<": ";
  ss << min_mass_initial_propellant << endl;Output(ss.str());  

  // power exponent
  ss.str("");
  ss << left << setw(w) << "SC: power exponent"<<": ";
  ss << p_exponent << endl;Output(ss.str());  

  // propulsion system type
  if (psystype != SOLARSAIL)
  {
    ss.str("");
    ss << left << setw(w) << "SC: propulsion system type"<<": ";
  }
  switch (psystype)
  {
  case NEP:ss<<"nuclear electric powered";break;
  case NSTAR:ss<<"NSTAR";break;
  case T6:ss<<"T6";break;
  case RIT:ss<<"RIT";break;
  default:ss << "undefined";
  }
  if (psystype != SOLARSAIL)
  {
    ss << endl;Output(ss.str());  
  }

  switch (psystype)
  {
  case SOLARSAIL: printSolarSailSettings();
  default:break;
  }


  // throttle type
  ss.str("");
  ss << left << setw(w) << "SC: throttle type"<<": ";
  switch (throttletype)
  {
  case UNDEFINED_TT:ss<<"undefined";break;
  case BANGBANG:ss<<"bang bang";break;
  case VARIABLE:ss<<"variable";break;
  default:ss << "undefined";
  }
  ss << endl;Output(ss.str());

  // max. thrust
  ss.str("");
  ss << left << setw(w) << "SC: max. thrust"<<": ";
  ss << thrust_max << endl;Output(ss.str());

  // # of thrusters
  ss.str("");
  ss << left << setw(w) << "SC: number of thrusters"<<": ";
  ss << thrusters << endl;Output(ss.str());

  // moments of inertia
  ss.str("");
  ss << left << setw(w) << "SC: m.o.i (x,y,z) [kg*m^2]"<<": ";
  ss << moments_of_inertia[0] << ", "; 
  ss << moments_of_inertia[1] << ", ";
  ss << moments_of_inertia[2];
  ss << endl;Output(ss.str());

  // max. steering moments
  ss.str("");
  ss << left << setw(w) << "SC: max. steer. moments (x,y,z) [Nm]"<<": ";
  ss << max_steering_moments[0] << ", ";
  ss << max_steering_moments[1] << ", "; 
  ss << max_steering_moments[2];
  ss << endl;Output(ss.str());
}

//------------------------------------------------------------------------------
void CInputSC::printSolarSailSettings(void) const
{
  const unsigned int w(37);   // output field width
  ostringstream ss;

  // propulsion system type
  ss.str("");
  ss << left << setw(w) << "SC: propulsion system type"<<": solar sail";
  ss << endl; Output(ss.str());

  // solar sail type
  ss.str("");
  ss << left << setw(w) << "SC: solar sail type"<<":";
  ss << (ideal_or_real_sail==true?" ideal sail":" real sail") << endl; Output(ss.str());

  // lightness number
  ss.str("");
  ss << left << setw(w) << "SC: lightness number"<<": "
     << sailcraft_performance << endl; 
  Output(ss.str());

  if (ideal_or_real_sail == false) // real sail (false)
  {
    // emission coefficient backside aged
    ss.str("");
    ss << left << setw(w) << "SC: backside emiss. coeff(aged)"<<": ";
    ss << emission_coefficient_back_aged << endl;Output(ss.str());

    // emission coefficient backside half life
    ss.str("");
    ss << left << setw(w) << "SC: backside emiss. coeff(half)"<<": ";
    ss << emission_coefficient_back_halflife << endl;Output(ss.str());

    // emission coefficient backside new
    ss.str("");
    ss << left << setw(w) << "SC: backside emiss. coeff(new)"<<": ";
    ss << emission_coefficient_back_new << endl;Output(ss.str());

    // emission coefficient frontside aged
    ss.str("");
    ss << left << setw(w) << "SC: frontside emiss. coeff(aged)"<<": ";
    ss << emission_coefficient_front_aged << endl;Output(ss.str());

    // emission coefficient frontside half life
    ss.str("");
    ss << left << setw(w) << "SC: frontside emiss. coeff(half)"<<": ";
    ss << emission_coefficient_front_halflife << endl;Output(ss.str());

    // emission coefficient frontside new
    ss.str("");
    ss << left << setw(w) << "SC: frontside emiss. coeff(new)"<<": ";
    ss << emission_coefficient_front_new << endl;Output(ss.str());

    // max. temperature of sail
    ss.str("");
    ss << left << setw(w) << "SC: max. sail temperature [K]"<<": ";
    ss << max_temperature << endl;Output(ss.str());  

    // non lambertian coefficient of sail backside
    ss.str("");
    ss << left << setw(w) << "SC: non. lambert. coeff. backs."<<": ";
    ss << non_lambertian_coefficient_back << endl;Output(ss.str());  

    // non lambertian coefficient of sail frontside
    ss.str("");
    ss << left << setw(w) << "SC: non. lambert. coeff. fronts."<<": ";
    ss << non_lambertian_coefficient_front << endl;Output(ss.str());

    // reflection coefficient aged
    ss.str("");
    ss << left << setw(w) << "SC: reflec. coeff. aged"<<": ";
    ss << reflection_coefficient_aged << endl;Output(ss.str());  

    // reflection coefficient halflife
    ss.str("");
    ss << left << setw(w) << "SC: reflec. coeff. halflife"<<": ";
    ss << reflection_coefficient_halflife << endl;Output(ss.str());

    // reflection coefficient new
    ss.str("");
    ss << left << setw(w) << "SC: reflec. coeff. new"<<": ";
    ss << reflection_coefficient_new << endl;Output(ss.str());

  } // if (ideal_or_eal_sail == false)


}

