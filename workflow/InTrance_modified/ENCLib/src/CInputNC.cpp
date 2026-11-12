/*!
\file CInputNC.cpp
\brief An auxiliarary class that contains the values of the neurocontroller parameter file.

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
  andreas.ohndorf@dlr.de


\date

  01.12.2005 (AO)

\b Compiler:

  GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++

\b Copyright:

  (c) 2001-2004 by Bernd Dachwald, Andreas Ohndorf

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/

#if defined _MSC_VER             // Microsoft C++ Compiler
 #pragma warning (disable: 4786) // "identifier was truncated to '255' characters in the browser information"
 #pragma warning (disable: 4250) // "inherits ??? via dominance"
#endif

// STL
#include <algorithm>
#include <iomanip>
#include <sstream>

// SpaceLib
#include <Global.h>

// ENCLib
#include <CInputNC.h>


using namespace std;
using namespace Global;

//-----------------------------------------------------------------------------
CInputNC::CInputNC(void)
 : controller_type(DEFAULT_CONTROLLER_TYPE)
  ,input_type(IT_UNDEF)
  ,input_neurons(0)
  ,hidden_layers(0)
  ,hidden_neurons(0)
  ,output_neurons(0)
  ,activation_function(AF_UNDEF)
  ,extra_parameters(0)
{}

//-----------------------------------------------------------------------------
const CInputNC& CInputNC::operator =(const CInputNC& rhs)
{
  controller_type = rhs.controller_type;

  input_type          = rhs.input_type;
  input_neurons       = rhs.input_neurons;
  hidden_layers       = rhs.hidden_layers;

  SafeDeleteArray(hidden_neurons);
  if (hidden_layers > 0) 
  {
    hidden_neurons = new int[hidden_layers];
    if (hidden_neurons != 0) 
    {
      for (uint i = 0; i < hidden_layers; ++i) 
        hidden_neurons[i] = rhs.hidden_neurons[i];
    } else Global::FatalError("NCInput::operator=(const CInputNC&) => \
Error during allocation memory for hidden neuron information!");
  }
  output_neurons      = rhs.output_neurons;
  activation_function = rhs.activation_function;
  extra_parameters    = rhs.extra_parameters;
  return *this;
} // const CInputNC& CInputNC::operator =(const CInputNC& rhs)


//-----------------------------------------------------------------------------
CInputNC::CInputNC(
   ParamList& parameters
  ,bool extra_data/*==true*/
  ,const string& trailing_chars
)
: controller_type(DEFAULT_CONTROLLER_TYPE)
 ,input_type(IT_UNDEF)
 ,input_neurons(0)
 ,hidden_layers(0)
 ,hidden_neurons(NULL)
 ,output_neurons(0)
 ,activation_function(AF_UNDEF)
 ,extra_parameters(0)
{
  //check for file type and version
  string file_type;
  InputStringParameter("FILE_TYPE", file_type, parameters, 
    "CInputNC:CInputNC(ParamList&, bool, const std::string&)", trailing_chars);     
  upper(file_type);removeChars(file_type, " \"");
  if (file_type.compare("INTRANCE1.7NEUROCONTROLLERPARAMETERFILE")!=0)
    throw invalid_argument(string("CInputNC::CInputNC(): not an InTrance 1.7 neurocontroller parameter file\n"));


  // controller type (application specific and therefor must not be given
  try
  {
    InputParameter("CTRL_TYPE", controller_type, parameters, "", trailing_chars);
  }
  catch (std::invalid_argument&) {};


  // spacecraft type
  string sctype;
  if (extra_data == true)
  {
    InputStringParameter("SC_TYPE", sctype, parameters, 
      "CInputNC:CInputNC(ParamList&, bool, const std::string&)", trailing_chars);     
    upper(sctype);removeChars(sctype, " \"");
  }


  // read, if input and output neurons for gravity assist optimization are to be included
  string include_ga;
  if (extra_data == true)
  {
    InputStringParameter("INCLUDE_GA_NEURONS", include_ga, parameters, 
      "CInputNC:CInputNC(ParamList&, bool, const std::string&)", trailing_chars);     
    upper(include_ga);removeChars(include_ga, " \"");
  }


  // read input type
  string cart_yes_no;
  if (extra_data == true)
  {
    InputStringParameter("CARTESIAN_INPUT", cart_yes_no, parameters, 
      "CInputNC:CInputNC(ParamList&, bool, const std::string&)", trailing_chars);     
    upper(cart_yes_no);removeChars(cart_yes_no, " \"");

    if (cart_yes_no == "YES") {
      if (include_ga == "YES") input_neurons += 18;
      else input_neurons += 12;
      input_type = C_INPUT;
    }
    else if (cart_yes_no == "NO")
     ;
    else throw invalid_argument(string("CInputNC::CInputNC(): error in NC parameter file\n"));
  }



  // spherical coordinate input
  string sph_yes_no;
  if (extra_data == true)
  {
    InputStringParameter("SPHERICAL_INPUT", sph_yes_no, parameters, 
      "CInputNC:CInputNC(ParamList&, bool, const std::string&)", trailing_chars);     
    upper(sph_yes_no);removeChars(sph_yes_no, " \"");

    if (sph_yes_no == "YES") {
      if (include_ga == "YES") input_neurons += 18;
      else input_neurons += 12;
      if (input_type == C_INPUT) input_type = CP_INPUT;
      else // input_type == IT_UNDEF
        input_type = P_INPUT;
    }
    else if (sph_yes_no == "NO")
      ;
    else throw invalid_argument(string("CInputNC::CInputNC(): error in NC parameter file\n"));
  } // if (extra_data == true)



  // Kepler elements input
  string elem_yes_no;
  if (extra_data == true)
  {
    InputStringParameter("ORBITAL_ELEMENTS_INPUT", elem_yes_no, parameters, 
      "CInputNC:CInputNC(ParamList&, bool, const std::string&)", trailing_chars);     
    upper(elem_yes_no);removeChars(elem_yes_no, " \"");

    if (elem_yes_no == "YES") {
      if (include_ga == "YES") input_neurons += 18;
      else input_neurons += 12;
      if (input_type == C_INPUT) input_type = CE_INPUT;
      else if (input_type == P_INPUT) input_type = PE_INPUT;
      else if (input_type == CP_INPUT) input_type = CPE_INPUT;
      else // input_type == IT_UNDEF
        input_type = E_INPUT;
    }
    else if (elem_yes_no == "NO")
      ;
    else throw invalid_argument(string("CInputNC::CInputNC(): error in NC \
parameter file\n"));
  } //if (extra_data == true)



  // delta coding or full values
  string full_delta;
  if (extra_data == true)
  {
    InputStringParameter("NC_INPUT", full_delta, parameters, 
      "CInputNC:CInputNC(ParamList&, bool, const std::string&)", trailing_chars);     
    upper(full_delta);removeChars(full_delta, " \"");

    if (full_delta == "FULL")
     ;
    else if (full_delta == "DELTA") {
      if (include_ga == "YES") input_neurons = input_neurons * 2 / 3;
      else input_neurons /= 2;
    }
    else throw invalid_argument(string("CInputNC::CInputNC(): error in NC \
parameter file\n"));
  } // if (extra_data == true)


  if (extra_data == true)
  {
    if (sctype == "SOLARSAIL")
      ;
    else if (sctype == "NEP" || sctype == "SEP") {
      input_neurons += 1;
      output_neurons += 2;
    }
    else throw invalid_argument((string("CInputNC::InputNC(): unknown spacecraft type ") + sctype + string("\n")).c_str());
  } // if (extra_data == true)


  // read number of hidden layers
  InputParameter("HIDDEN_LAYERS", hidden_layers,
    parameters,"CInputNC::CInputNC(ParamList&, bool, const std::string&)",
    trailing_chars);

  // read numbers of hidden layer neurons
  if (hidden_layers) {
    hidden_neurons = new int[hidden_layers];
    Assure(hidden_neurons, "CInputNC::CInputNC()");
    for (uint idx = 0; idx < hidden_layers; idx++)
    {
      ostringstream t("");
      t << (idx+1);             // set new index of hidden layer (string)

      // and try to read number of neurons in current hidden layer
      InputParameter("NEURONS_IN_HIDDEN_LAYER"+t.str(), hidden_neurons[idx],
        parameters,"CInputNC::CInputNC(ParamList&, bool, const std::string&)",
        trailing_chars);
    }
  } else hidden_neurons = NULL;


  // # of input and output neurons
  if (extra_data == false)
  {
     InputParameter("INPUT_NEURONS", input_neurons, parameters, 
       "CInputNC:CInputNC(ParamList&, bool, const std::string&)", trailing_chars);     

     InputParameter("OUTPUT_NEURONS", output_neurons, parameters, 
       "CInputNC:CInputNC(ParamList&, bool, const std::string&)", trailing_chars);
  } // if (extra_data = false)

  // read output type
  string output_type;
  if (extra_data == true) 
  {
    InputStringParameter("NC_OUTPUT", output_type, parameters, 
      "CInputNC:CInputNC(ParamList&, bool, const std::string&)", trailing_chars);     
    upper(output_type);removeChars(output_type, " \"");

    if (output_type == "LSLMINIMIZEMAXIMIZE") {
      output_neurons += 10;
    } else if (output_type == "LSLOPTIMIZE") {
      output_neurons += 5;
    } else if (output_type == "DIRECT") {
      output_neurons += 3;
    } else throw invalid_argument(string("CInputNC::CInputNC(): unknown output type ") + output_type+string("\n"));
  } // if (extra_data == true)

  // possible gravity assist neurons
  if (extra_data == true && include_ga == "YES") output_neurons += 2;

  // read neural activation function
  string af;
  InputStringParameter("TRANSFER_FUNCTION", af, parameters, 
    "CInputNC:CInputNC(ParamList&, bool, const std::string&)", trailing_chars);     
  upper(af);removeChars(af, " \"");

  if (af == "SIGMOID") activation_function = SIGMOID;
  else if (af == "TANH") activation_function = TANH;
  else if (af == "SIMPLE") activation_function = SIMPLE;
  else if (af == "SIGMOID,IDENTITYFOROUTPUTNEURONS") activation_function = SIGMOIDIDENTITY;
  else if (af == "TANH,IDENTITYFOROUTPUTNEURONS") activation_function = TANHIDENTITY;
  else if (af == "SIMPLE,IDENTITYFOROUTPUTNEURONS") activation_function = SIMPLEIDENTITY;
  else throw invalid_argument(string("CInputNC::InputNC(): unknown neural activation function ") + af+string("\n"));
};


//-----------------------------------------------------------------------------
void CInputNC::printSettings(void) const
{
  const unsigned int w(37);   // output field width
  ostringstream ss;

  // controller type (application specific)
  ss << left << setw(w) << "NC: controller type"<<": ";
  ss << controller_type << endl;Output(ss.str());

  // activation function
  ss.str("");ss << left << setw(w) << "NC: activation_function"<< ": ";
  switch (activation_function) 
  {
    case SIGMOID: ss << "SIGMOID";break;
    case TANH:    ss << "TANH";break;
    case SIMPLE:  ss << "SIMPLE";break;
    case SIGMOIDIDENTITY: ss << "SIGMOID,IDENTITYFOROUTPUTNEURONS";break;
    case TANHIDENTITY: ss << "TANH,IDENTITYFOROUTPUTNEURONS"; break;
    case SIMPLEIDENTITY: ss << "SIMPLE,IDENTITYFOROUTPUTNEURONS";break;
    default: ss << " unknown";break;
  }
  ss<<endl;Output(ss.str());

  // extra parameters
  ss.str("");ss << left << setw(w) << "NC: extra parameters"<< ": ";
  ss << extra_parameters << endl;
  Output(ss.str());

  // # of hidden layers
  ss.str("");ss << setw(w) << "NC: hidden layers" << ": " << hidden_layers << endl;
  Output(ss.str());

  // # of hidden neurons
  if (hidden_layers>0)
  {
    ss.str("");ss << setw(w) << "NC: hidden neurons" << ": (";
    for (unsigned int i=0;i<hidden_layers-1;++i) ss << hidden_neurons[i] <<",";
    ss << hidden_neurons[hidden_layers-1] << ")\n";
    Output(ss.str());
  }

  // # of input neurons
  ss.str("");ss << setw(w) << "NC: input neurons" << ": " << input_neurons << endl;
  Output(ss.str());

  // # of output neurons
  ss.str("");
  ss << setw(w) << "NC: output neurons" << ": " << output_neurons << endl;
  Output(ss.str());

  // input type
  ss.str("");ss << setw(w) << "NC: input type" << ": ";
  switch (input_type)
  {
    case IT_UNDEF: ss << "undefined";break;
    case C_INPUT:  ss << "cartesian coordinates";break;
    case P_INPUT:  ss << "polar coordinates";break;
    case E_INPUT:  ss << "orbital elements";break;
    case CP_INPUT: ss << "cartesian + polar coordinates";break;
    case CE_INPUT: ss << "cartesian coordinates + orbital elements";break;
    case PE_INPUT: ss << "polar coordinates + orbital elements";break;
    case CPE_INPUT:ss << "cartesian + polar coordinates + orbital elements";break;
  }
  ss << endl;
  Output(ss.str());
}


//-----------------------------------------------------------------------------
CInputNC::CInputNC( const NCInputType          p_input_type
                   ,const int                  p_input_neurons
                   ,const int                  p_hidden_layers
                   ,      int*                 p_hidden_neurons
                   ,const int                  p_output_neurons
                   ,const NCActivationFunction p_AcFunc
                   ,const int                  p_extra_parameters) 
 : controller_type(DEFAULT_CONTROLLER_TYPE)
{

  input_type = p_input_type;
  input_neurons = p_input_neurons;
  hidden_layers = p_hidden_layers;
  hidden_neurons = p_hidden_neurons;
  output_neurons = p_output_neurons;
  activation_function = p_AcFunc;
  extra_parameters = p_extra_parameters;

}

//-----------------------------------------------------------------------------
void CInputNC::reset(void)
{
  controller_type = DEFAULT_CONTROLLER_TYPE;
  input_type = IT_UNDEF;
  input_neurons = 0;
  hidden_layers = 0;
  SafeDeleteArray(hidden_neurons);
  output_neurons = 0;
  activation_function = AF_UNDEF;
  extra_parameters = 0;
} // CInputNC::reset(void)
