/*!
\file CNC.cpp
\brief A class that provides the functionality for a neurocontroller.

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

  03.06.2005 (AO)

\b Compiler:

  GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++

\b Copyright:

  (c) 2001-2005 by Bernd Dachwald, Andreas Ohndorf

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/
#include <cmath>
#include <cassert>

#include <Global.h>

#include <CInputNC.h>
#include <CNC.h>


using namespace std;
using namespace Global;
using namespace MathLib;


//*****************************************************************************
//
// C N C N e u r o n
//
//*****************************************************************************

//------------------------------------------------------------------------------
CNCNeuron::CNCNeuron(const NCLayerType ltype, CNCLayer* connected_layer, CNC* nc)
: m_number_of_links(0),
  m_weight(NULL),
  m_bias(0.0),
  m_temperature_exponent(0.0),
  m_input(0.0),
  m_output(0.0),
  m_layer_type(ltype),
  m_connected_layer(connected_layer),
  m_nc(nc)
{
  if (connected_layer)
  {
    nc->incNeuronNumber();

    m_number_of_links = connected_layer->getNumberOfNeurons();

    m_weight = new double[m_number_of_links];
    Assure(m_weight, "CNCNeuron::CNCNeuron()");

    for (int i = 0; i < m_number_of_links; ++i)
    {
      m_weight[i] = 0.0;
      /*nc->m_parameter[nc->m_parameter_index] = &m_weight[i];*/
      nc->setParameter(nc->getParameterIndex(), &m_weight[i]);
      /*nc->m_parameter_neuron_number[nc->m_parameter_index] = nc->m_neuron_number;*/
      nc->setParameterNeuronNumber();
      /*++nc->m_parameter_index;*/
      nc->incParameterIndex();
    }

    /*nc->m_parameter[nc->m_parameter_index] = &m_bias;*/
    nc->setParameter(nc->getParameterIndex(), &m_bias);
    /*nc->m_parameter_neuron_number[nc->m_parameter_index] = nc->m_neuron_number;*/
    nc->setParameterNeuronNumber();
    /*++nc->m_parameter_index;*/
    nc->incParameterIndex();

    /*nc->m_parameter[nc->m_parameter_index] = &m_temperature_exponent;*/
    nc->setParameter(nc->getParameterIndex(), &m_temperature_exponent);
    /*nc->m_parameter_neuron_number[nc->m_parameter_index] = nc->m_neuron_number;*/
    nc->setParameterNeuronNumber();
    /*++nc->m_parameter_index;*/
    nc->incParameterIndex();
  }
};


//------------------------------------------------------------------------------
//CNCNeuron::~CNCNeuron() {SafeDeleteArray(m_weight);}
//{
//  if (m_weight)
//  {
//    delete []m_weight;
//    m_weight = NULL;
//  }
//};


//------------------------------------------------------------------------------
double CNCNeuron::getWeight(const int index) const
{
#ifdef development_version
  if(index >= m_number_of_links) FatalError("CNCNeuron::getWeight(): index out of range");
  if(m_weight == NULL) FatalError("CNCNeuron::getWeight(): input layer neurons do not have a weight");
#endif

  return m_weight[index];
};


//------------------------------------------------------------------------------
void CNCNeuron::setWeight(const int index, const double value)
{
#ifdef development_version
  if(index >= m_number_of_links) FatalError("CNCNeuron::getWeight(): index out of range");
  if(m_weight == NULL) FatalError("CNCNeuron::getWeight(): weight can not be set for input layer neurons");
#endif

  m_weight[index] = value;
};


//------------------------------------------------------------------------------
void CNCNeuron::setOutput(const double output)
{
#ifdef development_version
  if(m_weight != NULL) FatalError("CNCNeuron::setOutput(): output can only be set for input layer neurons");
#endif

  m_output = output;
};


//------------------------------------------------------------------------------
void CNCNeuron::Propagate(void)
{
  // sum over all inputs

  m_input = 0.0;

  for (int i = 0; i < m_number_of_links; ++i)
    m_input += m_weight[i] * m_connected_layer->Neuron(i)->getOutput();

  m_input += m_bias;
  m_input *= m_temperature_exponent;

  // calculate output from input using the activation function

  switch(m_nc->getActivationFunction())
  {
    case SIGMOID:
      m_output = 1.0/(1.0+exp(-m_input));
      break;

    case TANH:
            m_output = tanh(m_input);
      break;

    case SIMPLE:
      m_output = (m_input < 0.0) ? ((m_input < -2.0) ? -1.0 : -1.0 + 0.25 * (m_input+2.0)*(m_input+2.0)) :
                                   ((m_input >  2.0) ?  1.0 :  1.0 - 0.25 * (m_input-2.0)*(m_input-2.0));
      break;

    case SIGMOIDIDENTITY:
      m_output = m_input;
      if(getLayerType() == HIDDEN_LAYER)
        m_output = 1.0/(1.0+exp(-m_input));
      break;

    case TANHIDENTITY:
      m_output = m_input;
      if(getLayerType() == HIDDEN_LAYER)
        m_output = tanh(m_input);
      break;

    case SIMPLEIDENTITY:
      m_output = m_input;
      if(getLayerType() == HIDDEN_LAYER)
        m_output = (m_input < 0.0) ? ((m_input < -2.0) ? -1.0 : -1.0 + 0.25 * (m_input+2.0)*(m_input+2.0)) :
                                     ((m_input >  2.0) ?  1.0 :  1.0 - 0.25 * (m_input-2.0)*(m_input-2.0));
      break;

    default:
      FatalError("CNCNeuron::Propagate(): unknown activation function");
      break;
  }
};


//------------------------------------------------------------------------------
CNCLayer::CNCLayer(int index
                  ,int number_of_neurons
                  ,const NCLayerType ltype
                  ,CNCLayer* prev
                  ,CNC* nc)
: m_index(index)
 ,m_number_of_neurons(number_of_neurons)
{

  m_neuron = new CNCNeuron*[number_of_neurons];
  Assure(m_neuron, "CNCLayer::CNCLayer()");

  for (int i = 0; i < number_of_neurons; ++i) 
  {
    m_neuron[i] = new CNCNeuron(ltype, prev, nc);
    Assure(m_neuron[i], "CNCLayer::CNCLayer()");
  }

  if(prev) {   // concatenate layer list
    prev->m_next_layer = this;
  }
  m_prev_layer = prev;
  m_next_layer = NULL;
};


//------------------------------------------------------------------------------
CNCLayer::~CNCLayer() 
{
  for (int i = 0; i < m_number_of_neurons; ++i) SafeDelete(m_neuron[i]);
  SafeDeleteArray(m_neuron);
};


//-----------------------------------------------------------------------------
CNC::CNC(const CInputNC& input_ncp)
: m_nc_type(input_ncp.controller_type),
  m_input_type(input_ncp.input_type),
  m_activation_function(input_ncp.activation_function),
  m_parameter(NULL),
  m_parameter_neuron_number(NULL),
  m_parameter_index(0),
  m_neuron_number(0),
  m_number_of_extra_parameters(input_ncp.extra_parameters)
{
  m_number_of_neurons = m_number_of_links = 0;

  // read number of neurons in the input layer

  m_number_of_input_neurons = input_ncp.input_neurons;
  m_number_of_neurons += m_number_of_input_neurons;

  int number_of_neurons_in_previous_layer = m_number_of_input_neurons;

  // read number of hidden layers
  m_number_of_hidden_layers = input_ncp.hidden_layers;

  // read number of neurons in the hidden layers
  if (m_number_of_hidden_layers) 
  {
    m_number_of_hidden_neurons = new int[m_number_of_hidden_layers];
    Assure(m_number_of_hidden_neurons, "CNC::CNC()");
    for (int i = 0; i < m_number_of_hidden_layers; ++i) 
    {
      m_number_of_hidden_neurons[i] = input_ncp.hidden_neurons[i];
      m_number_of_neurons += m_number_of_hidden_neurons[i];
      m_number_of_links += number_of_neurons_in_previous_layer * m_number_of_hidden_neurons[i];
      number_of_neurons_in_previous_layer = m_number_of_hidden_neurons[i];
    }
  }
  else 
  {
    m_number_of_hidden_neurons = NULL;
  }

  // read number of output layer neurons
  m_number_of_output_neurons = input_ncp.output_neurons;
  m_number_of_neurons += m_number_of_output_neurons;
  m_number_of_links += number_of_neurons_in_previous_layer * m_number_of_output_neurons;

  m_number_of_parameters = 2 * (m_number_of_neurons - m_number_of_input_neurons) + m_number_of_links;

  // create parameter array etc.
  m_parameter = new double*[m_number_of_parameters];
  m_parameter_neuron_number = new int[m_number_of_parameters];

  // create layers, begin with input layer
  m_number_of_layers = 0;
  m_input_layer = new CNCLayer(m_number_of_layers++, m_number_of_input_neurons, INPUT_LAYER, NULL, this);
  Assure(m_input_layer, "CNC::CNC()");

  // create hidden layers
  CNCLayer *prev_layer, *act_layer;
  prev_layer = m_input_layer;

  for (int j = 0; j < m_number_of_hidden_layers; ++j) 
  {
    act_layer = new CNCLayer(m_number_of_layers++, m_number_of_hidden_neurons[j], HIDDEN_LAYER, prev_layer, this);
    Assure(act_layer, "CNC::CNC()");
    prev_layer = act_layer;
  }

  m_first_hidden_layer = m_input_layer->getNextLayer();

  // create output layer
  m_output_layer = new CNCLayer(m_number_of_layers++, m_number_of_output_neurons, OUTPUT_LAYER, prev_layer, this);
  Assure(m_output_layer, "CNC::CNC()");

  if(m_parameter_index != m_number_of_parameters) throw runtime_error("CNC::CNC(const CInputNC&):\n");
}; // CNC:CNC(const CInputNC& input_ncp)


//-----------------------------------------------------------------------------
CNC::~CNC(void)
{
  CNCLayer *to_delete_layer
          ,*next_to_delete_layer;

  to_delete_layer = m_output_layer;

  do {
    next_to_delete_layer = to_delete_layer->getPrevLayer();
    delete to_delete_layer;
    to_delete_layer = next_to_delete_layer;
  }
  while(to_delete_layer);

  SafeDeleteArray(m_number_of_hidden_neurons);
  SafeDeleteArray(m_parameter);
  SafeDeleteArray(m_parameter_neuron_number);
};


//------------------------------------------------------------------------------
void CNC::propagate(const Vector& input_vector, Vector& output_vector)
{
#ifdef development_version
  // check input and output vector dimensions
  assert(input_vector.getDimension() == m_input_layer->getNumberOfNeurons());
  assert(output_vector.getDimension() == m_output_layer->getNumberOfNeurons());
#endif

  int input_dim = input_vector.getDimension();
  int output_dim = output_vector.getDimension();

  // load input vector
  for (int i = 0; i < input_dim; ++i)
    m_input_layer->Neuron(i)->setOutput(input_vector(i));

  // propagate input vector through the neurocontroller
  CNCLayer *act_layer = 
    m_first_hidden_layer ? m_first_hidden_layer : m_output_layer;

  do {
    act_layer->Propagate();
    act_layer = act_layer->getNextLayer();
  } while(act_layer);

  // create output vector
  for (int j = 0; j < output_dim; ++j)
    output_vector(j) = m_output_layer->Neuron(j)->getOutput();
} // void CNC::propagate(const Vector& input_vector, Vector& output_vector)


//------------------------------------------------------------------------------
template<class T, class U>
void CNC::propagate(
    const T& input
   ,const unsigned int first
   ,const unsigned int last
   ,      U& output
)
{
  // first check for valid indices
  if (last<=first)
    FatalError("CNC::propagate(const T& input,....): last<first !\n");

  // check if input scope
  if ((last-first+1)!=input_dim)
    FatalError("CNC::propagate(const T& input,....): input data dimension not\
 equal to input neurons count!\n");

  // check if end index is valid
  if (last>=input.size())
    FatalError("CNC::propagate(const T& input,....): last index >= input dimension\n");

  // check whether input scope described by start and end indices is valid
  if ((last-first+1)>input.size())
    FatalError("CNC::propagate(const T& input,....): last index - first index >=\
 input dimension\n");

  // check output data structure dimension
  if (output.size()<m_output_layer->getNumberOfNeurons())
    FatalError("CNC::propagate(const T& input,....): output dimension to small");

  // load input vector
  for (uint i=0,uint ii=first; i < input_dim; ++i,++ii)
  {
    m_input_layer->Neuron(i)->setOutput((ii<=last)?input(ii):0.0);
  }

  // get neural net output and store in given output data structure
  propagate();                      // compute output neurons
  output.resize(output_dim);        // resize output data structure
  for (i = 0; i < output_dim; ++i)
    output_vector(i) = m_output_layer->Neuron(i)->getOutput();   
}



//------------------------------------------------------------------------------
double CNC::getOutput(const Vector& input_vector
                     ,const int output_neuron)
{
  // check input vector and output dimension

  int input_dim = input_vector.getDimension();
#ifdef development_version
  if(input_dim != m_input_layer->getNumberOfNeurons()) FatalError("CNC::getOutput(): neurocontroller input layer has different dimension than input vector");
#endif

#ifdef development_version
  if(output_neuron >= m_output_layer->getNumberOfNeurons()) FatalError("CNC::getOutput(): neurocontroller output layer has less neurons");
#endif

  // load input vector
  for (int i = 0; i < input_dim; i++)
    m_input_layer->Neuron(i)->setOutput(input_vector(i));

  // propagate input vector through the neurocontroller

  CNCLayer *act_layer = m_first_hidden_layer ? m_first_hidden_layer : m_output_layer;

  do {
    act_layer->Propagate();
    act_layer = act_layer->getNextLayer();
  } while(act_layer);

  return m_output_layer->Neuron(output_neuron)->getOutput();
};


//------------------------------------------------------------------------------
void CNC::load(const Vector& parameter_vector) 
{
#ifdef development_version
  if(parameter_vector.getDimension() != m_number_of_parameters + m_number_of_extra_parameters)
    FatalError("CNC::Load(): input NC does not fit the size of actual NC");
#endif

  m_parameter_index = 0;
  do 
  {
    *m_parameter[m_parameter_index] = parameter_vector(m_parameter_index);
  } while (++m_parameter_index != m_number_of_parameters);
};


//------------------------------------------------------------------------------
void CNC::load(const MathLib::Vector& param
              ,const int first
              ,const int last)
{
  const int sz(last-first+1); // number of parameters to load

#ifdef development_version
  if (sz<0) FatalError("CNC::Load(const Vector&, const int, const int): sz<0\n");
  if (first<0) FatalError("CNC::Load(const Vector&, const int, const int): first<0\n");
  if (last<=0) FatalError("CNC::Load(const Vector&, const int, const int): last<=0\n");
  if (last>param.getDimension()) 
    FatalError("CNC::Load(const Vector&, const int, const int): last> vector dimension\n");
  if(sz != m_number_of_parameters + m_number_of_extra_parameters)
    FatalError("CNC::Load(): input NC does not fit the size of actual NC");
#endif
  
  for (int i=0;i<sz;++i) *m_parameter[i]=param[i];
} // void CNC::Load(const MathLib::Vector& parameter_vector,const in first)


//------------------------------------------------------------------------------
void CNC::load(const char* filename) {
  ifstream in(filename);
  Assure(in, filename);

  // check if filetype is correct

  string filetype = ReadString(in);
  if(filetype != "neurocontroller") FatalError("CNC::Load(): wrong input file type");

  // check number of neurons in the input layer
  if(static_cast<int>(ReadNumber(in)) != m_number_of_input_neurons) FatalError("CNC::Load(): input NC does not fit the size of actual NC");

  // check number of hidden layers
  if(static_cast<int>(ReadNumber(in)) != m_number_of_hidden_layers) FatalError("CNC::Load(): input NC does not fit the size of actual NC");

  // check number of neurons in the hidden layers
  for (int i = 0; i < m_number_of_hidden_layers; i++) {
    if(static_cast<int>(ReadNumber(in)) != m_number_of_hidden_neurons[i]) FatalError("CNC::Load(): input NC does not fit the size of actual NC");
  }

  // check number of neurons in the output layer
  if(static_cast<int>(ReadNumber(in)) != m_number_of_output_neurons) FatalError("CNC::Load(): input NC does not fit the size of actual NC");

  // load values
  m_parameter_index = 0;

  do {
    *m_parameter[m_parameter_index] = ReadNumber(in);
  } while (++m_parameter_index != m_number_of_parameters);
};


//------------------------------------------------------------------------------
void CNC::save(const char* filename) {
  ofstream out(filename);
  Assure(out, filename);

  out << "neurocontroller" << endl;

  // save neurocontroller topology
  out << m_number_of_input_neurons << endl;
  out << m_number_of_hidden_layers << endl;
  for (int i = 0; i < m_number_of_hidden_layers; i++) {
    out <<  m_number_of_hidden_neurons[i] << endl;
  }
  out << m_number_of_output_neurons << endl;

  // save values
  m_parameter_index = 0;

  do {
    out << *m_parameter[m_parameter_index] << endl;
  } while (++m_parameter_index != m_number_of_parameters);
};
