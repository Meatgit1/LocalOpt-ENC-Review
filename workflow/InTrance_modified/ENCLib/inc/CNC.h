/*!
\file CNC.h
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
  andreas.ohndorf@gmx.de

\date

  09.05.2006 (AO)

\b Compiler:

  GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++

\b Copyright:

  (c) 2001-2006 by Bernd Dachwald, Andreas Ohndorf

  This software is protected by national and international copyright.
  Any unauthorized use, reproduction or modificaton is unlawful and
  will be prosecuted. Commercial and non-private application of the
  software in any form is strictly prohibited unless otherwise granted
  by the authors.

  The code is provided without any warranty; without even the implied
  warranty of merchantibility or fitness for a particular purpose.
*/
#ifndef CNC_H_INCLUDED
#define CNC_H_INCLUDED

// SpaceLib
#include <MLVector.h>

// ENCLib
#include <configENC.h>
#include <CNCTypes.h>

// forward declarations
class CInputNC;
class CNCLayer;
class CNC;


/*!\class CNCNeuron 
 * \brief neuron
 */
class CNCNeuron
// ----------------------------------------------------------------------------
{
  //-----------------------------------------------------------------
  int         m_number_of_links;
              ///< number of links to the connected (previous) layer

  /*!\brief pointer to the array of weights (dimension = number of neurons in the connected layer)
   *
   */
  double*     m_weight;

  double      m_bias; ///< bias value of the neuron
              
      public: 
        /// get neuron bias value
        double getBias(void) const { return m_bias; };
        /// set neuron bias value
      public: void setBias(const double value) { m_bias = value; };
private:

  double      m_temperature_exponent;  ///< temperature exponent of the neuron
  double      m_input; ///< neuron input value
  double      m_output; ///< neuron output value
public: 
  /// return neuron output value
  double getOutput(void) const { return m_output; };
private:

  NCLayerType m_layer_type; ///< type of layer in which the neuron is in

public: 
  /*! \brief return type of neuron layer
   * \return neuron layer type
   */
  NCLayerType getLayerType(void) const { return m_layer_type; };
private:

  /*! \brief pointer to the connected (previous) layer for this neuron
   * (i.e. the layer where the input comes from)
   */
  CNCLayer*   m_connected_layer;

  ///pointer to neurocontroller
  CNC*        m_nc;

//=============================================================================
// CNCNeuron CLASS INTERFACE
//=============================================================================
public:

  /*!\brief  layer type/neurocontroller constructor
   * Constructs a neuron from given information.
   * \param ltype type of the layer in which the neuron is in
   * \param connected_layer pointer to the connected (previous) neuron layer,
   * or NULL if an input neuron is to be constructed
   * \param nc pointer to neurocontroller
   */
  CNCNeuron(const NCLayerType ltype,
                  CNCLayer* connected_layer,
                  CNC* nc);

  /// destructor
  ~CNCNeuron() {SafeDeleteArray(m_weight);}

  /*! \brief  connection weight to the <index>-th neuron in the connected (previous)
   * neuron layer
   * \param index number of connected neuron from previous layer
   * \return desired connection weight value
   */
  double getWeight(const int index) const;

  /*!\brief set connection weight for given input connection 
   * \param index input connection index
   *  \param value new connection weight value
   */                   
  void setWeight(const int index,
                 const double value);

  /*!\brief set neuron output value (only possible for input neurons)
   * \param output new neuron output value */
  void setOutput(const double output);


  /*! \brief propagate neuron input through the neuron */
  void Propagate(void);

}; // end class CNCNeuron



/*! \brief neuron layer
 *
 */
class CNCLayer
{
  int         m_index; ///< layer index
public: 
  /// layer index
  int getIndex(void) const { return m_index; };
private:

  int         m_number_of_neurons; ///< number of layer neurons
public:
  /// layer neurons count
  int getNumberOfNeurons(void) const { return m_number_of_neurons; };
private:

  /// pointer to pointer array of neuron layers (dimension = number of neurons in the layer)
  CNCNeuron** m_neuron;
public:
  /*! \brief return neuron pointer
   * \param index index of neuron to get
   */
  CNCNeuron* Neuron(const int index) const { return m_neuron[index]; };
private:

  CNCLayer*   m_prev_layer; ///< connection (pointer) to previous layer
public: 
  /*! \brief get pointer to previous layer*/
  CNCLayer* getPrevLayer(void) const { return m_prev_layer; };
private:

  CNCLayer*   m_next_layer; ///< pointer to the next layer (NULL, if there is any)
public:
  /*! \brief get pointer of next layer */
  CNCLayer* getNextLayer(void) const { return m_next_layer; };

  //-----------------------------------------------------------------

//=============================================================================
// CNCLayer CLASS INTERFACE
//=============================================================================
public:

  /*! \brief constructor
   *  \param index layer index
   *  \param number_of_neurons number of layer neurons
   *  \param ltype layer type (input, hidden, output)
   *  \prev previous neuron layer
   *  \param nc neurocontroller
   */
  CNCLayer(int index
          ,int number_of_neurons
          ,const NCLayerType ltype
          ,CNCLayer* prev
          ,CNC* nc
          );

  /// destructor
  ~CNCLayer();

  /// propagate neuron layer input through the neuron layer
  //void Propagate(void) const;
  inline void Propagate(void) const 
  {
    for (int i = 0; i < m_number_of_neurons; ++i) m_neuron[i]->Propagate();
  }

}; // end class CNCLayer




class CNeuralNet
{
  
public:
  CNeuralNet(void) {};
  CNeuralNet(const CInputNC& input_ncp) throw(std::runtime_error&);
  virtual ~CNeuralNet(void) {};

  virtual NCInputType getInputType(void) const = 0;
  
  virtual int getNumberOfInputNeurons(void) const = 0;
  virtual int getNumberOfLinks(void) const = 0;
  virtual int getNumberOfNeurons(void) const = 0;
  virtual int getNumberOfOutputNeurons(void) const = 0;
  virtual int getNumberOfParameters(void) const = 0;
  virtual int getParameterNeuronNumber(const int i) const = 0;
  virtual unsigned int getType(void) const = 0;
  
  virtual void load(const char* filename) = 0;
  virtual void load(const MathLib::Vector& parameters) = 0;
  virtual void load(const MathLib::Vector& parameters, const int first, const int last) = 0;
  virtual void propagate(const MathLib::Vector& input_vector, MathLib::Vector& output_vector) = 0;
  virtual void setParameter(const int index, double* value) = 0;
  virtual std::vector<unsigned int>& getNeuronParameterStartIndices(void) = 0;
  virtual std::vector<unsigned int>& getConnectionWeightsIndices(void) = 0;
  virtual std::vector<unsigned int>& getNeuronActFctParameterIndices(void) = 0;
  

  virtual NCActivationFunction getActivationFunction(void) const = 0;
  virtual int getParameterIndex(void) const = 0;
  virtual int getNumberOfExtraParameters(void) const = 0;
  virtual void setParameterNeuronNumber(void) = 0;
  virtual void incNeuronNumber(void) = 0;
  virtual void incParameterIndex(void) = 0;

  template<class T, class U>
  void propagate(
    const T& input
   ,const unsigned int first
   ,const unsigned int last
   ,      U& output);

  virtual void save(const char* filename) = 0;  
};

/*! Tbd.*/
class CLayeredNeuralNet : virtual public CNeuralNet
{
public:
  CLayeredNeuralNet(void);
  CLayeredNeuralNet(const CLayeredNeuralNet& copy);
  CLayeredNeuralNet& operator=(const CLayeredNeuralNet& rhs);
  ~CLayeredNeuralNet(void);
};


/*! Tbd.*/
class CFeedForwardNeuralNet : virtual public CNeuralNet
{
public:
  CFeedForwardNeuralNet(void);
  CFeedForwardNeuralNet(const CFeedForwardNeuralNet& copy);
  CFeedForwardNeuralNet& operator=(const CFeedForwardNeuralNet& rhs);
  ~CFeedForwardNeuralNet(void);
};





/*! \brief neurocontroller
 *
 */
class CNC
{

  unsigned int         m_nc_type;    ///< application specific nc type identifier

public:
  /// get neurocontroller type
  unsigned int getType(void) const {return m_nc_type;};
private:

  NCInputType          m_input_type; ///< input neuron data type
                       // neurocontroller input type
      public:          
        /// neuron input type
        NCInputType getInputType(void) const { return m_input_type; };
      private:

  NCActivationFunction m_activation_function; ///< neuron activiation function

      public:          
        /// neuron activation function
        NCActivationFunction getActivationFunction(void) const { return m_activation_function; };
      private:

  int                  m_number_of_hidden_layers; ///< number of hidden layers                       

  int                  m_number_of_input_neurons; ///< number of input neurons
                       
      public:          
        /// number of input neurons
        int getNumberOfInputNeurons(void) const { return m_number_of_input_neurons; };
      private:

  /// pointer to an array containing the number of hidden neurons
  int*                 m_number_of_hidden_neurons;
                       // pointer to an array containing the number of hidden neurons

  int                  m_number_of_output_neurons; ///< number of output neurons
                       
      public:          
        /// number of output neurons
        int getNumberOfOutputNeurons(void) const { return m_number_of_output_neurons; };
      private:

  int                  m_number_of_layers; ///< total number of layers

  int                  m_number_of_neurons; ///< total number of neurons
                       
      public:          
        /// total number of neurons
        int getNumberOfNeurons(void) const { return m_number_of_neurons; };
      private:

  /// extra parameters coded at the end of chromosome
  int                  m_number_of_extra_parameters;
                       
      public:          
        /// extra parameters coded at the end of chromosome
        int getNumberOfExtraParameters(void) const {return m_number_of_extra_parameters;};
        private:

  int                  m_number_of_links; ///< total number of links
                       
      public:         
        /// total number of links
        int getNumberOfLinks(void) const { return m_number_of_links; };
      private:

  /// total number of parameters (non-input layer-neurons + their weights)
  int                  m_number_of_parameters;
                       
      public:          
        /// total number of parameters
        int getNumberOfParameters(void) const { return m_number_of_parameters; };
      private:

  /// pointer to the input layer
  CNCLayer*            m_input_layer;
                       
  /// pointer to the first (and maybe only) hidden layer
  CNCLayer*            m_first_hidden_layer;                       

  /// pointer to the output layer
  CNCLayer*            m_output_layer;


  /// apply neural net to input neurons
  inline void propagate(void)
  {
    CNCLayer *act_layer = m_first_hidden_layer ? m_first_hidden_layer : m_output_layer;

    do {
      act_layer->Propagate();
      act_layer = act_layer->getNextLayer();
    } while(act_layer);
  }

   //-----------------------------------------------------------------

//=============================================================================
// CNC CLASS INTERFACE
//=============================================================================
private:

  /// pointer to an array of pointers to the neurocontroller parameters
  double** m_parameter;
           
  /// pointer to an array containing the neuron number to each neurocontroller
  /// parameter
  int*     m_parameter_neuron_number;

  /// auxiliary index to move through the array of pointers above
  /// to be global in CNC
  int      m_parameter_index;

  /// auxiliary variable to be global in CNC
  int      m_neuron_number;
public:

  virtual inline int getParameterIndex(void) const
  {
    return m_parameter_index;
  }

  virtual inline void incNeuronNumber(void) 
  {
    ++m_neuron_number;
  }

  virtual inline void incParameterIndex(void)
  {
    ++m_parameter_index;
  }
  
  virtual inline int getParameterNeuronNumber(const int i) const 
  {
    return m_parameter_neuron_number[i];
  }

  virtual inline void setParameterNeuronNumber(void)
  {
    m_parameter_neuron_number[m_parameter_index] = m_neuron_number;
  }

  virtual inline void setParameter(const int index, double* value)
  {
    m_parameter[index]=value;
  }


  //-----------------------------------------------------------------

  /// configuration input constructor
  CNC(const CInputNC& input_ncp) throw(std::runtime_error&);

  ~CNC(void); ///< destructor

  // __functionality__

  /*! \brief apply neural net to input values
   *
   * Propagate <input_vector> through the neurocontroller, i.e.
   * create <output_vector> from <input_vector>.
   *
   * \param input_vector
   * \param output_vector
   */
  void propagate(const MathLib::Vector& input_vector
                ,      MathLib::Vector& output_vector);


  /*! \brief apply neural net to input values (given data scope)
   *
   */
  template<class T, class U>
  void propagate(
    const T& input
   ,const unsigned int first
   ,const unsigned int last
   ,      U& output);

  

  /*!\brief get neurocontroller output data (vector)
   * \param input_neuron
   * \param output_neuron
   * \return
   */
  double getOutput(const MathLib::Vector& input_vector
                  ,const int              output_neuron);

  /*! \brief set neurocontroller input data (vector)
   *
   * Load neurocontroller parameters from vector <parameter_vector>.
   */
  void load(const MathLib::Vector& parameter_vector);

  /*!\brief Set neurocontroller with parameters from given vector.
   * 
   * The data of given vector (starting from index first to the one
   * at index last) are used to set the neurocontroller's neural
   * net
   *
   * \param parameter_vector Data vector with new neurocontroller settings.
   * \param first Start index of data to use from given data vector.
   * \param last  End index of data to use from given data vector.
   *
   */
  void load(const MathLib::Vector& parameter_vector
           ,const int first
           ,const int last);

  /*! \brief set neural net with input file content
   *
   * Load neurocontroller parameters from file
   *\param filename name of data file
   */
  void load(const char* filename);


  /*!\brief save neurocontroller parameters to file <filename>
   * 
   * \param filename name of file to save neural net parameters to
   */
  void save(const char* filename);

}; // end class CNC


#endif // CNC_H_INCLUDED
