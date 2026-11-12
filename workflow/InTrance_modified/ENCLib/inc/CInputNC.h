/*!
\file CInputNC.h
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

  30.05.2005 (AO)

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
#ifndef CINPUTNC_H_INCLUDED
#define CINPUTNC_H_INCLUDED

// SpaceLib
#include <Global.h>

// ENCLib
#include <configENC.h>
#include <CNCTypes.h>

const unsigned int DEFAULT_CONTROLLER_TYPE(0);

/// neurocontroller input object class
class CInputNC 
{
//=============================================================================
// CInputNC CLASS INTERFACE
//=============================================================================
public:

  unsigned int         controller_type;      ///< application specific type

  NCInputType          input_type;           ///< neurocontroller input type
  unsigned int         input_neurons;        ///< number of input neurons
  unsigned int         hidden_layers;        ///< number of hidden layers
  int*                 hidden_neurons;       ///< pointer to an array containing the number of hidden neurons
  unsigned int         output_neurons;       ///< number of output neurons
  NCActivationFunction activation_function;  ///< neural activation function
  unsigned int         extra_parameters;     ///< additional values at the end of chromosomes (they won't be used for NC)
  
  /// standard constructor
  CInputNC(void);

  /// assign operator
  const CInputNC& operator=(const CInputNC& to_copy);

  /*!
   * \brief Construct object from parameter list.
   *
   * Constructs a CInputNC object from a parameter vector (STL list)
   * which contains parameters to describe the neurocontroller
   * topology
   */
  CInputNC(
     ParamList& parameters
    ,bool extra_data = true
    ,const std::string& trailing_chars = " \t"
  );

  /*!
   * \brief direct constructor for use without parameter list.
   */
  CInputNC( const NCInputType          p_input_type
           ,const int                  p_input_neurons
           ,const int                  p_hidden_layers
           ,      int*                 p_hidden_neurons
           ,const int                  p_output_neurons
           ,const NCActivationFunction p_AcFunc
           ,const int                  p_extra_parameters);

  /// Sets number of additional input neurons.
  void setExtraParameters(const int p_extra_parameters) { extra_parameters = p_extra_parameters;};

  /// Returns  number of additional input neurons.
  int  getExtraParameters(void) const { return extra_parameters;};

  ///brief prints actual settings through given output stream object
  void printSettings(void) const;

  /// resets input object
  void reset(void);

  /// destructor
  ~CInputNC(void) {SafeDeleteArray(hidden_neurons);};

}; // end class CInputNC


///// Stream operator for CInputNC objects.
//inline std::ostream& operator << (std::ostream& s, const CInputNC& i_nc)
//{
//  i_nc.printSettings(s);
//  return s;
//};

#endif // CINPUTNC_H_INCLUDED
