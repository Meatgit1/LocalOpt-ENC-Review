/*!
\file CNCTypes.h
\brief A header that provides type declarations for the neurocontroller.

\version

  1.7

\author

  Bernd Dachwald (AO)\n
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

  03.06.2004 (AO)

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
#ifndef CNCTYPES_H_INCLUDED
#define CNCTYPES_H_INCLUDED


enum NCInputType // available NC input types
{
  IT_UNDEF,  // input type undefined
  C_INPUT,   // cartesian coordinates
  P_INPUT,   // polar coordinates
  E_INPUT,   // orbital elements
  CP_INPUT,  // cartesian + polar coordinates
  CE_INPUT,  // cartesian coordinates + orbital elements
  PE_INPUT,  // polar coordinates + orbital elements
  CPE_INPUT, // cartesian + polar coordinates + orbital elements
  Q_INPUT    // quaternion/q.derivative input
};

enum NCActivationFunction // available NC activation functions
{
  AF_UNDEF,            // activation function undefined
  SIGMOID,             // sigmoid activation function
  TANH,                // hyperbolic tangens
  SIMPLE,              // simple activation function (ref. Diss. Goerke)
  SIGMOIDIDENTITY,     // sigmoid activation function in the hidden layers and identity in the output layer
  TANHIDENTITY,        // hyperpolic tangens in the hidden layers and identity in the output layer
  SIMPLEIDENTITY       // simple activation function (ref. Diss. Goerke) in the hidden layers and identity in the output layer
};

enum NCLayerType // NC layer types
{
  INPUT_LAYER,  // input layer
  HIDDEN_LAYER, // hidden layer
  OUTPUT_LAYER  // output layer
};


#endif // CNCTYPES_H_INCLUDED
