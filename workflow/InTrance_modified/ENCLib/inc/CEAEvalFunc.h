/*!
\file CEAEvalFunc.h
\brief A virtual class for the evaluation function of an evolutionary algorithm.

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

  06.01.2005 (AO)

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
#ifndef CEAEVALFUNC_H_INCLUDED
#define CEAEVALFUNC_H_INCLUDED

// ENCLib
#include <configENC.h>

//forward declaration
class CEAFitness;
class CIndividual;


// evolutionary algorithm evaluation function class declaration
// ----------------------------------------------------------------------------
class CEAEvalFunc 
{

public:

  // virtual base function
  virtual CEAFitness* evaluate(     CIndividual&  solution
                              ,const bool         final    = false
                              ,const bool         monitor  = false
                              ,const double       time     = 0.0)  = 0;

}; // end class CEAEvalFunc


#endif // CEAEVALFUNC_H_INCLUDED
