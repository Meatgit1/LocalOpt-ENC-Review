/*!
\file MathLib.h
\version 1.7
\brief The header file that you can use to include the whole MathLib.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 04.02.2005 (BD)

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.

\b Copyright:
   (c) 2001-2004 Bernd Dachwald, Andreas Ohndorf
*/

#ifndef MATHLIB_H
#define MATHLIB_H

#include <math.h>
#include <fstream>
#include <iomanip>

#include <MLConst.h>
#include <MLGlobal.h>
#include <MLPolynomial.h>
#include <MLQuaternion.h>
#include <MLDE.h>

/** \namespace MathLib
    \brief The namespace of the MathLib functions.
    All the classes and functions of MathLib are encapsulated in the namespace MathLib.
    So if you want to use the MathLib it is necessary to include the "using"-statement
	  in your implementation file. */

using namespace MathLib;

#endif // MATHLIB_H