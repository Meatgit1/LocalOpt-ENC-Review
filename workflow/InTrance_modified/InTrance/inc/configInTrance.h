/*!
\file configInTrance.h
\version 1.7
\brief InTrance and InTranceLib compilation configuration
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 06.01.2005 (AO)

\b Purpose:
    Configuration in this central header file should ease compilation
    and linking of either a singlethreaded executable InTrance file
    or a multithreaded DLL version of InTrance for use with a JAVA
    graphical user interface (GUI).

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.

\b Copyright:
   (c) 2001-2005 Bernd Dachwald, Andreas Ohndorf
*/

#ifndef CONFIGINTRANCE_H_INCLUDED
#define CONFIGINTRANCE_H_INCLUDED


#ifdef _MSC_VER    // Microsoft C++ Compiler

 #pragma warning (disable: 4786) // "identifier was truncated to '255' characters in the browser information"
 #pragma warning (disable: 4250) // "inherits ??? via dominance"

 #ifndef WINDOWS
  #define WINDOWS
 #endif

#endif //_MSC_VER

#if defined (__BORLANDC__) & !(defined BORLAND)
 #define BORLAND 1
#endif //__BORLANDC__

#ifdef JAVA_GUI
 #if defined (WIN32) || defined (_WIN32)    // Windows version
   #define COMMONCPP 1  // multithreading library for windows
 #endif // WIN32
#endif //JAVA_GUI

#if defined _DEBUG
 #define DEBUG 1
#endif //_DEBUG

#endif //CONFIGINTRANCE_H_INCLUDED
