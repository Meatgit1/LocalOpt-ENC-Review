/*!
\file InTrance.cpp
\version 1.7
\brief InTrance main function.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 29.10.2006 (AO)

\b Purpose:
   This file is the main file for InTrance.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.

\b Copyright:
   (c) 2001-2006 Bernd Dachwald, Andreas Ohndorf
*/
//#define _CRTDBG_MAP_ALLOC  
//#include <stdlib.h>  
//#include <crtdbg.h>
// InTrance
#include <configInTrance.h>
#include <CInTrance.h>
#include <ExitDiagnostics.h>
#include <NativeBridgeProbe.h>
#include <iostream>
#include <string>
#include <memory>

// The following identifier should be used later to distinguish between the
// public standard version and the confidential research version. 
// The standard version will not contain the commands and the 
// functionalities that are part of the authors' current research.
#define RESEARCH_VERSION

// The launcher can identify this build without opening mission inputs or data.
static void printUnifiedCapabilities()
{
#if defined(_WIN32) && defined(_M_IX86)
  const char* native_build = "x86_pe32";
#else
  const char* native_build = "unsupported";
#endif
  std::cout << "{\"schema\":\"intrance_unified_v1\","
    "\"modes\":[\"search\",\"export\",\"schedule\"],"
    "\"profiles\":[\"mercury_body_rendezvous_v1\",\"mercury_body_rendezvous_precision_v3\"],"
    "\"observer_schema\":1,\"native_build\":\"" << native_build << "\","
    "\"pointer_bits\":" << 8 * sizeof(void*) << ","
    "\"fresh_process_per_stage\":true}\n";
}

int main(int argc, char* argv[])
{
  if (argc == 2 && std::string(argv[1]) == "--capabilities")
  {
    printUnifiedCapabilities();
    return 0;
  }
  InstallInTranceExitDiagnostics();
  cout << "                          ___  ";
#ifdef RESEARCH_VERSION
  cout << "Version 1.7.6" << endl;
#else
  cout << "Release 1.7" << endl;
#endif
  cout << "  ___                     \\  | build "<< __DATE__ << endl;
  cout << " |   |___________________/ /\\|       "<< __TIME__ << endl;
  cout << " |   |_______   ________ \\/ ____  _____  ____"        << endl;
  cout << " |   |/    \\|   |  _ \\/ \\ \\/    \\/  ___\\/ __ \\" << endl;
  cout << " |   |   |  \\   |  |\\/ /_\\ \\  |  \\  \\__\\  ___/" << endl;
  cout << " |___|___|__/___|__|(_______)_|__/\\_____\\_____>"     << endl << endl;
  cout << " INtelligent spacecraft TRAjectory optimization"       << endl;
  cout << "        using NeuroController Evolution"               << endl << endl;
  cout << " (c) 2001-2006    B e r n d   D a c h w a l d"         << endl;
  cout << "                  German Aerospace Center (DLR)"       << endl;
  cout << "                  Institute of Space Simulation"       << endl;
  cout << "                  Linder Hoehe, 51170 Cologne"         << endl;
  cout << "________________________________________________"      << endl << endl;

  int result = -1;

  if (argc == 2)
  {

    try 
    {
      // create the InTrance object using the input file
      std::auto_ptr<CInTrance> intrance(new CInTrance(argv[1]));
      Assure(intrance.get(), "main()");

      intrance->run(); // run InTrance for the specified input file
      NativeBridgeProbe::finish();
      std::cerr << "EXIT_TRACE main: before delete intrance" << std::endl;
      intrance.reset();
      std::cerr << "EXIT_TRACE main: after delete intrance" << std::endl;
      result = 0;
    }
    catch (std::exception& e)
    {
      std::cout << e.what() << std::endl;
    }

  }
  else
  {
    cout << endl << "usage: InTrance <input_file>" << endl;
  }

  std::cerr << "EXIT_TRACE main: returning " << result << std::endl;
  return result;
};
