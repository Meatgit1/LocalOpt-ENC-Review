/*!
\file CSolarSystem.cpp
\version 1.7
\brief A class that provides the functionality of the whole solar system (CSolarSystem).
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Tanja D. Schmidt <tanja.schmidt@dlr.de> (TS)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 07.11.2006(AO)

\b Purpose:
   A CSolarSystem object is intended to keep the main bodies of the solar
   system in its internal data structures and provide easy access to all the
   bodies data. Additionaly there are some useful functions implemented, e.g.
   the HILL-Radius of a solar system body.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.
   This code is based on the code from:
   Montenbruck, Gill: Satellite Orbits. Springer, 2000 and
   Montenbruck, Pfleger: Astronomie mit dem Personal Computer (3. Aufl.). Springer, 1999

\b Copyright:
   (c) 2001-2006 Bernd Dachwald, Tanja D. Schmidt, Andreas Ohndorf
*/

//SpaceLib
#include <configSpaceLib.h>

//STL
#include <sstream>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <new>

// standard C library
#include <math.h>

//SpaceLib
#include <MLDE.h>
#include <Kepler.h>
#include <CBody.h>
#include <CBodyData.h>
#include <CSolarSystem.h>
#include <CSolarSystemData.h>

using namespace std;
using namespace Global;
using namespace MathLib;


namespace SpaceLib 
{

  const std::string S_NOT_APPLICABLE = "n/a";

  //----------------------------------------------------------------------------
  CSolarSystem::CSolarSystem(const string& DE405FileName)throw (std::bad_alloc&)
  {
    const std::string method_name = "CSolarSystem::CSolarSystem(const string&)";
    const std::string err_msg_01 = ": failed to allocate memory for solar \
                                   system data!\n";
    const std::string error_message;
    
    solarSystemData = 0;                            // creation of internal data
    solarSystemData = new CSolarSystemData;                        // structures
    if (solarSystemData == nullptr) {                                 // error handling
        std::string error_message = method_name + err_msg_01;
        //throw std::make_exception_ptr(std::bad_alloc(error_message.c_str()));
        cout << error_message << endl;
    }

    solarSystemData->DE405File=DE405FileName;   // store name of DE405 data file

    bool result=false;

    solarSystemData->myDE->setDataFile(DE405FileName);
    try
    {
      result=solarSystemData->myDE->initialize();
    }
    catch (FileErrorEXC& err)
    {
      err.PrintEXC();
    }

    if (result==true)
    {

      CBody*ssBary=0;
      ssBary=new CBody;
      if (ssBary!=0) {
        ssBary->InitializeDE405(SS_BARY, SS_BARY, ssBary, DE405_SS_BARY,
             SS_BARY_NAME, this, SPECIAL, solarSystemData->myDE);
        insertBody(ssBary);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (ss_bary)");

      CBody*sun=new CBody;
      if (sun!=0) {
        sun->InitializeDE405(SUN, SS_BARY, ssBary, DE405_SUN, SUN_NAME, this,
             BODY, solarSystemData->myDE);
        insertBody(sun);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (sun)");


      CBody*temp=new CBody;
      if (temp!=0) {
        temp->InitializeDE405(MERCURY, SUN, sun, DE405_MERCURY, MERCURY_NAME,
             this, BODY, solarSystemData->myDE);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (mercury)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->InitializeDE405(VENUS, SUN, sun, DE405_VENUS, VENUS_NAME, this,
             BODY, solarSystemData->myDE);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (venus)");


      CBody*earth=0;
      earth=new CBody;
      if (earth!=0) {
        earth->InitializeDE405(EARTH, SUN, sun, DE405_EARTH, EARTH_NAME, this,
             BODY, solarSystemData->myDE);
        insertBody(earth);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (earth)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->InitializeDE405(LUNA, EARTH, earth, DE405_MOON, LUNA_NAME, this,
             BODY, solarSystemData->myDE);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (moon)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->InitializeDE405(MARS, SUN, sun, DE405_MARS, MARS_NAME, this,
             BODY, solarSystemData->myDE);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (mars)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->InitializeDE405(JUPITER, SUN, sun, DE405_JUPITER, JUPITER_NAME,
             this, BODY, solarSystemData->myDE);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (jupiter)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->InitializeDE405(SATURN, SUN, sun, DE405_SATURN, SATURN_NAME,
             this, BODY, solarSystemData->myDE);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (saturn)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->InitializeDE405(URANUS, SUN, sun, DE405_URANUS, URANUS_NAME,
             this, BODY, solarSystemData->myDE);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (uranus)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->InitializeDE405(NEPTUNE, SUN, sun, DE405_NEPTUNE, NEPTUNE_NAME,
             this, BODY, solarSystemData->myDE);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (neptune)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->InitializeDE405(PLUTO, SUN, sun, DE405_PLUTO, PLUTO_NAME, this,
             BODY, solarSystemData->myDE);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (pluto)");


    } else {

      CBody*ssBary=0;

      ssBary=new CBody;
      if (ssBary!=0) {
        ssBary->CState::SetState(Vector(6),
             CARTESIAN, MJD_J2000, SS_BARY, AUD, SPECIAL, 0,
             DE405_SS_BARY, ssBary);
        ssBary->planetData->m_index=SS_BARY;
        ssBary->setName(SS_BARY_NAME);
        ssBary->setRefSolSys(this);
        insertBody(ssBary);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (ss_bary)");


      CBody*sun=0;

      sun=new CBody;
      if (sun!=0) {
        sun->Initialize(SUN, SUN_NAME, 0, MJD_J2000, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0);
        sun->CState::SetState(Vector(6),
             CARTESIAN, MJD_J2000, SUN, AUD, BODY, 0, DE405_SUN, sun);
        sun->setName(SUN_NAME);
        sun->setRefSolSys(this);
        sun->setGM(GM_SUN_AUD, AUD);
        insertBody(sun);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (sun)");


      CBody*temp=0;

    temp=new CBody;
      if (temp!=0) {
        temp->Initialize(MERCURY, MERCURY_NAME, 3.302e23, MJD_J2000,
             0.38709893, 0.20563069, 7.00487*DEG2RAD, 48.33167*DEG2RAD,
             77.45645*DEG2RAD, 252.25084*DEG2RAD, 0.00000066*DAY2JCY,
             0.00002527*DAY2JCY, -23.51*ARCSEC2RAD*DAY2JCY,
             -446.30*ARCSEC2RAD*DAY2JCY, 573.57*ARCSEC2RAD*DAY2JCY,
             538101628.29*ARCSEC2RAD*DAY2JCY, ELEMENTS, sun);
        temp->setRefSolSys(this);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (mercury)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->Initialize(VENUS, VENUS_NAME, 48.685e23, MJD_J2000, 0.72333199,
             0.00677323, 3.39471*DEG2RAD, 76.68069*DEG2RAD, 131.53298*DEG2RAD,
             181.97973*DEG2RAD, 0.00000092*DAY2JCY, -0.00004938*DAY2JCY,
             -2.86*ARCSEC2RAD*DAY2JCY, -9196.89*ARCSEC2RAD*DAY2JCY,
             -108.80*ARCSEC2RAD*DAY2JCY, 210664136.06*ARCSEC2RAD*DAY2JCY,
             ELEMENTS, sun);
        temp->setRefSolSys(this);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (venus)");


      CBody*earth=0;

      earth=new CBody;
      if (earth!=0) {
// GTOC2
        earth->Initialize(EARTH, EARTH_NAME, 59.73698991e23, 54000,
          0.999988049532578
         , 0.01671681163160
         , 0.8854353079654e-3*DEG2RAD
         , 175.40647696473*DEG2RAD
         , (175.40647696473+287.61577546182)*DEG2RAD
         , (175.40647696473+287.61577546182+257.60683707535)*DEG2RAD
         , 0, 0, 0, 0, 0,
          129597740.63*ARCSEC2RAD*DAY2JCY, ELEMENTS, sun);
// end of GTOC2

        //earth->Initialize(EARTH, EARTH_NAME, 59.73698991e23, MJD_J2000,
        //     1.00000011, 0.01671022, 0.00005*DEG2RAD, 348.73936*DEG2RAD,
        //     102.94719*DEG2RAD, 100.46435*DEG2RAD, -0.00000005*DAY2JCY,
        //     -0.00003804*DAY2JCY, -46.94*ARCSEC2RAD*DAY2JCY,
        //     -18228.25*ARCSEC2RAD*DAY2JCY, 1198.28*ARCSEC2RAD*DAY2JCY,
        //     129597740.63*ARCSEC2RAD*DAY2JCY, ELEMENTS, sun);
        earth->setRefSolSys(this);
        insertBody(earth);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (earth)");

      // data for luna from Schmitt, Tanja; unprooved
      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->Initialize(LUNA, LUNA_NAME, 7.348e22, MJD_J2000, 384388631.0/AU,
             0.0549, 5.1454*DEG2RAD, 125.1228*DEG2RAD, 83.1862*DEG2RAD,
             198.5516*DEG2RAD, 0.0*DAY2JCY, 0.0*DAY2JCY, 0.0*DAY2JCY,
             -0.0529538083*DEG2RAD, 0.111403514*DEG2RAD, 13.1763964649*DEG2RAD,
             ELEMENTS, earth);

        temp->CState::SetState(384388631.0/AU, //a
             0.0549, //e
             5.1454*DEG2RAD, //i
             125.1228*DEG2RAD, //Omega
             83.1862*DEG2RAD, //CurliPi
             198.5516*DEG2RAD, //L
             ELEMENTS, MJD_J2000, EARTH, AUD);
        temp->setRefSolSys(this);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (luna)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->Initialize(MARS, MARS_NAME, 6.4185e23, MJD_J2000, 1.52366231,
             0.09341233, 1.85061*DEG2RAD, 49.57854*DEG2RAD, 336.04084*DEG2RAD,
             355.45332*DEG2RAD, -0.00007221*DAY2JCY, 0.00011902*DAY2JCY,
             -25.47*ARCSEC2RAD*DAY2JCY, -1020.19*ARCSEC2RAD*DAY2JCY,
             1560.78*ARCSEC2RAD*DAY2JCY, 68905103.78*ARCSEC2RAD*DAY2JCY,
             ELEMENTS, sun);
        temp->setRefSolSys(this);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (mars)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->Initialize(JUPITER, JUPITER_NAME, 18986.0e23, MJD_J2000,
             5.20336301, 0.04839266, 1.30530*DEG2RAD, 100.55615*DEG2RAD,
             14.75385*DEG2RAD, 34.40438*DEG2RAD, 0.00060737*DAY2JCY,
             -0.00012880*DAY2JCY, -4.15*ARCSEC2RAD*DAY2JCY,
             1217.17*ARCSEC2RAD*DAY2JCY, 839.93*ARCSEC2RAD*DAY2JCY,
             10925078.35*ARCSEC2RAD*DAY2JCY, ELEMENTS, sun);
        temp->setRefSolSys(this);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (jupiter)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->Initialize(SATURN, SATURN_NAME, 5684.6e23, MJD_J2000, 9.53707032,
             0.05415060, 2.48446*DEG2RAD, 113.71504*DEG2RAD, 92.43194*DEG2RAD,
             49.94432*DEG2RAD, -0.00301530*DAY2JCY, -0.00036762*DAY2JCY,
             6.11*ARCSEC2RAD*DAY2JCY, -1591.05*ARCSEC2RAD*DAY2JCY,
             -1948.89*ARCSEC2RAD*DAY2JCY, 4401052.95*ARCSEC2RAD*DAY2JCY,
             ELEMENTS, sun);
        temp->setRefSolSys(this);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (saturn)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->Initialize(URANUS, URANUS_NAME, 868.32e23, MJD_J2000,
             19.19126393, 0.04716771, 0.76986*DEG2RAD, 74.22988*DEG2RAD,
             170.96424*DEG2RAD, 313.23218*DEG2RAD, 0.00152025*DAY2JCY,
             -0.00019150*DAY2JCY, -2.09*ARCSEC2RAD*DAY2JCY,
             -1681.40*ARCSEC2RAD*DAY2JCY, 1312.56*ARCSEC2RAD*DAY2JCY,
             1542547.79*ARCSEC2RAD*DAY2JCY, ELEMENTS, sun);
        temp->setRefSolSys(this);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (saturn)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->Initialize(NEPTUNE, NEPTUNE_NAME, 1024.3e23, MJD_J2000,
             30.06896348, 0.00858587, 1.76917*DEG2RAD, 131.72169*DEG2RAD,
             44.97135*DEG2RAD, 304.88003*DEG2RAD, -0.00125196*DAY2JCY,
             0.00002514*DAY2JCY, -3.64*ARCSEC2RAD*DAY2JCY,
             -151.25*ARCSEC2RAD*DAY2JCY, -844.43*ARCSEC2RAD*DAY2JCY,
             786449.21*ARCSEC2RAD*DAY2JCY, ELEMENTS, sun);
        temp->setRefSolSys(this);
        insertBody(temp);
      } else
        FatalError("CSolarSystem::CSolarSystem error during memory allocation (saturn)");


      temp=0;
      temp=new CBody;
      if (temp!=0) {
        temp->Initialize(PLUTO, PLUTO_NAME, 0.15e23, MJD_J2000, 39.48168677,
             0.24880766, 17.14175*DEG2RAD, 110.30347*DEG2RAD,
             224.06676*DEG2RAD, 238.92881*DEG2RAD, -0.00076912*DAY2JCY,
             0.00006465*DAY2JCY, 11.07*ARCSEC2RAD*DAY2JCY,
             -37.33*ARCSEC2RAD*DAY2JCY, -132.25*ARCSEC2RAD*DAY2JCY,
             522747.90*ARCSEC2RAD*DAY2JCY, ELEMENTS, sun);
        temp->setRefSolSys(this);
        insertBody(temp);
      } // of if (temp != 0)
    } // end of else
  } // end CSolarSystem::CSolarSystem


  //---------------------------------------------------------------------------
  DE405*CSolarSystem::getEphemObj(void) 
  { 
    return solarSystemData->myDE; 
  }



  //---------------------------------------------------------------------------
  void CSolarSystem::setCometDataFile(const std::string& file_name) {

    ifstream data_file(file_name.c_str(), ios::in);

    if (data_file)
    {
      solarSystemData->cometFile = file_name;

#ifdef DEBUG
      cout << "CSolarSystem::setCometDataFile(std::string&)"
           << " =>using comet data file (" << file_name << ")" << endl;
      cout.flush();
#endif //DEBUG

    }
    else
    {
      throw invalid_argument("CSolarSystem::setCometDataFile(std::string&) =>comet data file (" + file_name + ") not available!\n");
    }

  }


  //-----------------------------------------------------------------------------
  /** Tbd. */
  DataField::DataType estimate_data_type(const string& s) {

    string::size_type pos=string::npos;
    DataField::DataType myType=DataField::NOTYPE;
    string myString(s);

    // remove possible spaces and convert to capitel letters
    while ((pos=myString.find(BLANK))!=string::npos) myString.erase(pos, 1);
    for (string::iterator j=myString.begin(); j<myString.end(); ++j) (*j)=toupper(*j);

    pos=string::npos;
    pos=myString.find_first_of("ABCDFGHIJKLMNOPQRSTUVWXYZ");
    int cntE=0;
    int cntP=0;
    if (pos!=string::npos) myType=DataField::TEXT; else {

      cntE=count(myString.begin(), myString.end(), 'E');
      cntP=count(myString.begin(), myString.end(), '.');

      switch (cntE) {
        case 0: // there's no E, like 1.0E-1 that indicates a floating point number
          switch (cntP) {
          case 0:
          myType=DataField::INTEGER;
          break; // if there's no decimal point => integer
          case 1:
          myType=DataField::FLOAT;
          break; // there's only one => floating point number
          default:
          myType=DataField::TEXT;
          break; // else it must be text
          }
          break;
        case 1: // ok, now we've got exactly one 'E' or 'e'
          switch (cntP) {
          case 1:
          myType=DataField::FLOAT;
          break; // is there a decimal point? =>float
          default:
          myType=DataField::TEXT;
          break; // anything else => text
          }
          break;
        default:
          myType=DataField::TEXT; // more than one 'e' or 'E' => text
      }
    }
    return myType;

  }

  // tries to parse data structure of a text file
  //-----------------------------------------------------------------------------
  /** Tbd. */
  FileStructure get_data_structure(const string& s1, const string& s2, const string& s3) {

    // remove possibly leading and trailing spaces and
    // estimate number of columns
    string ws = s2;
    while(*(ws.begin()) == BLANK) ws.erase(ws.begin());
    while(*(ws.end()-1) == BLANK) ws.erase(ws.end()-1);
    int noc = count(ws.begin(), ws.end(), ' ') + 1; // number of columns

    vector<DataField>myVec;               // initialize data structure

    if (ws.find_first_not_of(" -") == string::npos) {

      if (noc>0) {
        DataField myDataField;
        string::size_type pos=0, pos_old;

        for (int i=0; i<noc; i++) {
          pos_old=pos;             // save beginning of keyword
          pos=s2.find(BLANK, pos); // look for a space character

          // estimate position in string
          if ((pos==string::npos)&&(pos_old!=0)) pos=s2.size();
          myDataField.pos=pos_old;

          myDataField.width=pos-pos_old;  // set size of keyword

          // and copy string (the actual keyword)
          myDataField.name.assign(s1, pos_old, myDataField.width);

          // convert to capitel letters
          for (string::iterator j=myDataField.name.begin();
               j<myDataField.name.end(); ++j) (*j)=toupper(*j);

          // remove spaces
          while (myDataField.name.find(BLANK)!=string::npos)
            myDataField.name.erase(myDataField.name.find(BLANK), 1);

          // estimate data type of parameter value
          myDataField.type= estimate_data_type(s3.substr(pos_old, myDataField.width));

          // and insert the complete parameter
          myVec.insert(myVec.end(), myDataField);

          // reinitialize
          myDataField=DataField();
          ++pos;
        }
      }
    }
    return myVec;
  }


  void CSolarSystem::setAsteroidDataFile(const std::string& file_name)
  {
    FileStructure myStruc;
    ifstream input(file_name.c_str(), ios::in);
    if (input) {
      string rows[3];
      // read the first three lines from data file
      for(unsigned int i=0;i<(sizeof(rows)/sizeof(string));++i) getline(input, rows[i]);
      input.close();

      //try to analyze file structure
      myStruc = get_data_structure(rows[0], rows[1], rows[2]);

      if (myStruc.size() != 0) {
        // if the size of resulting data vector is not zero
        // we insert it into our file describing structure
        solarSystemData->bodyDataFiles.push_back(TBodyDataItem(file_name,myStruc));

#ifdef DEBUG
        cout << "CSolarSystem::setAsteroidDataFile(std::string&)"
             << " => successfully inserted asteroid data file ("
             << file_name << ")" << endl;
        cout.flush();
#endif //DEBUG

      }
      else
      {
		throw invalid_argument("CSolarSystem::setAsteroidDataFile(std::string&) => there has been an error during analyzing file ("+
			file_name + ")\n");
      }
    }
    else
    {
	  throw invalid_argument("CSolarSystem::setAsteroidDataFile(std::string&) =>asteroid data file (" + file_name + ") not available!\n");
    }
  }


  //---------------------------------------------------------------------------
  BodyID CSolarSystem::insertBody(CBody*body, bool force, int attempts)
  {

    pair<BodyMap::iterator, bool>retVal;

    BodyID startID=body->getId(); // ID from argument body data as first guess
    BodyID tempID=startID;

    do {
      // try to insert the given body into the solar system
      retVal=solarSystemData->bodies.insert(pair<int, CBody*>(tempID, body));

      ++tempID; // increase the ID; just for the case "force" is set to true
    } 
    while ((force==true)&&(retVal.second==false)&&((tempID-startID)<=attempts));

    if (retVal.second==true)
    {
      body->setId(retVal.first->first);
      body->setRefSolSys(this);
      std::cerr << "EXIT_TRACE registry: inserted key=" << retVal.first->first
                << " pointer=" << static_cast<const void*>(body) << std::endl;
      return retVal.first->first;
    } else return 0;
  } // BodyID CSolarSystem::insertBody(CBody*body, bool force, int attempts)


  //---------------------------------------------------------------------------
  void CSolarSystem::removeBody(const BodyID body)
  {
    BodyMap::const_iterator trace_body = solarSystemData->bodies.find(body);
    std::cerr << "EXIT_TRACE registry: remove key=" << body << " pointer="
              << (trace_body == solarSystemData->bodies.end() ? 0 : static_cast<const void*>(trace_body->second)) << std::endl;
    solarSystemData->bodies.erase(body);
  } // void CSolarSystem::removeBody(const BodyID body)


  //----------------------------------------------------------------------------
  const int CSolarSystem::size(void) const
  {
    return solarSystemData->bodies.size();
  } // const int CSolarSystem::size(void) const


  //---------------------------------------------------------------------------  
  struct CometFileDataSet {
    string numName;   ///< number and/or name of the comet
    string epoch;     ///< epoch of the elements
    string q;         ///< perihelion distance (AU)
    string e;         ///< eccentricity of the comets orbit
    string i;         ///< inclination of the orbit (deg) with respect to J2000 ecliptic plane
    string omega;     ///< argument of perihelion (J2000 - ecliptic)
    string Omega;     ///< longitude of ascending node (J2000 - ecliptic)
    string Tp;        ///< time of perihelion passage (YYYYMMDD.DDD)
    string reference; ///< orbit solution reference
  };


  //----------------------------------------------------------------------------
  BodyID CSolarSystem::insertComet(std::string stringID)
  {
    const int WIDTH_FIELD1=38; ///< name
    const int WIDTH_FIELD2=6;  ///< epoch
    const int WIDTH_FIELD3=11; ///< perihelion distance
    const int WIDTH_FIELD4=11; ///< eccentricity
    const int WIDTH_FIELD5=10; ///< inclination
    const int WIDTH_FIELD6=10; ///< argument of perihelion
    const int WIDTH_FIELD7=10; ///< longitude of ascending node
    const int WIDTH_FIELD8=15; ///< time of perihelion passage
    const int WIDTH_FIELD9=13; ///< data reference

    int result=0; // default value for an unsuccessful insertation
    char*idx;
    char buf[128];
    struct CometFileDataSet data_set;

    std::string::size_type i;

    RemoveBlanksUpcase(stringID);
    if (solarSystemData->cometFile.size() != 0)
    {
      ifstream input(solarSystemData->cometFile.c_str());
      if (input.is_open())
      {

        while (!input.eof() && data_set.numName.compare(stringID)!=0)
        {
          input.getline(buf, sizeof(buf));
          if (!input.eof())
          {
            data_set.numName.assign(buf, WIDTH_FIELD1);
            while ((i=data_set.numName.find(' '))!=string::npos)
            {
              data_set.numName.erase(i, 1);
            } // while ((i=data_set.numName.find(' '))!=string::npos)
          } // if (!input.eof())
        } // while (!input.eof() && data_set.numName.compare(stringID)!=0)

        if (!input.eof()) 
        {
          // parsing of data strings
          idx=buf+WIDTH_FIELD1;
          data_set.epoch.assign(idx, WIDTH_FIELD2);
          idx+=WIDTH_FIELD2;
          data_set.q.assign(idx, WIDTH_FIELD3);
          idx+=WIDTH_FIELD3;
          data_set.e.assign(idx, WIDTH_FIELD4);
          idx+=WIDTH_FIELD4;
          data_set.i.assign(idx, WIDTH_FIELD5);
          idx+=WIDTH_FIELD5;
          data_set.omega.assign(idx, WIDTH_FIELD6);
          idx+=WIDTH_FIELD6;
          data_set.Omega.assign(idx, WIDTH_FIELD7);
          idx+=WIDTH_FIELD7;
          data_set.Tp.assign(idx, WIDTH_FIELD8);
          idx+=WIDTH_FIELD8;
          data_set.reference.assign(idx, WIDTH_FIELD9);
        } // if (!input.eof()) 
        else return 0;
        input.close();
      } // if (input!=NULL)

      // remove possible spaces
      while ((i=data_set.Tp.find(' '))!=string::npos)
      { 
        data_set.Tp.erase(i, 1); 
      }

      // now the string variables have to be converted to double values

      double temp_epoch=0;
      double temp_q    =0;
      double temp_a    =0;
      double temp_e    =0;
      double temp_i    =0;
      double temp_omega=0;
      double temp_Omega=0;
      double temp_Tp   =0;
      int temp_year    =0;
      int temp_month   =0;
      int temp_day     =0;
      double temp_hour =0;

      istringstream(data_set.epoch)>>temp_epoch;
      istringstream(data_set.q)>>temp_q;
      istringstream(data_set.e)>>temp_e;
      istringstream(data_set.i)>>temp_i;
      temp_i*=DEG2RAD;
      istringstream(data_set.omega)>>temp_omega;
      temp_omega*=DEG2RAD;
      istringstream(data_set.Omega)>>temp_Omega;
      temp_Omega*=DEG2RAD;

      // convert the strings to double values
      istringstream(data_set.Tp.substr(0, 4))>>temp_year;
      istringstream(data_set.Tp.substr(4, 2))>>temp_month;
      istringstream(data_set.Tp.substr(6, 2))>>temp_day;
      istringstream(data_set.Tp.substr(9, 5))>>temp_hour;
      temp_hour*=24.0/100000.0;

      // and calculate a modified julian date out of it
      temp_Tp = date2MJD(temp_year, temp_month, temp_day, temp_hour);

      // create a new planet object
      CBody*temp_planet=0;

      temp_planet=new CBody;
      if (temp_planet!=0) {

        temp_a=OrbitParameter(temp_q, temp_e)/(1.0-temp_e*temp_e);

        double mean_motion=0;

        mean_motion=sqrt(getBody(SUN)->getGM(AUD)/pow(fabs(temp_a), 3.0));
        result=insertBody(temp_planet, true);

        if (result!=0) {
          temp_planet->Initialize(result, // BodyIndex (will be overwritten)
            stringID, 1.0, // default mass of 1 kg
            temp_Tp, // time of perihelion passage for epoch
            temp_a, // the calculated semimajor axis
            temp_e, // the given eccentricity
            temp_i, // inclination of the orbit
            temp_Omega, // longitude of ascending node
            temp_Omega+temp_omega, // curly pi
            temp_Omega+temp_omega, // mean longitude (perihelion passage)
            0, // da
            0, // de
            0, // di
            0, // dOmega
            0, // domega
            mean_motion, // dL
            ELEMENTS, // set state with Keplerian elements
            getBody(SUN) // pointer to center body
            );

          temp_planet->CState::SetState(
            Vector(temp_a, temp_e, temp_i, temp_Omega, temp_Omega+temp_omega,
            temp_Omega+temp_omega),
            ELEMENTS, temp_Tp, SUN, AUD, BODY, false, DE405_NO_BODY,
            getBody(getPlanetIndex(SUN_NAME)));
        } else
          Output(string("CSolarSystem::InsertComet Error with inserting planet class into data structure\n"));

      } else
        Output(string("CSolarSystem::InsertComet Error while allocating CBody object\n"));
    }
    return result;
  }


  /** Functional object used with STL's find - algorithm to obtain
      a certain data member in file structure description. The key
      which is tested to decide, whether it is the right element or
      not, is a data member named name (STL string). So, to use this
      functional object, it has to be assured that this particular
      data member is part of class T, which is the template parameter.
      Here it is used for class "Datafield" (typedef struct TDatafield Datafield). */
  template<class T>
  class name_equal : public unary_function<T, bool> {
    string _s;
    public:
    /** Tbd. */
    explicit name_equal(const string& s) : _s(s) {};
    /** Tbd. */
    bool operator()(T x) const {return (x.name.compare(_s)==0);};
  };

  //---------------------------------------------------------------------------
  BodyID CSolarSystem::insertAsteroid(std::string stringID) {

    BodyID result = 0;                     // default value (unsuccessful)
    string ws,ws2;                         // work strings

    RemoveBlanksUpcase(stringID);    // preprocess ID

#ifdef EXTRA_OUTPUT
    cout << "Looking for asteroid (" << stringID << ") in " << endl;
    cout.flush();
#endif

    istringstream iws;  // input work string
    bool found = false; // flag whether the desired line in data file has been found
    ifstream input;     // input file stream
    unsigned int data_files;  // index in vector of data file descriptors
    int lc;             // line counter in data file
    bool key_found;     // indicates, whether desired data column could have been
    // found. Usually this is the one named "NAME" or "DESCRIPTION"
    // It has to be determined one time for each data structure (data file) to identify
    // the body names in data file and compare them to given name of desired body.

    // number of stored asteroid data file names (and it's data describing structures)
    const unsigned int N_data_files = solarSystemData->bodyDataFiles.size();

    // loop over data files
    for (data_files = 0; data_files < N_data_files; ++data_files)
    {
      key_found=false;  // initialze flag to "desired body name not found"

      // get file describing structure of actual data file to process
      FileStructure* fileStruc =
        &solarSystemData->bodyDataFiles[data_files].second;

      // index in structure
      FileStructure::iterator j;

      FileStructure::iterator first = fileStruc->begin(); // first element of data structure
      FileStructure::iterator last  = fileStruc->end();   // and last (for loop purposes)

      // look, whether there's a data field named "NAME" or "DESIGNATION"
      // in the particular file structure description
      // case 1 (NAME)       : numbered asteroids data file
      // case 2 (DESIGNATION): unnumbered asteroids data file
      j = find_if(first, last, name_equal<DataField>("NAME"));

      if (j==last)
      {
        // the desired column name of data file could not been found. Then we have to try
        // another possible column heading which is "DESIGNATION".
        j = find_if(first, last, name_equal<DataField>("DESIGNATION"));

        // set flag to indicate, that we found the right column
        if (j!=last) key_found=true;
      } else key_found = true;

      if (key_found)
      {
        // If we found the right data field in file structure description,
        // we may open the data file with the stored filename.
        // But first we have to clear a possibly occured error flag in the
        // input file stream. This may happen only if another file has been
        // processed before.
        input.clear();

        // open data file
        input.open(solarSystemData->bodyDataFiles[data_files].first.c_str());

        if (input.is_open()) 
        {

#ifdef EXTRA_OUTPUT
          cout << solarSystemData->bodyDataFiles[data_files].first << endl;
          cout.flush();
#endif

          lc=0;                       // initialize line counter
          while (!input.eof()) 
          {

            getline(input, ws);++lc;  // read line and increase line counter
            if (!input.eof())
            if ( (lc > 2) && (ws.length() > 0) ) {

              ws2.assign(ws.substr((*j).pos, (*j).width)); // get name
              // iterator j still points to data element "NAME" or "DESIGNATION"
              // so we can use it here to obtain the actual asteroid name of
              // current line (string read from data file).
              RemoveBlanksUpcase(ws2);               // preprocess

              if (ws2.compare(stringID)==0) {         // I it the desired one?
                // Yes, it is, so we may read the information.
                found = true;

#ifdef EXTRA_OUTPUT
                cout << "Found asteroid!" << endl; cout.flush();
#endif

                double temp_epoch=0;
                double temp_a    =0;
                double temp_e    =0;
                double temp_i    =0;
                double temp_omega=0;
                double temp_Omega=0;
                double temp_M    =0;


                // epoch
                j = find_if(first, last, name_equal<DataField>("EPOCH"));
                if (j!=last) istringstream(ws.substr((*j).pos, (*j).width)) >> temp_epoch;

                // semi major axis
                j = find_if(first, last, name_equal<DataField>("A"));
                if (j!=last) istringstream(ws.substr((*j).pos, (*j).width)) >> temp_a;

                // eccentricity
                j = find_if(first, last, name_equal<DataField>("E"));
                if (j!=last) istringstream(ws.substr((*j).pos, (*j).width)) >> temp_e;

                // inclination
                j = find_if(first, last, name_equal<DataField>("I"));
                if (j!=last) istringstream(ws.substr((*j).pos, (*j).width)) >> temp_i;
                temp_i*=DEG2RAD;

                // argument of pericenter
                j = find_if(first, last, name_equal<DataField>("W"));
                if (j!=last) istringstream(ws.substr((*j).pos, (*j).width)) >> temp_omega;
                temp_omega*=DEG2RAD;

                // longitude of ascending node
                j = find_if(first, last, name_equal<DataField>("NODE"));
                if (j!=last) istringstream(ws.substr((*j).pos, (*j).width)) >> temp_Omega;
                temp_Omega*=DEG2RAD;

                // mean anomaly
                j = find_if(first, last, name_equal<DataField>("M"));
                if (j!=last) istringstream(ws.substr((*j).pos, (*j).width)) >> temp_M;
                temp_M*=DEG2RAD;

                CBody* temp_planet=0;
                temp_planet=new CBody;
                if (temp_planet!=0) {

                  double mean_motion=0;
                  mean_motion=sqrt(getBody(SUN)->getGM(AUD)/pow(temp_a, 3.0));
                  result=insertBody(temp_planet, true);

                  if (result!=0) {
                    temp_planet->Initialize(result, // BodyIndex (will be overwritten by the method
                                            stringID, 1.0, // default mass of 1 kg
                                            temp_epoch,    // epoch
                                            temp_a, temp_e, temp_i, temp_Omega, temp_Omega+temp_omega, // curly pi
                                            temp_Omega+temp_omega+temp_M, // mean longitude
                                            0, // da
                                            0, // de
                                            0, // di
                                            0, // dOmega
                                            0, // domega
                                            mean_motion, // dL
                                            ELEMENTS, getBody(SUN));

#ifdef EXTRA_OUTPUT
                    cout << "Initialized asteroid!" << endl; cout.flush();
#endif


                    temp_planet->CState::SetState(Vector(temp_a,
                                                         temp_e,
                                                         temp_i,
                                                         temp_Omega,
                                                         temp_Omega+temp_omega,
                                                         temp_Omega+temp_omega+temp_M),
                                                  ELEMENTS,
                                                  temp_epoch,
                                                  SUN,
                                                  AUD,
                                                  BODY,
                                                  false,
                                                  DE405_NO_BODY,
                                                  getBody(getPlanetIndex(SUN_NAME)));
#ifdef EXTRA_OUTPUT
  cout << "Set elements of asteroid!" << endl; cout.flush();
#endif

                  } else FatalError("CSolarSystem::InsertAsteroid; Error with inserting planet class into data structure");
                } else FatalError("CSolarSystem::InsertAsteroid; Error while allocating CBody object");
              }
            }
            if (found == true) break;
          }
        }
        input.close();
      }
      if (found == true) break;
    }
    return result;
  }


  // returns a pointer to the desired planet
  //---------------------------------------------------------------------------
  CBody* CSolarSystem::getBody(const BodyID bodyID) const {
    BodyMap::iterator pos;

    pos=solarSystemData->bodies.find(bodyID);
    if (pos!=solarSystemData->bodies.end())
      return pos->second; else
      return 0;
  }


  //---------------------------------------------------------------------------
  // returns the name of a planet stored in the CSolarSystem - class
  const char*CSolarSystem::Name(const BodyID planet) const {
    CBody*temp=getBody(planet);

    if (temp!=0) return temp->getName().c_str(); else
      return 0;
  };

  //----------------------------------------------------------------------------
  const std::string& CSolarSystem::getBodyName(const BodyID body_id) const
  {
    const CBody* body = getBody(body_id);
    if (body != 0)
    {
      return body->getName();
    }
    else
    {
      return S_NOT_APPLICABLE;
    }
  } // const std::string& getBodyName(const BodyID body) const


  //---------------------------------------------------------------------------
  // returns the number of planets actually stored in the solar system class
  int CSolarSystem::NumberOfPlanets(void)const {
    return solarSystemData->bodies.size();
  }


  Vector target_cart(6), center_cart(6);
  //----------------------------------------------------------------------------
  // returns a vector filled with the cartesian coordinates and velocity
  // componentes in a cartesian reference frame centered at the given "center"-
  // body at the given "date"
  Vector& CSolarSystem::Cart(const BodyID a_target, const BodyID a_center,
    const double date, const UnitType unit_type) 
  {

      CBody*target=0;
      CBody*center=0;

      target=getBody(a_target);      
      if (target==0) FatalError("CSolarSystem::Cart(const BodyID, \
 const BodyID, const double) => invalid first body index");
      target_cart = target->State(date).Cart(unit_type);

      center=getBody(a_center);
      if (target==0) FatalError("CSolarSystem::Cart(const BodyID, \
 const BodyID, const double) => invalid second body index");
      center_cart = center->State(date).Cart(unit_type);

      while (target!=target->getCenterState())
      {
        target=dynamic_cast<CBody*>(target->getCenterState());
        target_cart += target->State(date).Cart(unit_type);        
      } // while (target!=target->getCenterState())

      while (center!=center->getCenterState())
      {
        center=dynamic_cast<CBody*>(center->getCenterState());
        center_cart += center->State(date).Cart(unit_type);        
      }

      if (target!=center)
        FatalError("CSolarSystem::Cart(target, center, date, unit type)=>\
Resulting base object is not equal");

      return target_cart-=center_cart;
  } // Vector CSolarSystem::Cart(const BodyID, const BodyID, const double
  //, const UnitType unit_type)



  //---------------------------------------------------------------------------
  Vector CSolarSystem::Cart(const BodyID target, const BodyID center,
       const double date, const UnitType unit_type,
       const CoordSysType ref_frame, const double ref_frame_date) {

         // get the initial data in bodycentric, eclliptical
         // nonrotating reference frame (J2000) with the origin
         // at the center of the body that is indicated by the parameter "center"
         Vector result=Cart(target, center, date, unit_type);

         // define some variables, which will be used later
         Matrix T;
         Vector position;
         Vector velocity;
         double z_rot_angle_1;
         double x_rot_angle_1;
         double z_rot_angle_2;

         // do transformation according to the desired coodinate system
         switch (ref_frame) {
           case HEN:
             break;

           case HON:

             /* T = R_x(EPS0_J2000*ARCSEC2RAD);   // set up rotation matrix

             position = T * result.Slice(0,2); // do transformation
             velocity = T * result.Slice(3,5); // also for velocity

             return (result = Stack(position, velocity)); // return
             resulting values */
             cout<<"CSolarSystem::Cart(BodyID, BodyID, double, UnitType, CoordSysType) => invalid reference frame!"
                  <<endl;
             break;

           case BEN: // body centered, eclliptical, nonrotating
             break; // return, because the data mustn't be transformed

           case BER:
             // body centered, eclliptical, rotating with projection of position
             // vector onto reference plane

             position=result.Slice(0, 2);
             velocity=result.Slice(3, 5);

             // determine the rotation angle from the y- and x- coordinate of the position vector
             // because this is the projection of the position vector in the reference plane
             z_rot_angle_1=atan2(position(1), position(0));

             // rotate around the z-axis
             T=R_z(z_rot_angle_1);

             return (result=Stack(T*position, T*velocity));

           case BON:
             z_rot_angle_1=omega(center, ref_frame_date); // z-axis rotation 1
             // (arg. of periapsis)
             x_rot_angle_1=i(center, ref_frame_date); // x-axis rotation   (inclination)
             z_rot_angle_2=Omega(center, ref_frame_date); // z-axis rotation 2
             // (argument of line of nodes)

             T=R_z(z_rot_angle_1)*R_x(x_rot_angle_1)*R_z(z_rot_angle_2);

             position=T*result.Slice(0, 2);
             velocity=T*result.Slice(3, 5);

             return (result=Stack(position, velocity));

           case BORT:
             // !!!! true anomaly is of date, NOT of ref_frame_date !!!!
             z_rot_angle_1=omega(center, ref_frame_date)+f(center, date); //
             // z-axis rotation 1 (arg. of periapsis + true anomaly)
             x_rot_angle_1=i(center, ref_frame_date); // x-axis rotation   (inclination)
             z_rot_angle_2=Omega(center, ref_frame_date); // z-axis rotation 2
             // (argument of line of nodes)

             T=R_z(z_rot_angle_1)*R_x(x_rot_angle_1)*R_z(z_rot_angle_2);

             position=T*result.Slice(0, 2);
             velocity=T*result.Slice(3, 5);

             return (result=Stack(position, velocity));

           default:
             cout<<"CSolarSystem::Cart(BodyID, BodyID, double, UnitType, \
           CoordSysType) => CoordSysType not known!" <<endl;
             break;
         }

         return result;
  }



  //------------------------------------------------------------------------------
  // structs and functions do perform the gravity assist
  //------------------------------------------------------------------------------


  // data structure that holds some information for the differential
  // equation function f(..) and the stop - function
  //------------------------------------------------------------------------------
  /** Tbd. */
  typedef struct SAuxiliaryData
  {
    double gm;        ///< the gravitational parameter of the central body during the gravity assist manuever
    double stopValue; ///< the value that stops the numerical intregration
  } AuxiliaryData;

  //------------------------------------------------------------------------------
  // equation of motion for the 2 - body problem
  // vector equation has been transformed into first order diff. equation system
  //------------------------------------------------------------------------------
  /** Tbd. */
  void f(double t, const Vector&y, Vector&yp, void*_paux) {

    yp=Vector(y.getDimension());
    // every distance component of the six dimensional vector is multiplied with
    // a factor built by the negative gravitational parameter of the central planet
    // which is divided by the cubed norm of the extracted distance vector r
    double factor=-1.0*((AuxiliaryData*) _paux)->gm/(pow(Norm(y.Slice(0, 2)), 3));
    int i;

    for (i=0; i<3; i++) {
      yp(i)=y(i+3);
      yp(i+3)=y(i)*factor;
    }
  }


  //------------------------------------------------------------------------------
  // the function that stops the numerical integration, when the absolute value
  // of vector y is greater than the stop value provided by the data structure
  // the argument 'pA' points to
  //------------------------------------------------------------------------------
  /** Tbd. */
  bool sf(double t, Vector& y, void* pA)
  {
    if (((AuxiliaryData*) pA)->stopValue<=Norm(y.Slice(0, 2)))
      return true;
    return false;
  }


  //------------------------------------------------------------------------------
  // gravity assist maneuver check (sphere of influence vs. sphere of alert)
  //------------------------------------------------------------------------------
  CSolarSystem::GACResult CSolarSystem::GravityAssistPossible(const BodyID planet, //
  // index of center body for gravity assist
  CState&sc, // s/c state
  double&dT, // last time step size
  const double epsDist // maximum distance from sphere of influence
  // as a fraction of the sphere of influence radius
  ) {

    if (epsDist==0)
      return FAILURE_NO_EPS_DIST;
    // actual time of s/c
    const double t=sc.MJD();

    // sphere of influence
    const double SOI=CSolarSystem::RSOI(planet, t);

    // sphere of alert
    const double SOA=SOI*(1.0+sc.Velocity()*dT/SOI);

    // six dimensional state vector from target planet/body to s/c
    const Vector sv =sc.Cart()-CSolarSystem::Cart(planet, sc.MJD());
    const Vector pos=sv.Slice(0, 2);
    const Vector vel=sc.Cart().Slice(3, 5);

    // distance between s/c and target planet
    const double norm_pos=Norm(pos);

    // return if s/c is outside the sphere of alert
    if (norm_pos>SOA)
      return OUTSIDE_SOA;
    // return also if the s/c is inside the sphere of influence
    if ((SOI-norm_pos)>=0.0)
      return INSIDE_SOI;
    // create a test CState object to get the radiakl velocity with respect
    // to the target body
    CState test(sv, CARTESIAN, sc.MJD(), planet, AUD);

    // provide the still missing data to the test CState-object
    test.SetState(sv, CARTESIAN, sc.MJD(),
         planet, AUD, BODY, 0, DE405_NO_BODY,
         &CSolarSystem::State(planet, sc.MJD()));
    const double dr=test.dr(); // get the radial velocity of s/c with respect to the target planet

    // return if the s/c is not approaching the "target" - planet
    if (dr>0)
      return INSIDE_SOA_MOVING_AWAY;
    // return if we are not close enough to the SOI
    if (fabs(SOI-norm_pos)/SOI>=epsDist)
      return INSIDE_SOA_MOVING_TOWARDS;
    // if all checks failed then we must be at the sphere of influence
    // within the desired accuray
    return SUCCESS;
  } // CSolarSystem::GravityAssistPossible


  // this method tries to compute the end state (radius and velocity) by numerical
  // integration using the Dormand & Prince method RK8(7)-13M (1981)
  //------------------------------------------------------------------------------
  bool CSolarSystem::GravityAssistNumerical(const BodyID planet, // center body
  // in CSolarSystem for gravity assist
  const double end_time, // assumed maximum of time
  CState&dest // state of the s/c, that is doing the gravity assist
  ) {

    Vector y; // the vector function to integrate
    double t0           =0; // starting time
    double t            =t0; // independent variable of integration
    double t_end        =end_time; // initial end time
    double t_old        =t0; // temp variable to store the last end time
    double t_temp       =t0; // temporary time variable
    long maxCalls       =2500; // maximum calls of f
    double relAcc       =1E-8; // relative accuracy
    double absAcc       =1E-10; // absolute accuracy
    const double distEps=1E-9; // the maximal absolute difference between
    // the desired and computed end distance
    double diff; // the actual difference between the norm of the
    // computed distance vector and the desired distance
    double t_inc; // increment/decrement for calculation of new end time
    AuxiliaryData AuxDat; // data structure for passing variables to f and sf
    int result; // result value of integration step (see 'MLDE.h')
    CState center_state; // the center of s/c at the beginning of maneuver
    Vector input; // vector at the beginning of integration
    Vector center; // s/c's center of motion during manuever
    Vector output; // final vector

    center_state=CSolarSystem::State(planet, dest.MJD()); // actual center of s/c coordinate
    AuxDat.gm=center_state.getGM(AUD); // get the  gravitational parameter of it
    // is the center of both bodies the same (usually the sun) ?
    if (center_state.getCenterBodyId()!=dest.getCenterBodyId())
      return false;
    // initialize position and velocity of the destination state for the epoch of the
    // beginning of gravity assist maneuver
    input=dest.Cart();
    center=CSolarSystem::Cart(planet, dest.MJD());
    // transform into local coordinate system for maneuver
    input=input-center;
    DOPRI8 solver(SpaceLib::f, 6, &AuxDat); // init numerical integrator

    solver.Init(maxCalls, relAcc, absAcc, sf); // additional initialisation (stop function)
    // for gravity assist maneuvers the end distance is equal to the starting distance
    // store this value in the data block that is given to f and the stop function
    AuxDat.stopValue=Norm(input.Slice(0, 2));
    do {
      // initial 6-dimensional vector y is built by the given source vectors
      // of position and velocity
      y=input;
      // integration of ODES (equations of motion)
      result=solver.Step(t, y, t_end-t0);

      diff=Norm(y.Slice(0, 2))-AuxDat.stopValue; // actual distance error ?
      t_inc=fabs(t_old-t)/2; // calculate new time step size
      if (diff<0) { // depending on the sign of the difference,
        t_temp=t0+t+t_inc; // increase end time
      } else {
        t_temp=t0+t-t_inc; // or decrease end time
      }

      t_old=t; // make new setup
      t=t0;
      t_end=t_temp;
    } while ((result==RK_STOP_FUNCTION)||(fabs(diff)>distEps));
    center=CSolarSystem::Cart(planet, dest.MJD()+t_old);
    output=y+center;
    // set new s/c state
    dest.SetState(output, CARTESIAN, dest.MJD()+t_old
      , dest.getCenterBodyId(), dest.Unit(), dest.getType()
      , 0, DE405_NO_BODY, dest.getCenterState());
    return true;
  }


  // this method calculates a gravity assist maneuver analytically
  // reference: Ian Canelli (2003)
  //---------------------------------------------------------------------------
  bool CSolarSystem::GravityAssistAnalytical(const BodyID planet, CState&dest,
       double i, double j) {
         /*ofstream matlab; matlab.open("matlab.csv", ios::app);
         ofstream report; report.open("report.txt", ios::app);
         report << setiosflags(ios::fixed);
         Vector Rad1(dest.Cart().Slice(0,2));
         Vector Vel1(dest.Cart().Slice(3,5)); report << "Angle(R1,V1)=" <<
         acos(Dot(Rad1,Vel1)/(Norm(Rad1)*Norm(Vel1)))*180.0/PI << " ";*/

         const double Rsph=RSOI(planet, dest.MJD());
         const double rel_planet_radius=0.0005/Rsph;

         // i = min(i*i*i,0.8); j = min(j*j*j,0.8);
         double delta=sqrt(i*i+j*j);
         if (delta>1.0) // delta is larger than the sphere of influence
         {
           i/=delta;
           j/=delta;
         }
         i*=Rsph; j*=Rsph;
         if (i==0.0&&j==0.0) // only for the case that both i and j are zero
         {
           i=rel_planet_radius;
           delta=i;
         }
         if (delta<rel_planet_radius) // delta is smaller than the planet radius
         {
           i/=delta/rel_planet_radius;
           j/=delta/rel_planet_radius;
         }
         //matlab << i << ", " << j << ", ";
         //report << "i=" << i/Rsph << " " << "j=" << j/Rsph << endl;
         delta=sqrt(i*i+j*j);
         if (delta>Rsph) delta=Rsph;
         const double k=sqrt(Rsph*Rsph-delta*delta);

         Vector r1_bplane(i, j, k);

         //report << "r1_bplane=" << r1_bplane << endl;

         const double vx=dest.dx()-dx(planet, dest.MJD());
         const double vy=dest.dy()-dy(planet, dest.MJD());
         const double vz=dest.dz()-dz(planet, dest.MJD());
         const Vector v1_cart(vx, vy, vz);

         //report << "v1_cart=" << v1_cart << endl;

         const double pl_phi  =phi(planet, dest.MJD());
         const double pl_theta=theta(planet, dest.MJD());
         const double sp      =sin(pl_phi);
         const double cp      =cos(pl_phi);
         const double st      =sin(pl_theta);
         const double ct      =cos(pl_theta);

         Matrix A(3, 3);

         A.SetRow(0, Vector(cp*ct, sp*ct, st));
         A.SetRow(1, Vector(-sp, cp, 0.0));
         A.SetRow(2, Vector(-cp*st, -sp*st, ct));

         Vector v1_sph=A*v1_cart;

         Vector ezb=-v1_sph/Norm(v1_sph);

         //report << "ezb=" << ezb << endl;

         const double elev   =atan2(ezb(2), sqrt(ezb(0)*ezb(0)+ezb(1)*ezb(1)));
         const double azim   =atan2(ezb(1), ezb(0));
         const double co_elev=HALFPI-elev;
         const double co_azim=HALFPI-azim;

         Vector r1_sph(R_z(co_azim)*R_x(co_elev)*r1_bplane);

         //report << "r1_sph=" << r1_sph << endl;

         const double Rr    =r1_sph(0)+r(planet, dest.MJD());
         const double Rphi  =r1_sph(1);
         const double Rtheta=r1_sph(2);

         Matrix InvA=Inv(A);

         Vector R1_cart=InvA*Vector(Rr, Rphi, Rtheta);

         dest.SetState(R1_cart(0), R1_cart(1), R1_cart(2), dest.dx(), dest.dy(), dest.dz(),
              CARTESIAN, dest.MJD());

         /*report << resetiosflags(ios::fixed); report.close();
         matlab.close();*/

         return GravityAssistAnalytical(planet, dest);
  };


  // this method calculates a gravity assist maneuver analytically
  // reference: Ian Canelli (2003)
  //------------------------------------------------------------------------------
  bool CSolarSystem::GravityAssistAnalytical(const BodyID planet, CState&dest) {
    /*ofstream matlab; matlab.open("matlab.csv", ios::app); ofstream report;
    report.open("report.txt", ios::app); matlab << dest.Cart() << ", ";
    matlab << Cart(planet, dest.MJD()) << ", ";*/

    if (CSolarSystem::State(planet, dest.MJD()).Unit()==dest.Unit()) {

      const double t0=dest.MJD(); // initial time t0 (beginning of maneuver)

      Vector R1pl(CSolarSystem::x(planet, t0), // planets heliocentric position
      CSolarSystem::y(planet, t0), CSolarSystem::z(planet, t0));

      Vector V1pl(CSolarSystem::dx(planet, t0), // planets heliocentric velocity
      CSolarSystem::dy(planet, t0), CSolarSystem::dz(planet, t0));

      Vector R1(dest.x(), dest.y(), dest.z()); // s/c heliocentric position

      Vector V1(dest.dx(), dest.dy(), dest.dz()); // s/c heliocentric velocity

      Vector r1=R1-R1pl;

      Vector v1=V1-V1pl;

      const double k=CSolarSystem::State(planet, t0).getGM(AUD); //
      // gravitational parameter of central body
      const double Rsph=CSolarSystem::RSOI(planet, t0); // radius of the sphere of influence

      /*matlab << k << ", " << Rsph << ", "; report << "Rsph=" << Rsph << " ";
      report << "Angle(R1,V1)=" << acos(Dot(R1,V1)/(Norm(R1)*Norm(V1)))*180.0/PI
      << " "; report << "Angle(r1,v1)=" <<
      acos(Dot(r1,v1)/(Norm(r1)*Norm(v1)))*180.0/PI << " ";*/
      const double delta=Norm(r1)*sin(acos((Dot(r1, v1)/(Norm(r1)*Norm(v1)))));
      //report << "Norm(v1)=" << Norm(v1) << " Norm(r1)=" << Norm(r1) << " ";
      const double vinf=sqrt(Dot(v1, v1)-2.0*k/Norm(r1));
      //report << "vinf=" << vinf*AU/DAY2SEC/1000.0 << " ";
      const double Phi    =2.0*atan(k/(delta*vinf*vinf));
      const double PhiV   =-Phi;
      const double Gamma  =asin(delta/Rsph);
      const double PhiStar=PI-Phi+2.0*Gamma;
      const double Phi0   =atan2(r1(1), r1(0));
      const double Theta0 =atan2(r1(2), sqrt(r1(0)*r1(0)+r1(1)*r1(1)));

      const Matrix Ry=R_y(-Theta0);
      const Matrix Rz=R_z(Phi0);
      const Matrix CP=Ry*Rz;

      const Vector v1polar=CP*v1;
      //report << "v1polar=" << v1polar << " ";
      const double Xsi=atan2(v1polar(2), v1polar(1));

      const Matrix Rr=R_x(Xsi); // rotation about the radial vector
      const Matrix R3=R_z(PhiV); // rotation of the velocity vector

      // rotation of the position vector about the axis of angular momentum
      const Matrix R3Star=R_z(PhiStar);

      // complete transformation matrices
      const Matrix OmegaV=Transp(Rz)*Transp(Ry)*Transp(Rr)*R3*Rr*Ry*Rz;
      const Matrix OmegaR=Transp(Rz)*Transp(Ry)*Transp(Rr)*R3Star*Rr*Ry*Rz;

      const Vector r2=OmegaR*r1;
      const Vector v2=OmegaV*v1;

      const double a     =1.0/(Dot(v1, v1)/k-2.0/Rsph);
      const double H     =MathLib::acosh((1.0+Rsph/a)*sin(Phi/2.0));
      const double deltaT=2.0*sqrt(pow(a, 3)/k)*(1.0/sin(Phi/2)*sinh(H)-H); //
      // gravity assist duration

      const Vector R2pl(CSolarSystem::x(planet, t0+deltaT), // planets heliocentric position
      CSolarSystem::y(planet, t0+deltaT), CSolarSystem::z(planet, t0+deltaT));

      const Vector V2pl(CSolarSystem::dx(planet, t0+deltaT), // planets heliocentric velocity
      CSolarSystem::dy(planet, t0+deltaT),
           CSolarSystem::dz(planet, t0+deltaT));

      /*report << "Angle(v1,v2)=" << acos(Dot(v1,v2)/(Norm(v1)*Norm(v2)))*180.0/PI
      << " "; report << "Angle(r2,v2)=" <<
      acos(Dot(r2,v2)/(Norm(r2)*Norm(v2)))*180.0/PI << " ";
      report << "Angle(r1,v2)=" <<
      acos(Dot(r1,v2)/(Norm(r1)*Norm(v2)))*180.0/PI << " ";*/

      const Vector R2=R2pl+r2; // transformation back to heliocentric coordinates
      const Vector V2=V2pl+v2;

      //report << "Angle(V2,V1)=" <<
      // acos(Dot(V2,V1+V1pl)/(Norm(V2)*Norm(V1+V1pl)))*180.0/PI << " ";
      //const double olda = dest.a();

      dest.SetState(Stack(R2, V2), //our final heliocentric position and velocity
      CARTESIAN, t0+deltaT, dest.getCenterBodyId(), dest.Unit(), dest.getType(),
           0, DE405_NO_BODY, dest.getCenterState());

      /*matlab << dest.Cart() << ", "; matlab << Cart(planet, dest.MJD()) <<
      endl; report << "DeltaV=" << (Norm(V2)-Norm(V1))*AU/DAY2SEC/1000.0 << "
      "; report << "a1=" << olda << " " << "a2=" << dest.a() << endl;

      report << resetiosflags(ios::fixed); report.close();*/
      return true;
    } else {
      /* report << resetiosflags(ios::fixed); report.close();
      matlab.close(); */
      return false;
    }
  }


  //---------------------------------------------------------------------------
  double CSolarSystem::Distance(const BodyID planet,
       const StateType target_state, const double a, const double e,
       const double i, const double Omega, const double omega,
       const CState&body) {
         double distance;

         if (target_state==FREE_ORBIT_STATE) {
           CState virt_state(a, e, i, Omega, omega,
                body.phi()-Omega-omega, KEPLERF);
           const Vector
                planet_r(virt_state.x(), virt_state.y(), virt_state.z());
           const Vector body_r(body.x(), body.y(), body.z());

           distance=Norm(planet_r-body_r);
         } else { distance=getBody(planet)->Distance(body, target_state); }
         return distance;
  };


  //---------------------------------------------------------------------------
  double CSolarSystem::Distance_km(const BodyID planet,
       const StateType target_state, const double a, const double e,
       const double i, const double Omega, const double omega,
       const CState&body) {

         return Distance(planet, target_state, a, e, i, Omega,
              omega, body)*AU2KM;
  };


  //---------------------------------------------------------------------------
  double CSolarSystem::Velocity(const BodyID planet,
       const StateType target_state, const double a, const double e,
       const double i, const double Omega, const double omega,
       const CState&body) {
         double velocity;

         if (target_state==FREE_ORBIT_STATE) {
           CState virt_state(a, e, i, Omega, omega,
                body.phi()-Omega-omega, KEPLERF);
           const Vector
                planet_v(virt_state.dx(), virt_state.dy(), virt_state.dz());
           const Vector body_v(body.dx(), body.dy(), body.dz());

           velocity=Norm(planet_v-body_v);
         } else { velocity=getBody(planet)->Velocity(body, target_state); }
         return velocity;
  };



  //---------------------------------------------------------------------------
  double CSolarSystem::Velocity_kms(const BodyID planet,
       const StateType target_state, const double a, const double e,
       const double i, const double Omega, const double omega,
       const CState&body) {

         return Velocity(planet, target_state, a, e, i, Omega, omega,
              body)*AUD2KMS;
  };


  //----------------------------------------------------------------------------
  BodyID CSolarSystem::getPlanetIndex(const string& name) const 
  {
    BodyMap::iterator entry;
    entry = solarSystemData->bodies.begin();
    std::string _name = name;
    upper(_name);
    while (entry != solarSystemData->bodies.end())
    {
      if (entry->second->getName() == _name) return entry->first;
      ++entry;
    } // while (entry != solarSystemData->bodies.end())
    return 0;
  } // BodyID CSolarSystem::getPlanetIndex(const string&planet_name) const 

} // namespace SpaceLib

// end of file
