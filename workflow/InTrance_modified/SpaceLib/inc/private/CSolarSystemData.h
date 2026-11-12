// file     : CSolarSystemData.h
// author   : Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
// compiler : GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++
// date     : 16.08.2006 (AO)
// purpose  : Definition of the data structure where the data of
//            CSolarSystem class is stored in.

#ifndef CSOLARSYSTEMDATA_H
#define CSOLARSYSTEMDATA_H

//SpaceLib
#include <configSpaceLib.h>

//STL
#include <cstddef>
#include <map>
#include <vector>
#include <string>

#include <CState.h>
#include <DE405.h>

namespace SpaceLib {

  // some constants
  const char* DEFAULT_NAME = "NONAME";
  const char BLANK = ' ';


  // data structure item
  // This data field structure is used in conjunction with the parsing
  // of data files in text format like the asteroid and comet information
  // data files provided by NASA (http://ssd.jpl.nasa.gov/sb_elem.html).
  // Those files are constructed of fixed width columns and the first
  // row contains the columns headers whereas row two is used to determine
  // the width of the data columns
  // example:
  //
  //
  // 1  Num     Name       Epoch     a           e         i
  // 2------- ------------ ----- ---------- ---------- ---------
  // 31       Ceres        52400  2.7664122 0.07911582  10.58348
  //
  //--------------------------------------------------------------------------
  struct TDataField {

    enum DataType { TEXT=0, FLOAT=1, INTEGER=2, NOTYPE=3 };


    // data fields
    std::string name;
    int pos;
    int width;
    DataType type;

    // destructor
    ~TDataField(void) { };

    // standard constructor
    TDataField(void):name(DEFAULT_NAME), pos(0), width(0), type(NOTYPE) { };

    // copy constructor
    TDataField(const struct TDataField&toCopy):name(toCopy.name), pos(toCopy.pos),
      width(toCopy.width), type(toCopy.type) { };

    // assign operator
    TDataField& operator=(const struct TDataField& toAssign) {
      name=toAssign.name;
      pos=toAssign.pos;
      width=toAssign.width;
      type=toAssign.type;
      return (*this);
    };

    // compare operator (equal)
    bool operator==(const struct TDataField& toTest) {
      if ((name==toTest.name)   &&
        (pos==toTest.pos)     &&
        (width==toTest.width) &&
        (type==toTest.type))       return true;
      return false;
    }

    // compare operator (less)
    bool operator<(const struct TDataField& toTest) {
      if (name<toTest.name) return true;
      return false;
    }

  };

  // stream operator to be able to make some debug output
  ostream& operator<<(std::ostream& to, const struct TDataField& from) {
    to<<from.name.c_str()<<"," <<from.width<<",";
    switch (from.type) {
        case TDataField::TEXT:
          to<<"TEXT";
          break;
        case TDataField::FLOAT:
          to<<"FLOAT";
          break;
        case TDataField::INTEGER:
          to<<"INTEGER";
          break;
        case TDataField::NOTYPE:
          to<<"NOTYPE";
          break;
    }
    return to;
  }

  typedef struct TDataField DataField;

  // logical file structure (subsequent data fields)
  typedef std::vector<struct TDataField> FileStructure;

  // a pair structure to hold the file name and file structure
  typedef std::pair<string, FileStructure> TBodyDataItem;

  // a vector of TBodyDataItems (to handle more than one file)
  typedef std::vector<TBodyDataItem> TBodyDataFiles;

  // type of data structure to store solar system bodies
  typedef std::map<int,CBody*> BodyMap;


  struct CSolarSystemData 
  {
    std::string DE405File;           // filename of JPL's DE405 data source file
    TBodyDataFiles bodyDataFiles;              // data files (asteroids, comets)
    std::string cometFile;
    DE405* myDE;  ///< DE405 object pointer to get the planets state information
    BodyMap bodies;     ///< STL map of CBody objects incorporating space bodies

    CSolarSystemData(void);                            ///< standard constructor
    ~CSolarSystemData(void);                                     ///< destructor

  }; // struct CSolarSystemData


  //----------------------------------------------------------------------------
  inline CSolarSystemData::CSolarSystemData(void)
    : DE405File(), bodyDataFiles(), cometFile(), bodies()
  {
    myDE = new DE405();
  } //CSolarSystemData(void)

  //----------------------------------------------------------------------------
  inline CSolarSystemData::~CSolarSystemData(void)
  {
    while (!bodies.empty())       // loop until map of celestial bodies is empty
    {
      std::cerr << "EXIT_TRACE solar system: deleting body key=" << bodies.begin()->first
                << " pointer=" << static_cast<const void*>(bodies.begin()->second) << std::endl;
      delete ((bodies.begin())->second);                   // delete body object
      std::cerr << "EXIT_TRACE solar system: body deleted" << std::endl;
      bodies.erase(bodies.begin());                        // and it's map entry
    }
    std::cerr << "EXIT_TRACE solar system: before DE405 delete" << std::endl;
    SafeDelete(myDE);                 // delete DE405 ephemeris interface object
    std::cerr << "EXIT_TRACE solar system: destructor body complete" << std::endl;
  } // ~CSolarSystemData(void)

} //  namespace SpaceLib

#endif // CSOLARSYSTEMDATA_H
