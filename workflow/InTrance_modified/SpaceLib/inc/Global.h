/*!
\file Global.h
\version 1.7
\brief Some global functions that are widely used throughout the SpaceLib.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 08.06.2006 (AO)

\b Purpose:
   This file provides the interface for some global functions that are used
   within the SpaceLib and that can be used in programs using the SpaceLib.

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

#ifndef GLOBAL_H
#define GLOBAL_H

// STL
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <functional>
#include <stdexcept>
#include <valarray>

// SpaceLib

#include <CPhysUnit.h>

#pragma warning (disable:4786) 
#pragma warning (disable:4150) // deletion of non complete object pointers; encountered with template "SafeDel"
#pragma warning (disable:4290) // ignored exception declaration

#define development_version

// typedefs
typedef unsigned int uint;

const int KEY_ESCAPE =  27;
const int KEY_PLUS   =  43;
const int KEY_MINUS  =  45;
const int KEY_A      =  97;
const int KEY_B      =  98;
const int KEY_E      = 101;
const int KEY_H      = 104;
const int KEY_I      = 105;
const int KEY_M      = 109;
const int KEY_P      = 112;
const int KEY_S      = 115;

/// Some global functions that are used within the SpaceLib
/// and that can be used in programs using the SpaceLib.
namespace Global
{

  /** 
    \brief Returns the greater one of two given values.
    \param a first value
    \param b second value
    \return the greater of the two given values
    \sa min 
  */
  //template<class T> inline T max(const T lhs,T rhs) {return lhs<rhs?rhs:lhs;};

  /** 
    \brief Returns the smaller one of two given values.
    \param a first value
    \param b second value
    \return the smaller of the two given values
    \sa max 
  */
  //template<class T> inline T min(const T lhs,T rhs) {return lhs>rhs?rhs:lhs;};

  /// text ouput base class
  class SLOutput
  {

  private:

    std::string log_filename;

  public:

    /// constructor
    SLOutput(void) {};

    /// destructor
    virtual ~SLOutput(void) {};

    /// get log messages filename
    const char* getLogFilename(void) const
    {
      return log_filename.c_str();
    };

    /// set log messages filename
    void setLogFilename(const std::string& filename)
    {
      log_filename = filename;
    }

    // pure virtual output streaming function (char pointer argument)
    virtual SLOutput& operator<<(const char* message);

    // pure virtual output streaming function (string argument)
    virtual SLOutput& operator<<(const std::string& message)
    {
      return ((*this) << message.c_str());
    }

    /// assign operator
    SLOutput& operator=(const SLOutput& rhs)
    {
      log_filename = rhs.log_filename;
    }
  };

  /// set filename of global log messages file
  void setOutputLogFilename(const std::string filename);

  /// return used global log messages filename
  const char* getOutputLogFilename(void);


  // Derived standard console output class.
  class SLStandardOutput: virtual public SLOutput
  {
  public:

    /// constructor.
    SLStandardOutput(void):SLOutput() {};

    /// destructor.
    virtual ~SLStandardOutput(void) {};

    /// const char* stream operator.
    virtual SLOutput& operator<<(const char* message)
    {
      if (message != 0)
      {
        std::cout << message;
        return (SLOutput::operator <<(message));
      }
      return *this;
    }
  };


  /// global output function (string argument).
  void Output(const std::string& message);

  /// global output function (const char* argument).
  void Output(const char* message);

  /// global output function (output stringstream object).
  inline void Output(const std::ostringstream& message)
  {
    Output(message.str());
  }


  /*! \brief global output function (std::valarray<>)
   *
   * \param d std::valarray to put out
   */
  template<class T>
    inline void Output(const std::valarray< T>& d)
  {
    if (d.size()>0)
    {
      std::ostringstream s;
      for (uint i=0;i<d.size()-1;++i) s << d[i]<<", ";
      s << d[i] << std::endl;
      Output(s.str());
    }
  }


  /*! \brief global output function (std::vector<>)
   *
   * \param d vector to put out
   */
  template<class T>
    inline void Output(const std::vector< T>& d)
  {
    if (d.size()>0)
    {
      std::ostringstream s;
      for (uint i=0;i<d.size()-1;++i) s << d[i]<<", ";
      s << d[i] << std::endl;
      Output(s.str());
    }
  }


  /// Sets the output object; all further output data
  /// which is being transfered through Output(std::string)
  /// goes to this object
  void setOutput(SLOutput* output);

  /// Reset output to standard C++ console output.
  void restoreStandardOutput(void);


  // file reading functions

  /// Reads a line from the input file stream (read until endline).
  std::string ReadString(std::ifstream& in);

  /// Reads a number (comma or not) from the input file stream in a double value (read until endline).
  double ReadNumber(std::ifstream& in);

  // error handling functions

  /// Streams an error message to the error stream.
  inline void Warning(const std::string& warning_message)
  {
    std::cerr << "\nWARNING: " << warning_message << "\n\n";
  }

  /// Streams an error message to the error stream and exits the program with error code 1-
  inline void FatalError(const std::string& error_message)
  {
    flush(std::cerr << "\nERROR: " << error_message << "\n\n");
    exit(1);
  }

  /// Checks whether the given input file exists and calls FatalError if not.
  inline void Assure(std::ifstream& file, const char* filename)
  {
    if (!file)
      FatalError("could not open file " + std::string(filename));
  }

  /// Checks whether the given output file exists and calls FatalError if not.
  inline void Assure(std::ofstream& file, const char* filename)
  {
    if (!file)
      FatalError("could not open file " + std::string(filename));
  }

  /// Checks whether the given memory address is not NULL and calls FatalError if it is.
  inline void Assure(void* memory, const std::string& error_function)
  {
    if (!memory)
      FatalError(error_function + ": memory allocation failed");
  }

  /*! \brief check whether file is available/accessible (const char*)
   * \param name name of file that should be checked for availability
   * \return true if file exists, false if not
   */
  inline bool InputFileExistent(const char* name)
  {
    std::ifstream f(name);
    if (f.is_open()==true) return true;
    else return false;
  }

  /*! \brief check whether file is available/accessible (const std::string)
   * \param name name of file that should be checked for availability
   * \return true if file exists, false if not
   * \sa InputFileExistent(const char*)
   */
  inline bool InputFileExistent(const std::string& name)
  {
    return InputFileExistent(name.c_str());
  }


  /// converts a string to upcase
  std::string& upper(std::string& s);

  /// removes all occurences of the given character
  void removeChar(std::string& s, const char character);

  /// deletes all blank characters within given string
  inline std::string& removeBlanks(std::string& s)
  {
    removeChar(s, ' ');
    return s;
  }

    /** Like the function name already says, it removes space characters (BLANKS)
      possibly contained in the given string and converts all string letters
      to upcase. It is just a support function to improve the readability of
      CSolarSystem::insertAsteroid. */
  inline std::string& RemoveBlanksUpcase(std::string& s)
  {
    return upper(removeBlanks(s));
    //string::size_type i = 0;
    //while ((i=s.find(BLANK))!=string::npos) s.erase(i,1);
    //for (i=0;i<s.size();++i) s[i]=toupper(s[i]);
    //return s;
  } // void remove_blanks_and_upcase(string& s)

  // removes all occurences of the given characters
  void removeChars(std::string& s, const std::string& characters);

  /// removes all characters after given seperator string until end of line
  void removeComments(std::string& s, const std::string& sep);

  /// removes all characters from string beginning at the end until
  /// a character is not equal one of the given string characters
  void removeTrailingChars(std::string& s, const std::string& trail_chars);


  /** The base class for all Spacelib exceptions. It is able to put out an error message.
      In general, it is only the base class for more types of exceptions in the SpaceLib. */
  class SLEXC : public std::exception
  {
  private:

    std::string message; ///< exception message

  public:

    /** \brief Standard constructor. */
    SLEXC(void) : message() {};

    /** \brief Copy constructor. */
    SLEXC(const SLEXC& exc) : message(exc.message) {};

    /** \brief String constructor. */
    explicit SLEXC(const std::string& msg) :message(msg) {};

    /** \brief Print exception message. */
    virtual void PrintEXC(void) throw()
	  {
  	  Output(message);
  	};

    /// display exception message
    virtual const char* what() const throw() { return message.c_str(); };

    /** \brief Destructor. */
    virtual ~SLEXC(void) throw() {};

    /** \brief Assign operator. */
    SLEXC& operator=(const SLEXC& exc) {
      message = exc.message;
      return *this;
    };

  }; // end class SLEXEC

  /// seperate path name from fully qualified file name
  void StripPathFromFilename(std::string& name);

  /*! \brief Parameter item. 
   *
   * First string of this pair is the parameter name and the second one the
   * value (also expressed as string).
   */
  typedef std::pair<std::string, std::string> Param;

  /*! \brief List of parameters. */
  typedef std::vector<Param> ParamList;



  /** Functional object used with STL's find-algorithm to obtain a certain data
      member in file structure description. The key, which is tested to decide,
      whether it is the right element or not, is a data member named name (STL
      string). So, to use this functional object, it has to be assured that this
      particular data member is part of class T, which is the template parameter.
      Here it is used for class "Datafield" (typedef struct TDatafield Datafield). */
  template<class T>
  class name_equal : public std::unary_function<T, bool>
  {
    std::string _s; ///< tbd
  public:
    /** Tbd */
    explicit name_equal(std::string s): _s(s) {};
    /** Tbd */
    bool operator()(T x) const { return x.first==_s; };
  };

  /// Exception that is thrown in the case of a syntax error.
  class syntax_error_in_config_file : public std::exception
  {
    std::string _msg; ///< tbd
  public:
    /** Tbd */
    explicit syntax_error_in_config_file(const std::string& msg):_msg(msg) {};
    /** Tbd */
    virtual ~syntax_error_in_config_file() throw() {};
    /** Tbd */
    virtual const char* what() const throw() { return _msg.c_str(); };
  };


  /*!
   *\brief Reads parameter from given input file stream.
   *
   * All lines of the given input file stream will be read and
   * processed as follows:\
   * 1. remove all characters after a possibly given comment seperator
   * (inclusive comment seperator); such a comment
   * seperator may be '#' or '//'\n
   * 2. remove possible spaces or tab stops between the end of parameter value 
   * and removed comments\n
   *
   * The parameters of input file must be in the follwing syntax:\n
   * paramter name=parameter value possible comment
   *\param ifile The input file stream to read parameters from.
   *\param comment_seperator A string type identifier from where on
   *all following characters (until occurence of a new line) are  
   *treated as comments.
   *\exception syntax_error_in_config_file The function could not find
   *a assign character ('=') and assumes that the given input file stream
   *is not a valid configuration file.
   *\sa InputParameter,InputStringParameter
   */
  ParamList GetParameters(std::ifstream& ifile, const std::string& comment_seperator) 
    throw (syntax_error_in_config_file);



  /*! \brief extract and remove unit string from given string containing a number
  \return string with everything from the end of given string not belonging
    to a number*/
  std::string RemoveUnitString(std::string& s);


  /*!
   *\brief Assigns parameter value from parameter list to given reference.
   *\param name parameter name.
   *\param dest reference of data type to hold the parameter value
   *\param params parameters storage object
   *\param exception_text_prefix message which will be displayed before
   * actual exeception message (for tracing purposes)
   *\param trailing_chars string of characters which should be removed after
   * a data value (eg. tabs, blanks)
   *\return unit object of parsed (or not) physical unit
   *\exception std::invalid_argument The parameter "name" could not be found in
   *given paramter list.
   *\sa 
   */
  template<class T>
  inline SpaceLib::CPhysUnit* InputParameter(
    const std::string& name
   ,      T&           dest
   ,      ParamList&   p
   ,const std::string& exception_text_prefix
   ,const std::string& trailing_chars
   ) throw(std::invalid_argument)
  {
    ParamList::iterator i = find_if(p.begin(),p.end(),name_equal<Param>(name));
    if (i!=p.end()) 
    {
      // remove not needed chars
      removeTrailingChars((*i).second, trailing_chars); 
      removeChar((*i).second, '\"');                   // remove possible quotes
      std::string s_unit(RemoveUnitString((*i).second));
      istringstream((*i).second) >> dest;    // use stream operator on reference
      p.erase(i);                          // delete read parameter from storage
      return SpaceLib::GetUnit(s_unit);
    }
    else throw invalid_argument(
      exception_text_prefix+string(": missing parameter (")+name+string(")\n"));
    return SpaceLib::U_NO_UNIT;
  };


  /*!
   *\brief Assigns parameter value from parameter list to given reference.
   *\param name parameter name.
   *\param dest reference of data type to hold the parameter value
   *\param params parameters storage object
   *\param exception_text_prefix message which will be displayed before
   * actual exeception message (for tracing purposes)
   *\param trailing_chars string of characters which should be removed after
   * a data value (eg. tabs, blanks)
   * \param check_key_found throws exception in case of not found parameter key if
   * set to true, otherwise not
   *\return unit object of parsed (or not) physical unit
   *\exception std::invalid_argument The parameter "name" could not be found in
   *given paramter list (check_key_found must be set to true).
   *\sa InputStringParameter, InputBoolParameter
   */
  template<class T>
  inline SpaceLib::CPhysUnit* InputParameter(
    const std::string&      name
   ,      T&           dest
   ,      ParamList&   p
   ,const std::string& exception_text_prefix
   ,const std::string& trailing_chars
   ,bool check_key_found
   ) throw(std::invalid_argument)
  {
    ParamList::iterator i = find_if(p.begin(),p.end(),name_equal<Param>(name));
    if (i!=p.end()) 
    {
      removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
      removeChar((*i).second, '\"');                    // remove possible quotes
      std::string s_unit(RemoveUnitString((*i).second));
      istringstream((*i).second) >> dest;    // use stream operator on reference
      p.erase(i);                          // delete read parameter from storage
      return SpaceLib::GetUnit(s_unit);
    }
    else if (check_key_found == true) throw invalid_argument(
      exception_text_prefix+string(": missing parameter (")+name+string(")\n"));
    return SpaceLib::U_NO_UNIT;
  };


  /*!
  \brief Assigns string parameter from parameter list to given reference.
  \param name Parameter name.
  \param dest reference to string to hold the parameter value
  \param params the vector of parameters
  \param exception_text_prefix message which will be displayed before
  actual exeception message (tracing)
  \param trailing_chars string of characters which are not allowed after
 a data value (tabs, blanks)
  \exception std::invalid_argument The parameter "name" could not be found in
  given paramter list.
  */
  void InputStringParameter(
    const std::string& name
   ,std::string&       dest
   ,ParamList&         params
   ,const std::string& exception_text_prefix
   ,const std::string& trailing_chars
   ) throw(std::invalid_argument);

  /*!
   *\brief Assigns parameter value from parameter list to given bool variable reference.
   *\param name boolean parameter name.
   *\param dest reference of boolean variable data type to hold the parameter value
   *\param params parameters storage object
   *\param exception_text_prefix message which will be displayed before
   * actual exeception message (for tracing purposes)
   *\param trailing_chars string of characters which should be removed after
   * a data value (eg. tabs, blanks)
   * \param check_key_found throws exception in case of not found parameter key if
   * set to true, otherwise not
   *\exception std::invalid_argument The parameter "name" could not be found in
   *given paramter list (check_key_found must be set to true).
   *\sa InputParameter, InputStringParameter
   */
  void InputBoolParameter(
   const std::string&  name
   ,      bool&        dest
   ,      ParamList&   p
   ,const std::string& exception_text_prefix
   ,const std::string& trailing_chars
   ,bool check_key_found = false
   ) throw(std::invalid_argument);

  /*!
  \brief Assumes filename parameter to be parsed and possibly checks 
  its existence
  \param name Parameter name.
  \param dest reference to string to hold read filename
  \param params parameter storage structure
  \param exception_text_prefix message which will be displayed before
  actual exeception message (tracing)
  \param trailing_chars string of characters which are not allowed after
 a data value (tabs, blanks)
  \param path possible leading path string (will be added to filename)
  \param check_key_found throws exception in case of not found parameter key if
  set to true, otherwise not
  \param check_file_existent if set to true, the read filename will be checked for existence
  \exception std::invalid_argument The parameter "name" could not be found in
  given paramter list.
  */
  void InputFilenameParameter(
    const std::string& name
    ,std::string& dest
    ,ParamList& p
    ,const std::string& exc_prefix
    ,const std::string& trailing_chars
    ,const std::string& path
    ,bool check_key_found
    ,bool check_file_existent
    )throw(std::invalid_argument);

  /*!
   *\brief File output of given array of values.
   *
   * Data of given valarray will be written into file, each array member
   * into a new line. The output will be done with the current output precision
   * of std::cout object. So, if you want to use a certain precision for your
   * file output, all you've got to do is to set std::cout's precision with:\n
   * std::cout.precision(your value);\n
   *\param data std::valarray of data values to dump to file.
   *\param file A reference to std::ofstream object to dump data in.
   */
  template<class T>
  inline
  void FileOutputEndl(T& data, std::ofstream& file)
  {
    if (file.is_open() == true)
    {
      for (uint i=0; i<data.size(); ++i) 
        file << setprecision(cout.precision()) << data[i] << endl;
    }
  }


  /*! \brief Deletes pointer if it's not zero and sets to zero after deleting.
   *  \param p reference to a pointer of template argument type
   *  \return pointer (given by reference parameter p) will be set to zero
   */
  template<class T>
  inline void SafeDelete(T*& p)
  {
    if (p!=0)
    {
      delete p;
      p = 0;
    }
  } // void SafeDelete(T*& p)

  /*! \brief Deletes pointer to array if it's not zero and sets to zero after deleting.
   *  \param p reference of a pointer to template argument type array
   *  \return pointer (given by reference parameter p) will be set to zero
   */
  template<class T>
  inline void SafeDeleteArray(T*& p)
  {
    if (p!=0)
    {
      delete[] p;
      p = 0;
    }
  } // void SafeDelete(T*& p)

} // end of namespace Global

using namespace Global;

#endif // GLOBAL_H

