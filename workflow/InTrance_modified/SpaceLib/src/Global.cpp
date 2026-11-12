/*!
\file Global.cpp
\version 1.7
\brief Some global functions that are widely used throughout the SpaceLib.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
\date 29.10.2006 (AO)

\b Purpose:
   This file implements some global functions that are used
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


//SpaceLib
#include <configSpaceLib.h>

// STL
#include <algorithm>
#include <functional>
#include <sstream>
#include <cassert>

// SpaceLib
#include <Global.h>


using namespace std;


namespace Global
{

  const std::string NUMBER_ELEMENTS("1234567890.+-eE");

  SLStandardOutput stdOutput;    // standard output object

  SLOutput* myOut = &stdOutput;  // global output object

  //---------------------------------------------------------------------------
  SLOutput& SLOutput::operator<<(const char* message)
  {
    if (message != 0)
    {
      std::ofstream os(log_filename.c_str(), ios_base::app);
      if ( os.is_open() == true)
      {
        os << message;
      }
    }
    return *this;
  }

  //---------------------------------------------------------------------------
  void setOutput(SLOutput* output)
  {
    myOut = output;
  }

  //---------------------------------------------------------------------------
  void restoreStandardOutput(void)
  {
    myOut = &stdOutput;
  }

  //---------------------------------------------------------------------------
  void setOutputLogFilename(const std::string filename)
  {
    if (myOut != 0) myOut->setLogFilename(filename);
  }

  //---------------------------------------------------------------------------
  const char* getOutputLogFilename(void)
  {
    if (myOut != 0) return myOut->getLogFilename();
    return 0;
  }

  //---------------------------------------------------------------------------
  void Output(const std::string& message)
  {
    if (myOut != 0) (*myOut) << message;
  }

  //---------------------------------------------------------------------------
  void Output(const char* message)
  {
    if (myOut != 0 && message != 0) (*myOut) << message;
  }

  //---------------------------------------------------------------------------
  string ReadString(ifstream& in)
  {
    string str;
    getline(in, str);
    return str;
  }

  //---------------------------------------------------------------------------
  double ReadNumber(ifstream& in)
  {
    double number;
    in >> number;
    in.ignore(81, '\n');
    return number;
  }


  //---------------------------------------------------------------------------
  string& upper(string& s)
  {
    for (std::string::size_type i = 0;i<s.length();++i)
      s[i] = static_cast<char>(toupper(s[i]));
    return s;
  } // std::string& void upper(string& s)


  //---------------------------------------------------------------------------
  void removeChar(string& s, const char character)
  {
	string::size_type pos = string::npos;
	while ((s.empty() != true) && 
		   ((pos=s.find(character)) != string::npos)) s.erase(pos,1);
  }

  //---------------------------------------------------------------------------
  void removeChars(string& s, const string& characters)
  {
	string::size_type pos = string::npos;
	while ((s.empty() != true) &&
		((pos=s.find_first_of(characters)) != string::npos))
      s.erase(pos,1);
  }

  //---------------------------------------------------------------------------
  void removeComments(string& s, const string& sep)
  {
	  string::size_type pos = s.find(sep);
    if (pos!=std::string::npos)
    {
      std::string t(s.begin(), s.begin()+pos);
      s.assign(t);
    }
//	  if (pos != string::npos) s.erase(/*s.begin()+*/pos);//,s.end());
  }

  //---------------------------------------------------------------------------
  void removeTrailingChars(string& s, const string& trail_chars)
  {
    string::size_type i = 0;
    while ((i=s.find_last_of(trail_chars)) != string::npos && i == (s.length()-1) )
      s.erase(s.length()-1,1);
  }

  //----------------------------------------------------------------------------
  void removeLeftWhitespaces(
            string&            s
    ,       string::size_type& pos
    , const string&            whitespaces
)
  {
    if (pos>0 && pos!=string::npos)
    {
      string::size_type i(pos);
      while (whitespaces.find_first_of(s[i-1])!=string::npos && i>0) --i;
      if (i!=0 && i<pos) s.erase(i, pos-i);
      pos-=(pos-i);
    }
  }

  //----------------------------------------------------------------------------
  void removeLeadingWhitespaces(std::string& s,  const std::string& whitespaces)
  {
    if (s.size()>0)
    {
      string::size_type i(s.find_first_not_of(whitespaces));
      if (i>0 && i!=string::npos) s.erase(0, i);
    }
  }


  //----------------------------------------------------------------------------
  void removeRightWhitespaces(
            string&           s
    ,       string::size_type pos
    , const string&           whitespaces
)
  {
    if (pos>0 && pos!=string::npos)
    {
      string::size_type i(pos);
      while (whitespaces.find_first_of(s[i+1])!=string::npos && i<s.size()) ++i;
      if (i!=s.size() && i>pos) s.erase(pos+1, i-pos);

    }
  }


  //--------------------------------------------------------------------------
  void StripPathFromFilename(std::string& name)
  {
    std::string::size_type last=name.find_last_of("/\\");
    std::string::size_type first=name.find_first_of("/\\");
    if (last == std::string::npos) name.clear();
    else name.erase(last+1, name.length()-last-1);
  }


  //---------------------------------------------------------------------------
  ParamList GetParameters(ifstream& ifile, const string& comment_seperator)
    throw (syntax_error_in_config_file)
  {
    string::size_type pos;    // position in string; used for find operations
    string ws;                // work string
    ParamList tempList;       // final list of read parameter names and values
    int line = 0;             // line counter for error message
    ostringstream linestring; // stringstream to convert line counter (integer)
    const string whitespaces(" \t");

    while (!ifile.eof())
    {
      ws.resize(0);
      getline(ifile, ws);      // read line from file
      ++line;                  // line counter
      linestring.str("");      // initialize
      linestring << line;      // read from file

      // remove comments
	    removeComments(ws, comment_seperator);

	    // remove possible trailing spaces
	    removeTrailingChars(ws, whitespaces);

      // this block won't be entered, if the workstring is empty
      if (!ws.empty())
      {
        removeLeadingWhitespaces(ws, whitespaces);
        if ((pos=ws.find('=')) != string::npos)
        {
          removeLeftWhitespaces(ws, pos, whitespaces);
          removeRightWhitespaces(ws, pos, whitespaces);
          tempList.insert(tempList.end(),
                          Param(ws.substr(0,pos),
                                ws.substr(pos+1, ws.size()-pos-1)));
        }
        else
          throw syntax_error_in_config_file(string("syntax error in line ")+linestring.str());
      }
    }
    return tempList;
  }

  
  //---------------------------------------------------------------------------
  std::string RemoveUnitString(std::string& s)
  {    
    //std::string::iterator i =        
    //  find_first_of(s.rbegin(), s.rend(), NUMBER_ELEMENTS.begin()
    //  , NUMBER_ELEMENTS.end());//, std::not_equal_to<std::string::value_type>());
    using namespace std;
    string::iterator i = s.begin();

    /*string::iterator ref_begin = NUMBER_ELEMENTS.begin();
    string::iterator ref_end(NUMBER_ELEMENTS.end());*/

    while (i != s.end())
    {
      if (find(NUMBER_ELEMENTS.begin(), NUMBER_ELEMENTS.end(), *i)
        == NUMBER_ELEMENTS.end()) break;
      ++i;
    }

    if (i != s.end())
    {
      std::string s_unit(i, s.end());
      s.erase(i, s.end());
      upper(s_unit);
      return s_unit;
    }
    else return std::string("");
  } // std::string RemoveUnitString(std::string& s)

  //---------------------------------------------------------------------------
  void InputStringParameter(
    const std::string& name
   ,std::string& dest
   ,ParamList& params
   ,const std::string& exception_text_prefix
   ,const std::string& trailing_chars
   )throw(std::invalid_argument)
  {
    ParamList::iterator i = find_if(params.begin()
      , params.end(), name_equal<Param>(name));
    if (i!=params.end()) 
    {
      removeTrailingChars((*i).second, trailing_chars); // remove not needed chars
      removeChar((*i).second, '\"');                    // remove possible quotes
      dest=(*i).second;
    }
    else throw invalid_argument(
      exception_text_prefix+
      string(": missing parameter (")+
      name+
      string(")\n"));
      params.erase(i);
    }

  //---------------------------------------------------------------------------
  void InputBoolParameter(
   const std::string&  name
   ,      bool&        dest
   ,      ParamList&   p
   ,const std::string& exception_text_prefix
   ,const std::string& trailing_chars
   ,bool check_key_found
   ) throw(std::invalid_argument)
  {
    std::string s;
    try
    {
      InputStringParameter(name, s, p, exception_text_prefix, trailing_chars);
      upper(s);
      if ((s == "TRUE") || (s=="WAHR") || (s=="YES") || (s=="JA")) dest=true;
      else dest = false;
    }
    catch (std::invalid_argument& e)
    {
      if (check_key_found == true) throw e;
    }
  }


  //---------------------------------------------------------------------------
  void InputFilenameParameter(
    const std::string& name
   ,std::string& dest
   ,ParamList& p
   ,const std::string& exc_prefix
   ,const std::string& trailing_chars
   ,const std::string& path
   ,bool check_key_found
   ,bool check_file_existent
   )throw(std::invalid_argument)
  {
    try 
    {
      InputStringParameter(name, dest, p, exc_prefix, trailing_chars);
    }
    catch (std::invalid_argument& e)
    {
      if (check_key_found==true) throw e;
    }
    dest=path+dest;
    if (check_file_existent==true)
    {
      if (InputFileExistent(dest)!=true)
        throw std::invalid_argument(exc_prefix+"\""+dest+"\" file not found!\n");
    }
  }
} // end of namespace Global

