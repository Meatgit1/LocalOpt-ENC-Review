// file     : CBodyData.h
// author   : Andreas Ohndorf <andreas.ohndorf@gmx.de> (AO)
// compiler : GNU C++ 3.3.1, MS VC++ 6.0, MS VC++ .NET, Borland C++ 
// date     : 06.11.2006 (AO)
// purpose  : Definition of the data structure where the additional data of a
//            CBody class is stored in.

#ifndef CPLANETDATA_H
#define CPLANETDATA_H

#include <iostream>
#include <SpaceConst.h>

namespace SpaceLib
{
  class DE405;

  /// implementation of internal data structure of CBody
  class CBodyData 
  {
  public:

    DE405* m_DE;                      // a possible accessor class to obtain the
                                             // body's data from a JPL data file
    BodyID m_index;                       // planet's index (see "SpaceConst.h")
    std::string m_name;                                         // planet's name
    double m_epoch;                // epoch of osculatating elements (base time)
    Vector m_elements;// planetary mean orbital elements for epoch of osculation
    Vector m_delements;           // planetary orbital elements centennial rates

    double vrml_emissiveColor[3];      ///< emissive color of body for VRML vis.
    double vrml_diffuseColor[3];        ///< diffuse color of body for VRML vis.
    double vrml_transparency;           ///< value of transparency for VRML vis.

    // interface

    CBodyData(void);                                   ///< standard constructor
    CBodyData(const CBodyData& to_copy);                   ///< copy constructor
    virtual ~CBodyData(void) throw() {};                         ///< destructor
    const CBodyData &operator=(const CBodyData& _cp) throw();///<assign operator
    const bool operator==(const CBodyData& rhs) const;       ///< equal operator

    void setVrmlEmissiveColor(
      const double red, const double green, const double blue) throw();

    const double getVrmlEmissiveColorRed() const throw();
    const double getVrmlEmissiveColorGreen() const throw();
    const double getVrmlEmissiveColorBlue() const throw();

    void setVrmlDiffuseColor(
      const double red, const double green, const double blue) throw();

    const double getVrmlDiffuseColorRed() const throw();
    const double getVrmlDiffuseColorGreen() const throw();
    const double getVrmlDiffuseColorBlue() const throw();

    void setVrmlTransparency(const double transparency) throw();

    const double getVrmlTransparency(void) const throw();

  }; // of struct CBodyData

  //----------------------------------------------------------------------------
  inline CBodyData::CBodyData(void)
    : m_DE(0), m_index(SUN), m_name(), m_epoch(MJD_J2000), m_elements(6),
    m_delements(6) 
  {
    vrml_emissiveColor[0] = vrml_emissiveColor[1] = vrml_emissiveColor[2]=0;
    vrml_diffuseColor[0] = vrml_diffuseColor[1] = vrml_diffuseColor[2]=0;
    vrml_transparency = 0;
  } // standard constructor

  //----------------------------------------------------------------------------
  inline CBodyData::CBodyData(const CBodyData& to_copy)
    : m_DE(to_copy.m_DE),m_index(to_copy.m_index),m_name(to_copy.m_name)
    ,m_epoch(to_copy.m_epoch),m_elements(to_copy.m_elements)
    ,m_delements(to_copy.m_delements)
  {
    for (int i=0;i<3;++i)
    {
      vrml_emissiveColor[i] = to_copy.vrml_emissiveColor[i];
      vrml_diffuseColor[i] = to_copy.vrml_diffuseColor[i];
    }
    vrml_transparency = to_copy.vrml_transparency;
  } // copy constructor

  //----------------------------------------------------------------------------
  inline const CBodyData& CBodyData::operator=(const CBodyData& _cp)
  {
    m_DE = _cp.m_DE;
    m_index = _cp.m_index;
    m_name = _cp.m_name;
    m_epoch = _cp.m_epoch;
    m_elements = _cp.m_elements;
    m_delements = _cp.m_delements;
    for (int i=0;i<3;++i)
    {
      vrml_emissiveColor[i] = _cp.vrml_emissiveColor[i];
      vrml_diffuseColor[i] = _cp.vrml_diffuseColor[i];
    }
    vrml_transparency = _cp.vrml_transparency;
    return (*this);
  } // copy operator

  //----------------------------------------------------------------------------
  inline const bool CBodyData::operator==(const CBodyData& rhs) const
  {
    for (int i=0;i<3;++i)
    {
      if (vrml_emissiveColor[i]!=rhs.vrml_emissiveColor[i]) return false;
      if (vrml_diffuseColor[i]!=rhs.vrml_diffuseColor[i]) return false;
    }
    return (m_DE == rhs.m_DE && m_index==rhs.m_index &&
      m_name==rhs.m_name && m_epoch==rhs.m_epoch && m_elements==rhs.m_elements
      && m_delements==rhs.m_delements 
      && vrml_transparency==rhs.vrml_transparency);
  } // const bool CBodyData::operator=(const CBodyData& rhs)

  //----------------------------------------------------------------------------
  inline void CBodyData::setVrmlEmissiveColor(
    const double red, const double green, const double blue) throw()
  {
    vrml_emissiveColor[0] = red;
    vrml_emissiveColor[1] = green;
    vrml_emissiveColor[2] = blue;
  } // void CBodyData::setVrmlEmissiveColor(...)

  //----------------------------------------------------------------------------
  inline const double CBodyData::getVrmlEmissiveColorRed() const throw()
  {
    return vrml_emissiveColor[0];
  } // const double CBodyData::getVrmlEmissiveColorRed() const throw()

  //----------------------------------------------------------------------------
  inline const double CBodyData::getVrmlEmissiveColorGreen() const throw()
  {
    return vrml_emissiveColor[1];
  } // const double CBodyData::getVrmlEmissiveColorGreen() const throw()

  //----------------------------------------------------------------------------
  inline const double CBodyData::getVrmlEmissiveColorBlue() const throw()
  {
    return vrml_emissiveColor[2];
  } // const double CBodyData::getVrmlEmissiveColorBlue() const throw()

  //----------------------------------------------------------------------------
  inline void CBodyData::setVrmlDiffuseColor(
    const double red, const double green, const double blue) throw()
  {
    vrml_diffuseColor[0] = red;
    vrml_diffuseColor[1] = green;
    vrml_diffuseColor[2] = blue;
  } // void CBodyData::setVrmlDiffuseColor(...)

  //----------------------------------------------------------------------------
  inline const double CBodyData::getVrmlDiffuseColorRed() const throw()
  {
    return vrml_diffuseColor[0];
  } // const double CBodyData::getVrmlDiffuseColorRed() const throw()

  //----------------------------------------------------------------------------
  inline const double CBodyData::getVrmlDiffuseColorGreen() const throw()
  {
    return vrml_diffuseColor[1];
  } // const double CBodyData::getVrmlDiffuseColorGreen() const throw()

  //----------------------------------------------------------------------------
  inline const double CBodyData::getVrmlDiffuseColorBlue() const throw()
  {
    return vrml_diffuseColor[2];
  } // const double CBodyData::getVrmlDiffuseColorBlue() const throw()

  //----------------------------------------------------------------------------
  inline void CBodyData::setVrmlTransparency(const double transparency) throw()
  {
    vrml_transparency = transparency;
  } // void CBodyData::setVrmlTransparency(const double transparency) throw()

  //----------------------------------------------------------------------------
  inline const double CBodyData::getVrmlTransparency(void) const throw()
  {
    return vrml_transparency;
  } // const double CBodyData::getVrmlTransparency(void) const throw()

} // namespace SpaceLib

#endif // CPLANETDATA_H

// end of file

