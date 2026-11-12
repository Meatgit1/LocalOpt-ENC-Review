/*!
\file CLagrange.cpp
\version 1.7
\brief A class derived from class CBody for Lagrange equilibrium points (CLagrange).
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 23.08.2006 (AO)

\b Purpose:
   Class CLagrange makes it possible to treat a Lagrangian equlibrium point
   like any other body in solar system (planets, asteroids, comets). It uses
   some functionality of CSolarSystem to calculate position and velocity
   data of the L-point and sets its own data members accordingly. So you
   may get those information from like you get them from any other object
   which is stored within a CSolarSystem object.
   The name of a Lagrange point is rather simple and straightforward. It
   starts with the capitel letter "L" followed by the number of the L - point,
   then a underscore, the name of body one (usually sun), followed by another
   underscore and finally the name of body two (jupiter for example).
   As an example, the name of L4 in Sun-Jupiter system would be:
   "L4_sun_jupiter".
   The name may change, if the member function "getName" of CBody class changes.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.
   
\b Copyright:
   (c) 2005 Andreas Ohndorf
*/

//SpaceLib
#include <configSpaceLib.h>

#include <sstream>

#include <CLagrange.h>
#include <CStateData.h>

using namespace std;

namespace SpaceLib {

  // constructor
  // --------------------------------------------------------------------------
  CLagrange::CLagrange( CBody* b1,  CBody* b2, const int which) :
   CBody()
   ,L_Body_1(b1->getId())
   ,L_Body_2(b2->getId())
   ,m_LPoint(which) {

     // initialization of CState's data

     if ((which>=1 && which<=5) && (b1->getId()==b2->getCenterBodyId())) {

       stateData->useDE405   = false;
       stateData->indexDE405 = DE405_NO_BODY;
       CState::setType(SPECIAL);       // store
       CState::setRefSolSys(b1->getRefSolSys()==b2->getRefSolSys() ? b1->getRefSolSys() : 0);    // store the pointer to a solar system object
       CState::setGM(G_SI,SI);          // set only with gravitational constant to prevent a zero value
//       if (b1->Index()==b2->CentralBody()) {
         stateData->m_central_state = &(b1->State(b2->MJD()));
//       }
       CState::setCenterBodyId(b2->getCenterBodyId());
       stateData->unitType   = b2->Unit();
       stateData->m_mjd      = b2->MJD()-1;

       ostringstream t;
       t << which;

       CBody::setName("L"+t.str()+"_"+string(b1->getName())+"_"+string(b2->getName()));
       SetState(b2->MJD());

     }
  }



  // Sets the position and velocity data
  // --------------------------------------------------------------------------
  void CLagrange::SetState(const double mjd) {
    Vector s(stateData->refSys->getBody(L_Body_2)->LagrangePoint(m_LPoint, mjd));
    CState::SetState(
      s(0)
     ,s(1)
     ,s(2)
     ,s(3)
     ,s(4)
     ,s(5)
     ,CARTESIAN
     ,mjd
     ,stateData->m_central_body
     ,stateData->unitType);
  }


  // destructor
  // --------------------------------------------------------------------------
  CLagrange::~CLagrange(void)
  {
#ifdef DEBUG
	  cout << "CLagrange::~CLagrange()=>entering!"<< std::endl;std::cout.flush();
#endif

#ifdef DEBUG
	  cout << "CLagrange::~CLagrange()=>leaving!"<< std::endl;std::cout.flush();
#endif
  }

}
