/*!
\file CAttitude.cpp
\version 1.7
\brief An auxiliary class for attitude calculations (CAttitude).
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 03.09.2005 (AO)

\b Purpose:
   
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

\todo change data representation of quaternion to Vector   
*/

// standard C math library
#include <math.h>
#include <cassert>

//SpaceLib
#include <configSpaceLib.h>
#include <CAttitude.h>
#include <CAttitudeData.h>

namespace SpaceLib
{

  //---------------------------------------------------------------------------
  CAttitude::CAttitude( void ) : data( 0)
  {
    data = new CAttitudeData;
    if ( data == 0 )
      FatalError( "CAttitude::Attitude(void) => Error during memory allocation!" );
  }

  //----------------------------------------------------------------------------
  CAttitude::CAttitude(const CAttitude& to_copy) : data(0)
  {
    data = new CAttitudeData(*to_copy.data);
    if ( data == 0 )
      FatalError( "CAttitude::Attitude(void) => Error during memory allocation!" );
  }


  // Implementation of matrix constructor
  //---------------------------------------------------------------------------
  CAttitude::CAttitude( const MathLib::Matrix & input ) : data( 0 )
  {
    if ( (input.getNumberOfRows() == 3) && (input.getNumberOfColumns() == 3) )
    {
      data = new CAttitudeData( input );
      if ( data != 0 )
      {
        data->updateFromMatrix();
      }
      else
        FatalError( "CAttitude::Attitude(void) => Error during memory allocation!" );
    }
    else
      FatalError( "CAttitude::Attitude(const Matrix&) => Wrong matrix dimension!" );
  }


  // Implementation of quaternion constructor
  //---------------------------------------------------------------------------
  CAttitude::CAttitude( 
      const double q1
    , const double q2
    , const double q3
    , const double q4 ) : data( 0 )
  {
    data = new CAttitudeData( q1, q2, q3, q4 );
    if ( data != 0 )
    {
      data->updateFromQuaternion();
    }
    else
      FatalError( "CAttitude::Attitude(const double q1, const double q2, const double q3, const double q4) => Error during memory allocation!" );
  }


  // Implementation of axis/angle constructor
  //---------------------------------------------------------------------------
  CAttitude::CAttitude( const MathLib::Vector & axis, const double angle ) 
    : data( 0 )
  {
    data = new CAttitudeData( axis, angle );
    if ( data != 0 )
    {
      data->updateFromAxisAngle();
    }
    else
      FatalError( "CAttitude::Attitude(const Vector & axis, const double angle) \
 => Error during memory allocation!" );
  }


  // Implementation of unit vector constructor
  //---------------------------------------------------------------------------
  CAttitude::CAttitude( const MathLib::Vector & e1
                      , const MathLib::Vector & e2
                      , const MathLib::Vector & e3 ) : data( 0 )
  {
    data = new CAttitudeData( e1, e2, e3 );
    if ( data != 0 )
    {
      data->updateFromMatrix();
    }
    else
      FatalError( "CAttitude::Attitude(void) => Error during memory allocation!" );
  }


  // Quaternion class constructor
  //---------------------------------------------------------------------------
  CAttitude::CAttitude( const MathLib::CAttitudeQuaternion& q) : data(0)
  {
    data = new CAttitudeData(q);
    if ( data != 0 )
    {
      data->updateFromQuaternion();
    }
    else
    {
      FatalError( "CAttitude::Attitude(const CAttitudeQuaternion&) => Error during memory allocation!" );
    }
  }


  // Implementation of copy operator
  //---------------------------------------------------------------------------
  CAttitude & CAttitude::operator = ( const CAttitude & to_copy )
  {
    ( * data ) = ( * to_copy.data );
    return *this;
  }


  // Implementation of destructor
  //---------------------------------------------------------------------------
  CAttitude::~CAttitude( void )
  {
    SafeDelete(data);
  }


  // Read Access methods
  // Unit Vector one of coordinate system
  //---------------------------------------------------------------------------
  MathLib::Vector CAttitude::UnitVec1( void )const
  {
    return data->m_matrix.Row( 0 );
  }


  // Unit Vector two of coordinate system
  //---------------------------------------------------------------------------
  MathLib::Vector CAttitude::UnitVec2( void )const
  {
    return data->m_matrix.Row( 1 );
  }


  // Unit Vector three of coordinate system
  //---------------------------------------------------------------------------
  MathLib::Vector CAttitude::UnitVec3( void )const
  {
    return data->m_matrix.Row( 2 );
  }


  // Attitude matrix
  //---------------------------------------------------------------------------
  const MathLib::Matrix& CAttitude::AttitudeMatrix( void )const
  {
    return data->m_matrix;
  }


  // Rotation axis
  //---------------------------------------------------------------------------
  const MathLib::Vector& CAttitude::Axis( void )const
  {
    return data->m_axis;
  }


  // Rotation angle
  //---------------------------------------------------------------------------
  double CAttitude::Angle( void )const
  {
    return data->m_angle;
  }

  // Quaternion representation
  //---------------------------------------------------------------------------
  const MathLib::CAttitudeQuaternion& CAttitude::Quaternion( void )const
  {
    return data->m_q;
  }


  // Write Access methods

  // set whole attitude data structure
  //---------------------------------------------------------------------------
  void CAttitude::Set( const CAttitude & to_copy )
  {
    * this = to_copy;
  }


  // set attitude matrix
  //---------------------------------------------------------------------------
  void CAttitude::Set( const MathLib::Matrix & input )
  {
    if ( (input.getNumberOfRows() == 3) && (input.getNumberOfColumns() == 3) )
    {
      data->m_matrix = input;
      data->updateFromMatrix();
    }
    else
      FatalError( "CAttitude::Set(const Matrix&) => Wrong dimension of input matrix!" );
  }


  // set quaternion elements
  //---------------------------------------------------------------------------
  void CAttitude::Set( const double q1, const double q2, const double q3, const double q4 )
  {
    data->m_q[0] = q1;
    data->m_q[1] = q2;
    data->m_q[2] = q3;
    data->m_q[3] = q4;
    data->updateFromQuaternion();
  }


  // set quaternion elements from vector input
  //---------------------------------------------------------------------------
  void CAttitude::Set( const MathLib::Vector & quat )
  {
    data->m_q.setDataVector(quat);
    data->updateFromQuaternion();
  }

  // set with vector subscope data (quaternion assumed)
  //---------------------------------------------------------------------------
  void CAttitude::set(const MathLib::Vector& d, const unsigned int first)
  {
    data->m_q[0] = d[first];
    data->m_q[1] = d[first+1];
    data->m_q[2] = d[first+2];
    data->m_q[3] = d[first+3];
    data->updateFromQuaternion();
  }


  // set rotation axis and angle
  //---------------------------------------------------------------------------
  void CAttitude::Set( const MathLib::Vector& axis, const double angle )
  {
    data->m_axis = axis;
    data->m_angle = angle;
    data->updateFromAxisAngle();
  }


  // set attitude quaternion
  //---------------------------------------------------------------------------
  void CAttitude::Set( const MathLib::CAttitudeQuaternion& q)
  {
    data->m_q = q;
    data->updateFromQuaternion();
  }


  // set unit vectors
  void CAttitude::Set( 
      const MathLib::Vector & e1
    , const MathLib::Vector & e2
    , const MathLib::Vector & e3 )
  {
    if ( ( e1.getDimension() == 3 ) && ( e2.getDimension() == 3 ) && ( e3.getDimension() == 3 ) )
    {
      data->m_matrix.SetRow( 0, e1 );
      data->m_matrix.SetRow( 1, e2 );
      data->m_matrix.SetRow( 2, e3 );
      data->updateFromMatrix();
    }
    else
      FatalError( "CAttitude::Set(const Vector&, const Vector&, const Vector&) => Wrong dimension of input vector(s)!" );
  };


  // return derivative of attitude matrix
  MathLib::Matrix CAttitude::DerivativeMatrix( const MathLib::Vector & omega ) const {

    MathLib::Matrix Omega(3,3);

    if ( omega.getDimension() == 3 ) {

      Omega( 0, 0 ) = Omega( 1, 1 ) = Omega( 2, 2 ) = 0.0;
      Omega( 1, 0 ) = -( Omega( 0, 1 ) = omega( 2 ) );
      Omega( 0, 2 ) = -( Omega( 0, 2 ) = omega( 1 ) );
      Omega( 2, 1 ) = -( Omega( 1, 2 ) = omega( 0 ) );

    } else FatalError( "Matrix CAttitude::Derivative(const Vector&) => Wrong dimension of angular velocity vector!" );

    return Omega * data->m_matrix;

  }


  // return derivative of quaternion
  //---------------------------------------------------------------------------
  MathLib::CAttitudeQuaternion CAttitude::calcFirstDerivativeQuaternion( const MathLib::Vector & omega ) const {

    return data->m_q.calcFirstDerivative(omega);

  }



  //---------------------------------------------------------------------------
  MathLib::Vector CAttitude::AngularVelocity(const MathLib::Vector& quat_p) const
  {
    MathLib::Vector om(quat_p*2.0); // set right hand side of matrix/vector equation

    MathLib::Matrix Q2(4,4);        // matrix filled with quaternions components

    Q2(0,0)=  data->m_q[3]; Q2(0,1)= -data->m_q[2]; 
    Q2(0,2)=  data->m_q[1]; Q2(0,3)= 0.0;
    Q2(1,0)=  data->m_q[2]; Q2(1,1)=  data->m_q[3];
    Q2(1,2)= -data->m_q[0]; Q2(1,3)= 0.0;
    Q2(2,0)= -data->m_q[1]; Q2(2,1)=  data->m_q[0];
    Q2(2,2)=  data->m_q[3]; Q2(2,3)= 0.0;
    Q2(3,0)= -data->m_q[0]; Q2(3,1)= -data->m_q[1]; 
    Q2(3,2)= -data->m_q[2]; Q2(3,3)= 1.0;

    MathLib::Vector temp(4);        // needed by LU- decomposition

    LU_Decomp(Q2, temp);       // decomposition
    LU_BackSub(Q2, temp, om);  // backsubstitution

    return om.Slice(0,2);  // return angular velocity vector

  }

  //---------------------------------------------------------------------------
  void CAttitude::transformVectorToRefFrame(const MathLib::Vector& s, MathLib::Vector& d)
  {
    assert(s.size()==d.size() && s.size()==3);
    for (int j=0;j<3;++j)
    { 
      d[j]=0.0;
      for (int i=0;i<3;++i) d[j]+=data->m_matrix(i,j)*s[i];
    }
  }

}
