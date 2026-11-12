/*!
\file CAttitude.h
\version 1.7
\brief Declaration of an auxiliary class for astrodynamical 
       attitude calculations (CAttitude).
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 03.09.2005 (AO)

\b Purpose:\n

  Class CAttitude is designed to support attitude control tasks by providing
  a data representation of attitude data and functions to convert between
  them. Basically there are four kinds of attitude representation supported.\n
  These are:\n
  - Unit Vector (unit vector of new frame expressed in coordinates of 
  original frame)\n
  - Attitude matrix (transformation to get from a original frame to new one).\n
  - Euler Axis/Angle (Axis you have to rotate original frame and Angle to)\n
  - Quaternion (derived from Euler Axis/Angle representation.)\n
  do so.)

   
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

#ifndef CATTITUDE_H
#define CATTITUDE_H

// MathLib
#include <MLMatrix.h>
#include <MLQuaternion.h>


namespace SpaceLib
{

  struct CAttitudeData;


  /// orientation (attitude) object class
  //-----------------------------------------------------------------------------
  class CAttitude
  {
  private:

    CAttitudeData * data; ///< internal data

     public:

    /*!
    \brief standard constructor.
    */
    CAttitude( void );

    /*!
    \brief copy constructor
    */
    CAttitude(const CAttitude& to_copy);
    
    /*! 
    \brief Attitude matrix constructor. 
    */
    CAttitude( const MathLib::Matrix & input );
    
    /*!
    \brief Quaternion constructor.
    */
    CAttitude( const double q1, const double q2, const double q3
             , const double q4 );
    
    /*!
    \brief Axis/Angle constructor.
    */
    CAttitude( const MathLib::Vector & axis, const double angle );
    
    /*! 
    \brief Unit vectors constructor.
    */
    CAttitude( 
        const MathLib::Vector & e1
      , const MathLib::Vector & e2
      , const MathLib::Vector & e3 );

    /*!
    \brief Quaternion class constructor.
    */
    CAttitude( const MathLib::CAttitudeQuaternion& q);
    
    /*!
    \brief Destructor.
    */
    ~CAttitude( void );
    
    /*!
    \brief Copy operator.
    */
    CAttitude & operator = ( const CAttitude & to_copy );

    bool operator==(const CAttitude& rhs) const
    {
      return (data==rhs.data);
    }

    // Read Access methods
    /*!
    \brief Frame unit vector one.
    */
    MathLib::Vector UnitVec1( void )const;
    
    /*!
    \brief Frame unit vector two.
    */
    MathLib::Vector UnitVec2( void )const;
    
    /*!
    \brief Frame unit vector three.
    */
    MathLib::Vector UnitVec3( void )const;
    
    /*!
    \brief Attitude representating Attitude Matrix.
    */
    const MathLib::Matrix& AttitudeMatrix( void )const;
    
    /*!
    \brief Axis vector of axis/angle representation.

    The axis of rotation is an Eigenvector of the attitude matrix, 
    that means this vector is left unmodified by the transformation
    of the attitude matrix (A*e=e).
    \sa Angle
    */
    const MathLib::Vector& Axis( void )const;
    
    /*!
    \brief Angle to rotate inertial frame about a certain axis to get
     descibed frame.

     The mentioned angle is the the Eigenvalue of the attitude matrix
     which also fully describes the attitude of our system.
    \sa Axis
    */
    double Angle( void )const;
    

    /*! 
    \brief Attitude expressed as quaternion (CAttitudeQuaternion class). 
    */
    const MathLib::CAttitudeQuaternion& Quaternion(void) const;


    // Write Access methods
    /*!
    \brief Sets attitude with another CAttitude object.
    */
    void Set( const CAttitude & to_copy );
    
    /*!
    \brief Sets attitude with given attitude matrix.
    */
    void Set( const MathLib::Matrix & input );
    
    /*!
    \brief Sets attitude with given symmetric Euler paramters (quaternion).
    */
    void Set( const double q1, const double q2, const double q3, const double q4 );
    
    /*!
    \brief Sets attitude with given quaternion.
    */
    void Set( const MathLib::Vector & quat );

    /// set attitude with quaternion data (4 values) from given vector subscope
    void set(const MathLib::Vector& data, const unsigned int first_index);
    
    /*!
    \brief Sets attitude with given axis/angle data.
    */
    void Set( const MathLib::Vector& axis, const double angle );

    /*!
    \brief Sets attitude with given rotation quaternion.
    */
    void Set( const MathLib::CAttitudeQuaternion& q);
    
    /*!
    \brief Sets attitude with given unit vectors of new frame.

     The change from unit vector to attitude matrix representation
     is rather simple, because it's virtually the same. The perpendicular
     unit vectors of new frame are arranged as subsequent rows and
     that's already the desired attitude matrix.

    \param e1 Unit vector one of new frame.
    \param e2 Unit vector two of new frame.
    \param e3 Unit vector three of new frame.
    */
    void Set( 
       const MathLib::Vector & e1
     , const MathLib::Vector & e2
     , const MathLib::Vector & e3 );

    /*!
    \brief Returns derivative of attitude matrix with repect to time.

    \param omega Angular velocity vector expressed in transformed frame 
    (which is described by underlying quaternion).
    \return First derivative with respect to time of attitude Matrix.
    */
    MathLib::Matrix DerivativeMatrix( const MathLib::Vector & omega ) const;
    
    /*!
    \brief Returns derivative of quaternion with respect to time.

    \param omega Angular velocity vector expressed in transformed frame 
    (which is described by underlying quaternion).

    \return Vector filled with first derivatives of quaternion components.
    */
    MathLib::CAttitudeQuaternion calcFirstDerivativeQuaternion( 
      const MathLib::Vector & omega ) const;


    //-------------------------------------------------------------------------
    /*!
    \brief Calculates angular velocity vector (body frame coord.) from given 
    first derivative of quaternion


    \param quat_p dq/dt = relativ motion of this frame
    \return Angular velocity vector (in attitude frame coordinates)
    */
    MathLib::Vector AngularVelocity(const MathLib::Vector& quat_p) const;

    /*!
    \brief Transforms given vector s (which coordinates are assumed to
    be expressed in coordinates of frame described by this CAttitude object) into
    coordinates of this attitude's reference frame.

    \param s source vector which is intended to get transformed
    \param d reference to transformed vector
    */
    void transformVectorToRefFrame(const MathLib::Vector& s, MathLib::Vector& d);

  }; // class CAttitude

} // namespace SpaceLib

#endif
