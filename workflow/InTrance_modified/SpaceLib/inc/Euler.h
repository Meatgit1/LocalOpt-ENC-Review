/*!
\file Euler.h
\version 1.7
\brief A class to simplify the handling of rotation sequences/transformation matrices (Euler).
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 07.06.2004 (AO)

\b Purpose:
   Class Euler is intended to make handling of Euler angle rotation
   sequences/transformations easier. It provides functions to
   calculate transformation matrices and also to get the desired Euler angles
   back from a given transformation matrix. These operations are supported
   for every one of the 12 possible rotation sequences (three subsequent
   rotations about three different axis result in 12 possible rotation sequences).

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.
   
\b Copyright:
   (c) 2004 Andreas Ohndorf
*/

#ifndef EULER_H
#define EULER_H

#include <MLMatrix.h>

namespace SpaceLib {

  /*!
  \brief Maximum number of rotation sequences.

   The number of possible rotation sequences for three rotations and three
   possible rotation axis sum up to this.
  */
  const int EULER_N_SEQUENCES = 12;


  /** Number of maximal used Attitude matrix elements to estimate values
      for the Euler angles. They may be different for every rotation sequence. */
  const int EULER_USED_ELEMENTS = 7;

  /*!
    \brief Euler rotation sequence types.

    This is a data type to indicate the Euler rotation sequence, to use for
    estimating the desired three angles. The numbers after the capital E are
    in the order of subsequent rotations. The figures represent the axis of
    rotation. For example, E213 means a rotation about axis 2, followed by
    second rotation about axis 1 and a last one about axis 3. The used angles
    are always angle 1 for the first rotation, angle 2 for the second
    rotation and angle 3 for the third and last rotation. You can see, that
    in general there are two main types of rotation sequences. The ones
    where rotation axis one is the same as rotation axis three and the ones
    where they are different. Those are called type two and type one Euler
    rotation sequences.

    \sa GetEulerTransformation, getEulerAngles
  */
  enum EulerSequences {
    E123 = 0, ///< rotation about x-, y-, z-axis
    E132 = 1, ///< rotation about x-, z-, y-axis
    E231 = 2, ///< rotation about y-, z-, x-axis
    E213 = 3, ///< rotation about y-, x-, z-axis
    E312 = 4, ///< rotation about z-, x-, y-axis
    E321 = 5, ///< rotation about z-, y-, x-axis
    E121 = 6, ///< rotation about x-, y-, x-axis
    E131 = 7, ///< rotation about x-, z-, x-axis
    E212 = 8, ///< rotation about y-, x-, y-axis
    E232 = 9, ///< rotation about y-, z-, y-axis
    E313 = 10,///< rotation about z-, x-, z-axis
    E323 = 11 ///< rotation about z-, y-, z-axis
  };


  // Attitude matrix element.
  // Enumeration to identify a needed element of the Attitude matrix.
  // The content of the element may be different, depending on the
  // type of Euler angle representation.
  //
  // Type 1 representation: first and third rotation axis are not equal
  //  For the Euler sequences E123, C231, E312 there are also different
  //  elements of the Attitude matrix than for the sequences E213, E312, E321.
  //
  // Type 2 representation: first and third rotation axis are equal
  /** Tbd. */
  enum EulerElements {

         E1 = 0, ///< tbd
         // Type 1:sin(a1)/Type 2:cos(a1)

         E2 = 1, ///< tbd
         // Type 1:cos(a2)sin(a1) / Type 2:sin(a2)sin(a1)

         E3 = 2, ///< tbd
         // Type 1:cos(a2)cos(a1) / Type 2:sin(a2)cos(a1)

         E4 = 3, ///< tbd
         // Type 1:cos(a2)sin(a3) / Type 2:sin(a2)sin(a3)

         E5 = 4, ///< tbd
         // Type 1:cos(a2)cos(a3) / Type 2:sin(a2)cos(a3)

         E6 = 5, ///< tbd
         // Type 1, E123,E231,E312: sin(a3)cos(a1)+cos(a3)sin(a2)sin(a1)
         // Type 1, E321,E213,E132: sin(a3)sin(a1)+cos(a3)sin(a2)cos(a1)
         // Type 2,               : cos(a3)cos(a1)-sin(a3)cos(a2)sin(a1)

         E7 = 6 ///< tbd
         // Type 1, E123,E231,E312: cos(a3)sin(a1)+sin(a3)sin(a2)cos(a1)
         // Type 1, E321,E213,E132: cos(a3)cos(a1)+sin(a3)sin(a2)sin(a1)
         // Type 2,               : cos(a3)sin(a1)+sin(a3)cos(a2)cos(a1)
  };



  //----------------------------------------------------------------------------

  /*!
  \brief Supplementary class to handle Euler rotation sequences.
  */
  class Euler {

  public:

    Euler( void ) {};

    ~Euler() {};


    /*!
    \brief Attitude/Transformation matrix resulting from Euler rotations.

    The returned transformation matrix is calculated from the three angles
    given as arguments and the desired type of Euler rotation sequence. Therefor
    three subsequent rotation matrices will be multiplied, which for itself
    will be calculated with the according angle.
    \param seq Desired type of rotation sequence.
    \param a1 Angle of rotation 1
    \param a2 Angle of rotation 2
    \param a3 Angle of rotation 3
    \return Attitude/Transformation matrix
    \sa getEulerAngles, EulerSequences
    */
    MathLib::Matrix getEulerTransformation(
      const EulerSequences seq,
      const double a1,
      const double a2,
      const double a3
    );

    /** Tbd. */
    void getEulerAnglesType1(
      const MathLib::Matrix & T,
      double & a1, double & a2,
      double & a3,
      enum EulerSequences seq
    );

    /** Tbd. */
    void getEulerAnglesType2(
      const MathLib::Matrix & T,
      double & a1, double & a2,
      double & a3,
      enum EulerSequences seq
    );


    /*!
    \brief Tries to obtain Euler angles out of a given Attitude/transformation
    matrix.

    This function computes the three Euler angles from a given transformation
    matrix. Though these angles are usually different for every Euler rotation
    sequence, the desired sequence itself must be provided as parameter.\n
    If the resulting angles and the same rotation sequence are used as arguments
    for getEulerTransformation, one should obtain the same Attitude matrix.

    \param T Attitude/transformation matrix
    \param a1 Reference of angle 1 (return value)
    \param a2 Reference of angle 2 (return value)
    \param a3 Reference of angle 3 (return value)
    \param seq The desired Euler rotation sequence
    \sa getEulerTransformation, EulerSequences
    */
    void getEulerAngles(
      const MathLib::Matrix & T,
      double & a1,
      double & a2,
      double & a3,
      enum EulerSequences seq
    );

  }; // of class Euler

} // of namespace SpaceLib


#endif
