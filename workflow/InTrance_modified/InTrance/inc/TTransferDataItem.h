/*!
\file TTransferDataItem.h
\version 1.7
\brief Transfer class for data exchange between InTrance and Java GUI
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@dlr.de (@gmx.de)> (AO)
\date 04.01.2006 (AO)

\b Purpose:
   This class is used in conjunction with a double ended queue as an
   implementation to exchange date between the InTrance optimizer (as
   dynamic library) and a graphical user interface written in JAVA to
   display certain information about the current optimization process.

\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.

\b Copyright:
   (c) 2006 Bernd Dachwald, Andreas Ohndorf
*/

#ifndef TTRANSFERDATAITEM_H
#define TTRANSFERDATAITEM_H

// STL
#include <valarray>
#include <string>

// SpaceLib
#include <SpaceLib.h>

// JAVA Native Interface (JNI)
#include <jni.h>

// type definitions
typedef std::vector<MathLib::Vector> CTrajectory;
typedef CTrajectory*                 Trajectory;
typedef std::vector<double>          CDataArray;
typedef CDataArray*                  DataArray;


typedef
class TTransferDataItem
{

 public:

  Trajectory m_trajectory;          ///< plot of s/c position
  Trajectory m_approx_traj;         ///< plot of s/c pos due to approx. steering angles
  Trajectory m_target_trajectory;   ///< same for s/c's target

  DataArray  m_angle_alpha_history; ///< history of steering angle alpha

  /// Values from polynomial approximation of steering clock angle.
  DataArray  m_angle_alpha_history_approx;

  DataArray  m_angle_beta_history; ///< history of steering angle beta

  /// Values from polynomial approximation of steering beta angle.
  DataArray  m_angle_beta_history_approx;

  DataArray  m_thrust_history;         ///< thrust history values
  DataArray  m_thrust_history_approx;  ///< approximation of thrust values

  std::string m_reproduction_message; // output of actual best reproduction

  float m_repro_per_sec;

  std::vector<std::string> m_fitness_values;

  /*!
   * \brief deletes internal data structures
   */
  void destroy(void)
  {
    SafeDelete(m_trajectory);
    SafeDelete(m_approx_traj);
    SafeDelete(m_target_trajectory);
    SafeDelete(m_angle_alpha_history);
    SafeDelete(m_angle_alpha_history_approx);
    SafeDelete(m_angle_beta_history);
    SafeDelete(m_angle_beta_history_approx);
    SafeDelete(m_thrust_history);
    SafeDelete(m_thrust_history_approx);
  }

  /// destructor
  ~TTransferDataItem(void) {}

  /// standard constructor
  TTransferDataItem(void) :
      m_trajectory(0)
    , m_approx_traj(0)
    , m_target_trajectory(0)
    , m_angle_alpha_history(0)
    , m_angle_alpha_history_approx(0)
    , m_angle_beta_history(0)
    , m_angle_beta_history_approx(0)
    , m_thrust_history(0)
    , m_thrust_history_approx(0)
    , m_repro_per_sec(0.0)
  {};

  //---------------------------------------------------------------------------
  /// initialize constructor
  TTransferDataItem(
	 Trajectory trajectory,
   Trajectory approx_traj,
	 Trajectory target_trajectory,
   DataArray  angle_alpha_history,
   DataArray  angle_alpha_history_approx,
   DataArray  angle_beta_history,
   DataArray  angle_beta_history_approx,
   DataArray  thrust,
   DataArray  thrust_approx,
	 const std::string& repro_message,
   float repro_per_sec) :
    m_trajectory(trajectory),
    m_approx_traj(approx_traj),
    m_target_trajectory(target_trajectory),
    m_angle_alpha_history(angle_alpha_history),
    m_angle_alpha_history_approx(angle_alpha_history_approx),
    m_angle_beta_history(angle_beta_history),
    m_angle_beta_history_approx(angle_beta_history_approx),
    m_thrust_history(thrust),
    m_thrust_history_approx(thrust_approx),
    m_reproduction_message(repro_message),
    m_repro_per_sec(repro_per_sec)
  {};

  //---------------------------------------------------------------------------
  /// copy constructor
  TTransferDataItem(const TTransferDataItem& item)
  {
    *this = item;
  }

  //---------------------------------------------------------------------------
  /// assign operator
  TTransferDataItem& operator=(const TTransferDataItem& item)
  {
    m_trajectory                  = item.m_trajectory;
    m_approx_traj                 = item.m_approx_traj;
    m_target_trajectory           = item.m_target_trajectory;
    m_angle_alpha_history         = item.m_angle_alpha_history;
    m_angle_alpha_history_approx  = item.m_angle_alpha_history_approx;
    m_angle_beta_history          = item.m_angle_beta_history;
    m_angle_beta_history_approx   = item.m_angle_beta_history_approx;
    m_thrust_history              = item.m_thrust_history;
    m_thrust_history_approx       = item.m_thrust_history_approx;
    m_reproduction_message        = item.m_reproduction_message;
    m_repro_per_sec               = item.m_repro_per_sec;
    m_fitness_values              = item.m_fitness_values;
    return *this;
  }

  //---------------------------------------------------------------------------
  /// lower than operator (just for STL conformity)
  inline bool operator<(const TTransferDataItem& right)
  {
    return true;
  }

  //---------------------------------------------------------------------------
  /// equal operator (just for STL conformity)
  inline bool operator==(const TTransferDataItem& right)
  {
    return true;
  }

} TransferDataItem;

#endif // of #define TTransferDataItem
