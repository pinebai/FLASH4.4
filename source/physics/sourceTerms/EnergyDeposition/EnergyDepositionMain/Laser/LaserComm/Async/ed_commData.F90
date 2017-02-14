!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Async/ed_commData
!!
!! NAME
!!  ed_commData
!!
!! SYNOPSIS
!!
!!  use ed_commData
!!
!! DESCRIPTION 
!!  
!!  This is the data module for the energy deposition communication kernels
!!   
!!***

#include "constants.h"
#include "Flash.h"

module ed_commData
  !$ use omp_lib
  implicit none
  !$ integer (kind=omp_lock_kind) :: ed_commLock

  integer, save :: ed_commGlobalRays
  integer, save :: ed_commRayPosIndex(MDIM)
  integer, save :: ed_commChannelSize

  integer, save :: ed_commRaysBetweenMsgTest
  integer, save :: ed_commRaysUntilMsgTest

  integer, save :: ed_commPipelineState
  integer, parameter :: PIPELINE_ACTIVE = -100
  integer, parameter :: PIPELINE_CLOSING = -200
  integer, parameter :: PIPELINE_INACTIVE = -300

  integer, save :: ed_commCounterState
  integer, parameter :: COUNTER_ACTIVE = -400
  integer, parameter :: COUNTER_INACTIVE = -500

  integer, save :: ed_commLogUnit
  logical, save :: ed_commLog
  logical, save :: ed_commDebug
end module ed_commData
