!!****if* source/monitors/Timers/TimersMain/Tau/Timers_data
!!
!! NAME
!!  Timers_data
!!
!! SYNOPSIS
!!
!!  use Timers_data
!!
!! DESCRIPTION
!!
!!  Holds the data needed by the Timers Unit
!!
!!***


module Timers_data

#include "constants.h"

  integer, parameter :: tmr_MAX_CUSTOM_TIMERS = 500
  character(len=*), parameter :: tmr_customPrefix = "*** custom:"

  type tauTimerObj
     integer, dimension(2) :: tauSavedData
     character (len=MAX_STRING_LENGTH) :: tauString
     logical :: timerStarted
  end type tauTimerObj

  type (tauTimerObj), save, &
       dimension(tmr_MAX_CUSTOM_TIMERS) :: tmr_tauList

  integer, save :: tmr_freeSlot, tmr_globalMe, tmr_globalNumProcs, tmr_prefixLen

end module Timers_data
