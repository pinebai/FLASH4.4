!!****if* source/monitors/Timers/TimersMain/Tau/tmr_findTimerIndex
!!
!! NAME
!!   tmr_findTimerIndex - find a timer with a given name
!!
!! SYNOPSIS
!!   tmr_findTimerIndex(character(:), intent(IN) :: name, 
!!                      logical, intent(IN) :: createIfNone,
!!                      integer, intent(OUT) :: index)
!!
!! DESCRIPTION
!!   Find the integer key for a given name name, optionally 
!!   create it if it doesn't exist
!!
!! ARGUMENTS
!!   name --         a string containing the name of the timer to find
!!   createIfNone --  a logical determining if the routine is to create
!!                    a timer if it can't find it.
!!   index --        index of the timer
!!
!! PARAMETERS
!!
!!***

#include "constants.h"

subroutine tmr_findTimerIndex (name, createIfNone, index)

  use Timers_data, ONLY : tmr_prefixLen, &
       tmr_tauList, tmr_freeSlot, tmr_customPrefix, tmr_MAX_CUSTOM_TIMERS
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none  

  character(len=*), intent(in) :: name
  logical, intent(in)          :: createIfNone
  integer, intent(out) :: index

  character(len=MAX_STRING_LENGTH) :: storedTimerString
  integer :: iTmr, argStringLen, storedStringSize
  
  index = NONEXISTENT

  argStringLen = len(name)
  if (argStringLen + tmr_prefixLen > MAX_STRING_LENGTH) then
     call Driver_abortFlash("Following string name is too long:"//name)
  end if


  !Check if timer already exists.  If so, we return the timer's index.
  iTmr = 1
  do while (iTmr < tmr_freeSlot)

     storedStringSize = len(tmr_tauList(iTmr) % tauString)
     storedTimerString = tmr_tauList(iTmr) % &
          tauString(tmr_prefixLen+1:storedStringSize)
     
     if (trim(storedTimerString) == name) then
        index = iTmr
        return
     end if

     iTmr = iTmr + 1
  end do


  !If the timer does not exist, we optionally create the timer and initialise 
  !the tau specific data to 0 (as per examples in tau user-guide).
  if (createIfNone .eqv. .true.) then

     if ( tmr_freeSlot > tmr_MAX_CUSTOM_TIMERS) then
        call Driver_abortFlash("[tmr_findTimerIndex]: No space to add another timer!")
     end if
     
     iTmr = tmr_freeSlot
     tmr_tauList(iTmr) % tauSavedData(:) = 0
     tmr_tauList(iTmr) % tauString = tmr_customPrefix // name
     
     call TAU_PROFILE_TIMER ( &
          tmr_tauList(iTmr) % tauSavedData, &
          tmr_tauList(iTmr) % tauString )
     
     tmr_freeSlot = tmr_freeSlot + 1
     index = iTmr
  end if
     
end subroutine tmr_findTimerIndex
