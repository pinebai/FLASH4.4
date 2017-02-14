!!****f* source/RuntimeParameters/RuntimeParameters_getAll
!!
!!
!! NAME
!!  RuntimeParameters_getAll
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_getAll(          integer(inout) :: num, 
!!            character(len=MAX_STRING_LENGTH)(inout) :: names(num), 
!!                            real/int/str/log(inout) :: values(num), 
!!                                 logical(inout)     :: changed(num))
!!
!! DESCRIPTION
!!
!!  This function gets all the parameters of a given type from the
!!  RuntimeParameters database, and indicates whether they have changed
!!  during the current restart. 
!!  
!!  The routines below are implemented; they are not overloaded,
!!    so call them directly. 
!!  RuntimeParameter_getAllReal
!!  RuntimeParameter_getAllInt
!!  RuntimeParameter_getAllStr
!!  RuntimeParameter_getAllLog
!!
!! ARGUMENTS
!!
!!
!! num:        number of parameter to get
!! names:      names of parameters
!! values:     values of parameters
!! changed:    logicals indicating if parameter changed or not 
!!             from the initial run
!!
!!
!!***



   
subroutine RuntimeParameters_getAllReal (num, names, values, changed)

implicit none
#include "constants.h"

  integer, intent(inout)                                :: num
  character(len=MAX_STRING_LENGTH), intent(inout)          :: names(num)
  real, intent(inout)                                   :: values(num)
  logical, intent(inout)                                :: changed(num)
  
end subroutine RuntimeParameters_getAllReal


   
subroutine RuntimeParameters_getAllInt (num, names, values, changed)

implicit none
#include "constants.h"

  integer, intent(inout)                                :: num
  character(len=MAX_STRING_LENGTH), intent(inout)       :: names(num)
  integer, intent(inout)                                :: values(num)
  logical, intent(inout)                                :: changed(num)
  
end subroutine RuntimeParameters_getAllInt


   
subroutine RuntimeParameters_getAllStr (num, names, values, changed)

implicit none
#include "constants.h"

  integer, intent(inout)                                :: num
  character(len=MAX_STRING_LENGTH), intent(inout)       :: names(num), values(num)
  logical, intent(inout)                                :: changed(num)
  
end subroutine RuntimeParameters_getAllStr


   
subroutine RuntimeParameters_getAllLog (num, names, values, changed)

implicit none
#include "constants.h"

  integer, intent(inout)                                :: num
  character(len=MAX_STRING_LENGTH), intent(inout)          :: names(num)
  logical, intent(inout)                                   :: values(num)
  logical, intent(inout)                                :: changed(num)
  
end subroutine RuntimeParameters_getAllLog



