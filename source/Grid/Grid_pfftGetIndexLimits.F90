!!****f* source/Grid/Grid_pfftGetIndexLimits
!!
!! NAME 
!!
!!   Grid_pfftGetIndexLimits
!!
!! SYNOPSIS
!!
!!   Grid_pfftGetIndexLimits(integer(OUT) :: configLimits(MDIM),
!!                           integer(OUT) :: phaseLimits(MDIM))
!!
!! DESCRIPTION 
!!
!!  Return the starting and ending points of the local storage
!!  in configuration and phase spaces, relative to the global
!!  domain. The configuration phase usually refers to the data
!!  at the input, and phase space to the data after forward transform
!!  
!! ARGUMENTS
!!
!!  configLimits  - endpoints in configuration space
!!  phaseLimits  - endpoints in phase space
!!
!!***
subroutine Grid_pfftGetIndexLimits(configLimits,phaseLimits)
#include "constants.h"
  implicit none 
  integer,dimension(LOW:HIGH,MDIM),intent(OUT) :: configLimits, phaseLimits
  configLimits=1
  phaseLimits=1
  return
end subroutine Grid_pfftGetIndexLimits
