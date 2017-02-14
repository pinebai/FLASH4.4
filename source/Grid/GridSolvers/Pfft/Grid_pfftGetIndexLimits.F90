!!****if* source/Grid/GridSolvers/Pfft/Grid_pfftGetIndexLimits
!!
!! NAME 
!!
!!   Grid_pfftGetIndexLimits
!!
!! SYNOPSIS
!!
!!   Grid_pfftGetIndexLimits(integer(OUT) :: configLimits(LOW:HIGH,MDIM),
!!                           integer(OUT) :: phaseLimits(LOW:HIGH,MDIM))
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
#include "Flash.h"
  use gr_pfftData, ONLY : pfft_inLen, pfft_outLen,pfft_me

  implicit none
  integer,dimension(LOW:HIGH,MDIM),intent(OUT) :: configLimits, phaseLimits

  configLimits=1
  configLimits(LOW,1:NDIM)=pfft_me(1:NDIM)*pfft_inLen(1:NDIM) + 1
  configLimits(HIGH,1:NDIM)=configLimits(LOW,1:NDIM) + pfft_inLen(1:NDIM) - 1

  phaseLimits=1
  phaseLimits(LOW,1:NDIM)=pfft_me(1:NDIM)*pfft_outLen(1:NDIM) + 1
  phaseLimits(HIGH,1:NDIM)=phaseLimits(LOW,1:NDIM) + pfft_outLen(1:NDIM) - 1
  return
end subroutine Grid_pfftGetIndexLimits
