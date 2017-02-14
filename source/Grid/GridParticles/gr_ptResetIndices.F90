!!****if* source/Grid/GridParticles/gr_ptResetIndices
!!
!! NAME
!!
!!  gr_ptResetIndices
!!
!! SYNOPSIS
!!
!!  gr_ptResetIndices(integer(IN) :: index_list(:),
!!                  integer(IN) :: count)
!!
!! DESCRIPTION
!!
!!  GridParticles subunit is used for both, tracking the Lagrangian particles, 
!!  and for tracing rays. The attributes for both data structures, though similar,
!!  do differ some. So the routine within GridParticles cannot make assumptions about
!!  knowing the indices such as positions, tag etc. This routine resets all indiced to 
!!  zero. It should be called whenever the control goes out of scope for the indices
!!  set by gr_ptSetIndices to make sure that overlooked values in the indices, which 
!!  may have wrong but valid values for the current scope are not used.
!!  
!!
!! ARGUMENTS
!!
!! index_list - the list of indices for attributes within the data structures that
!!              need to be known by the GridParticles subunit
!! count      - the count of entries in the list carried by index_list
!!
!!***
#include "GridParticles.h"

subroutine gr_ptResetIndices(index_list,count)

  use gr_ptData, ONLY : gr_ptBlk, gr_ptProc, gr_ptTag, &
       gr_ptPosx, gr_ptPosy, gr_ptPosz,&
       gr_ptPos2x, gr_ptPos2y, gr_ptPos2z,&
       gr_ptVelx, gr_ptVely, gr_ptVelz,&
       gr_ptVel2x, gr_ptVel2y, gr_ptVel2z

  implicit none
  integer, intent(IN) :: count
  integer,dimension(count), intent(IN) :: index_list

  gr_ptPosx = GRPT_RESET
  gr_ptPosy = GRPT_RESET
  gr_ptPosz = GRPT_RESET
  gr_ptPos2x = GRPT_RESET
  gr_ptPos2y = GRPT_RESET
  gr_ptPos2z = GRPT_RESET
  gr_ptVelx = GRPT_RESET 
  gr_ptVely = GRPT_RESET 
  gr_ptVelz = GRPT_RESET 
  gr_ptVel2x = GRPT_RESET
  gr_ptVel2y = GRPT_RESET
  gr_ptVel2z = GRPT_RESET
  gr_ptBlk  = GRPT_RESET 
  gr_ptProc = GRPT_RESET 
  gr_ptTag  = GRPT_RESET 

end subroutine gr_ptResetIndices
