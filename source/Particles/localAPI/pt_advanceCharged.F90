!!****if* source/Particles/localAPI/pt_advanceCharged
!!
!! NAME
!!
!!  pt_advanceCharged
!!
!! SYNOPSIS
!!
!!  call pt_advanceCharged(real(in)    :: dtold,
!!                        real(in)    :: dtnew,
!!                        real(inout) :: particlesunused(NPART_PROPS,p_countUnused),
!!                        integer(in) :: p_countunused)
!!
!! DESCRIPTION
!!
!!   Advances particles in time
!!
!! ARGUMENTS
!!
!!   dtOld : previous time interval 
!!   dtNew : current time interval
!!   particlesunused -- particles on which to operate
!!   p_countunused - the number of particles in the list to advance
!!
!!
!!
!!***

subroutine pt_advanceCharged (dtOld,dtNew,particlesUnused,p_countUnused)

  implicit none

#include "Flash.h"
#include "constants.h"  
#include "Particles.h"

  real, intent(in)   :: dtOld, dtNew
  ! Cannot use these arguments since we want to add/delete particles
  integer,intent(in) :: p_countUnused
  real,dimension(NPART_PROPS,p_countUnused),intent(inout) :: particlesUnused

  
end subroutine pt_advanceCharged
