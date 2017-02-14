!!****if* source/Simulation/SimulationMain/unitTest/ParticlesAdvance/sim_getComputedError
!!
!! NAME
!!
!!  sim_getComputedError
!!
!! SYNOPSIS
!!
!!  sim_getComputedError(real(out) :: error)
!!
!! DESCRIPTION
!!
!!  Returns maximum Computed error from comparison of actual to analytic solution.
!!
!! ARGUMENTS
!!
!!   error -- them maximum error
!!  
!! SEE ALSO
!!
!!  sim_ptComputeError
!!***

!===============================================================================

subroutine sim_getComputedError (error)
    
  use Particles_data, ONLY: particles, pt_numLocal, useParticles
  
  implicit none

#include "Flash.h"
  
  real, INTENT(out)  :: error

!!------------------------------------------------------------------------------
  
  ! Don't do anything if runtime parameter isn't set
!!$      if (.not.useParticles ) return


  error = maxval(particles(MAXERRX_PART_PROP,1:pt_numLocal))
  if (NDIM > 1) then
     error = max(error,maxval(particles(MAXERRY_PART_PROP,1:pt_numLocal)))
  end if
  if (NDIM > 2) then
     error = max(error,maxval(particles(MAXERRZ_PART_PROP,1:pt_numLocal)))
  end if
     
  return
!!------------------------------------------------------------------------------
  
end subroutine sim_getComputedError


