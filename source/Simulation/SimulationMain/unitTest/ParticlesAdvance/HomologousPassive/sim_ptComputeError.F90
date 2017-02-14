!!****if* source/Simulation/SimulationMain/unitTest/ParticlesAdvance/HomologousPassive/sim_ptComputeError
!!
!! NAME
!!
!!  sim_ptComputeError
!!
!! SYNOPSIS
!!
!!  sim_ptComputeError(real(in) :: dtOld,
!!                           real(in) :: dtNew,
!!                           real(in) :: t)
!!
!! DESCRIPTION
!!
!!  Compute error, that is, compare actual to analytic solution.
!!
!!  Compares particles' POS{X,Y,Z}_PART_PROP and POSANAL{X,Y,Z}_PART_PROP
!!  properties.
!!
!! ARGUMENTS
!!
!!   dtOld -- previous time increment
!!   dtNew -- current time increment
!!   t     -- time at which solutions are compared
!!  
!!***

!===============================================================================

subroutine sim_ptComputeError (dtOld,dtNew,t)
    
  use Particles_data, ONLY: particles, pt_numLocal, useParticles
  
  implicit none

#include "Flash.h"
  
  real, INTENT(in)  :: dtOld, dtNew, t
  integer       :: i

!!------------------------------------------------------------------------------
  
  ! Don't do anything if runtime parameter isn't set
!!$      if (.not.useParticles ) return


  ! Update the particle positions.
  do i = 1, pt_numLocal
     particles(ERRX_PART_PROP,i) = particles(POSX_PART_PROP,i) - particles(POSANALX_PART_PROP,i)
     particles(ERRY_PART_PROP,i) = particles(POSY_PART_PROP,i) - particles(POSANALY_PART_PROP,i)
     particles(ERRZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) - particles(POSANALZ_PART_PROP,i)
     if (abs(particles(ERRX_PART_PROP,i)) > particles(MAXERRX_PART_PROP,i)) then
        particles(MAXERRX_PART_PROP,i) = abs(particles(ERRX_PART_PROP,i))
        particles(ERRWHENX_PART_PROP,i) = t
     end if
     if (abs(particles(ERRY_PART_PROP,i)) > particles(MAXERRY_PART_PROP,i)) then
        particles(MAXERRY_PART_PROP,i) = abs(particles(ERRY_PART_PROP,i))
        particles(ERRWHENY_PART_PROP,i) = t
     end if
     if (abs(particles(ERRZ_PART_PROP,i)) > particles(MAXERRZ_PART_PROP,i)) then
        particles(MAXERRZ_PART_PROP,i) = abs(particles(ERRZ_PART_PROP,i))
        particles(ERRWHENZ_PART_PROP,i) = t
     end if
  enddo
     
  return
!!------------------------------------------------------------------------------
  
end subroutine sim_ptComputeError


