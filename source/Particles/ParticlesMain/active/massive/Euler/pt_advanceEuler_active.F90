!!****if* source/Particles/ParticlesMain/active/massive/Euler/pt_advanceEuler_active
!!
!! NAME
!!
!!  pt_advanceEuler_active
!!
!! SYNOPSIS
!!
!!  pt_advanceEuler_active(real(IN) :: dtOld,
!!                  real(IN) :: dtNew,
!!                  real(inout):: particles(:,p_count),
!!                  integer(in):: p_count)
!!                  integer(in):: ind)
!!
!! ARGUMENTS
!!  
!!   dtOld : previous time interval 
!!   dtNew : current time interval
!!   particles -- particles to advance
!!   p_count  -- the number of particles in the list to advance
!!   ind    --- index for type into pt_typeInfo
!!  
!! DESCRIPTION
!!
!!  Time advancement routine for the particle module.
!!  
!!  This version is the forward Euler advancement for the active 
!!  submodule.
!!
!!***

!===============================================================================

subroutine pt_advanceEuler_active (dtOld,dtNew,particles,p_count, ind)
    
  use Particles_data, ONLY: useParticles
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Particles_interface, ONLY : Particles_longRangeForce, Particles_shortRangeForce
expensive
  implicit none

#include "Flash.h"
#include "constants.h"
    
  integer       :: i
  real, INTENT(in) :: dtOld, dtNew
  integer,intent(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(inout) :: particles
!-------------------------------------------------------------------------------
  integer :: mapType=WEIGHTED

  if (.not.useParticles ) return

  call Timers_start ("particles")

! Compute current components of the acceleration of each particle

  call Timers_start ("particle forces")

  do i = 1, p_count
     particles(ACCX_PART_PROP,i) = 0.
     particles(ACCY_PART_PROP,i) = 0.  
     particles(ACCZ_PART_PROP,i) = 0.
  enddo

  call Particles_longRangeForce(particle,p_count,mapType)
  call Particles_shortRangeForce

  call Timers_stop ("particle forces")

!-------------------------------------------------------------------------------

! Apply the Euler step to update particle positions and velocities

  call Timers_start ("move particles")
    
  do i = 1, p_count

     particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + &
          dtOld*particles(VELX_PART_PROP,i)

     particles(VELX_PART_PROP,i) = particles(VELX_PART_PROP,i) + &
          dtOld*particles(ACCX_PART_PROP,i)
      
     if (NDIM >= 2) then
        particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + &
             dtOld*particles(VELY_PART_PROP,i)

        particles(VELY_PART_PROP,i) = particles(VELY_PART_PROP,i) + &
             dtOld*particles(ACCY_PART_PROP,i)
     endif

     if (NDIM == 3) then
        particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + &
             dtOld*particles(VELZ_PART_PROP,i)
        
        particles(VELZ_PART_PROP,i) = particles(VELZ_PART_PROP,i) + &
             dtOld*particles(ACCZ_PART_PROP,i)
     endif
       
  enddo


  call Timers_stop ("move particles")

!-------------------------------------------------------------------------------

  call Timers_stop ("particles")

  return

end subroutine pt_advanceEuler_active

!===============================================================================

