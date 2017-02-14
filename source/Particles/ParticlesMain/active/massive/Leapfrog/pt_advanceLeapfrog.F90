!!****if* source/Particles/ParticlesMain/active/massive/Leapfrog/pt_advanceLeapfrog
!!
!! NAME
!!
!!  pt_advanceLeapfrog
!!
!! SYNOPSIS
!!
!!  call pt_advanceLeapfrog(real(IN) :: dtOld,
!!                  real(IN) :: dtNew,
!!                  real(inout):: particles(:,p_count),
!!                  integer(in):: p_count,
!!                  integer(in):: ind)
!!
!! ARGUMENTS
!!  
!!   dtOld : previous time interval 
!!   dtNew : current time interval
!!   particles -- particles to advance
!!   p_count  -- the number of particles in the list to advance
!!   ind -- index into pt_typeInfo
!!  
!! DESCRIPTION
!!
!!  Time advancement routine for the Particles unit.
!!  
!!  This version is the forward leapfrog advancement for the active 
!!  submodule.
!!
!!***

!======================================================================

subroutine pt_advanceLeapfrog (dtOld,dtNew,particles,p_count, ind)
    
  use Particles_data, ONLY: useParticles, pt_meshMe, pt_restart
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Particles_interface, ONLY : Particles_longRangeForce, Particles_shortRangeForce
  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"
    
  real, INTENT(in) :: dtOld, dtNew
  integer,intent(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(inout) :: particles

  integer       :: i
  real          :: wterm, woldterm
  real, parameter :: onethird = 1./3., onesixth = 1./6.
  logical, save :: firstCall = .true.
  integer,parameter :: mapType=WEIGHTED

!-------------------------------------------------------------------------------

  if (.not. useParticles) return

  if (firstCall .and. .not. pt_restart) then
     wterm = 0.5*dtNew
     woldterm = 0.
     firstCall = .false.
  else
     wterm = 0.5*dtNew + onethird*dtOld + onesixth*dtNew**2/dtOld
     woldterm = onesixth*(dtOld**2 - dtNew**2)/dtOld
  endif

  call Timers_start ("particles")

!-------------------------------------------------------------------------------
    
! Compute current components of the acceleration of each particle


  call Timers_start ("particle forces")

  do i = 1, p_count
     particles(ACCX_PART_PROP,i) = 0.
     particles(ACCY_PART_PROP,i) = 0.  
     particles(ACCZ_PART_PROP,i) = 0.
  enddo

  call Particles_longRangeForce(particles,p_count,mapType)
  call Particles_shortRangeForce

  call Timers_stop ("particle forces")

!-------------------------------------------------------------------------------

! Apply the Euler step to update particle positions and velocities

  call Timers_start ("move particles")

  do i = 1, p_count


     particles(VELX_PART_PROP,i) = particles(VELX_PART_PROP,i) + &
          wterm*particles(ACCX_PART_PROP,i) + &
          woldterm*particles(OACX_PART_PROP,i)

     particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + &
          dtNew*particles(VELX_PART_PROP,i)

     particles(OACX_PART_PROP,i) = particles(ACCX_PART_PROP,i)
      
     if (NDIM >= 2) then

        particles(VELY_PART_PROP,i) = particles(VELY_PART_PROP,i) + &
          wterm*particles(ACCY_PART_PROP,i) + &
          woldterm*particles(OACY_PART_PROP,i)

        particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + &
             dtNew*particles(VELY_PART_PROP,i)

        particles(OACY_PART_PROP,i) = particles(ACCY_PART_PROP,i)
     else
        particles(VELY_PART_PROP,i)=0.0
     endif

     if (NDIM == 3) then

        particles(VELZ_PART_PROP,i) = particles(VELZ_PART_PROP,i) + &
          wterm*particles(ACCZ_PART_PROP,i) + &
          woldterm*particles(OACZ_PART_PROP,i)

        particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + &
             dtNew*particles(VELZ_PART_PROP,i)
        
        particles(OACZ_PART_PROP,i) = particles(ACCZ_PART_PROP,i)
     else
        particles(VELZ_PART_PROP,i)=0.0
     endif
       
  enddo
  call Timers_stop ("move particles")

!-------------------------------------------------------------------------------

  call Timers_stop ("particles")

  return

end subroutine pt_advanceLeapfrog    

!===============================================================================

