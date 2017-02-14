!!****if* source/Particles/ParticlesMain/active/massive/LeapfrogCosmo/pt_advanceLeapfrog_cosmo
!!
!! NAME
!!
!!  pt_advanceLeapfrog_cosmo
!!
!! SYNOPSIS
!!
!!  call pt_advanceLeapfrog_cosmo(real(IN) :: dtOld,
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
!!   ind   -- index into pt_typeInfo
!!  
!! DESCRIPTION
!!
!!  Time advancement routine for the Particles unit.
!!  
!!  This version is the second-order leapfrog advancement for the active 
!!  submodule, including cosmological redshift terms.
!!
!!  The method used is that worked out by Scott Dodelson for COSMOS
!!  (2000, ApJ, 536, 122).
!!
!!***

!======================================================================

subroutine pt_advanceLeapfrog_cosmo (dtOld,dtNew,particles,p_count, ind)
    
  use Particles_data, ONLY: useParticles, pt_meshMe, pt_restart
  use Particles_interface, ONLY : Particles_longRangeForce, &
    Particles_shortRangeForce
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_getSimTime
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Cosmology_interface, ONLY : Cosmology_getOldRedshift, &
       Cosmology_getParams

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"

  real, INTENT(in) :: dtOld, dtNew
  integer,intent(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(inout) :: particles

  integer       :: i
  real          :: t, zOld, sOld
  logical, save :: firstCall = .true.
  real          :: alpha, dotalpha, aterm, bterm, wterm, woldterm
  real, save    :: Hubble, OmegaM, OmegaCurv, OmegaL, OmegaB
  real :: BNcoeff, CNcoeff, DNcoeff  
  integer,parameter :: mapType=WEIGHTED

  !-------------------------------------------------------------------------------
  
  if (.not. useParticles) return
  
  call Timers_start ("particles")
  
  !-------------------------------------------------------------------------------
  
  ! Initialize needed constants on the first call.
  
  if (firstCall) then
     call Cosmology_getParams(Hubble,OmegaM,OmegaB,OmegaL)
     OmegaCurv = OmegaM + OmegaL - 1.
  endif

  ! Get time, timestep, scale factor, and rate of change of scale factor.
  
  call Driver_getSimTime(t)
  call Cosmology_getOldRedshift(zOld)
  
  t        = t - dtNew
  sOld = 1./(1.+zOld)
  
  alpha    = 2*Hubble/sOld * sqrt(OmegaM/sOld-OmegaCurv+OmegaL*sOld**2)
  dotalpha = 2*Hubble**2 * (-3*OmegaM + 2*OmegaCurv*sOld)/sOld**3
  
  ! For the first step, velocity is defined at t_0, not t_{-1/2}, so we take
  ! a first-order step to begin.  Afterward velocity and position are staggered
  ! in time.
    
  if (firstCall .and. .not. pt_restart) then
     aterm      = dtNew
     bterm      = 1. - alpha*dtNew*0.5
     wterm      = 0.5*dtNew
     woldterm   = 0.
     firstCall = .false.
  else
     aterm      = dtNew
     bterm      = BNcoeff(dtNew, dtOld, alpha, dotalpha)
     wterm      = CNcoeff(dtNew, dtOld, alpha, dotalpha)
     woldterm   = DNcoeff(dtNew, dtOld, alpha, dotalpha)
  endif
  
  !-------------------------------------------------------------------------------
  
  ! Compute current components of the acceleration of each particle
  
  call Timers_start ("particle forces")
  
  
  do i = 1, p_count
     particles(ACCX_PART_PROP,i) = 0.
     particles(ACCY_PART_PROP,i) = 0.  
     particles(ACCZ_PART_PROP,i) = 0.    
  enddo
  
  call Particles_longRangeForce(particles,p_count,mapType)
  call Particles_shortRangeForce()
  
  call Timers_stop ("particle forces")
  
  !-------------------------------------------------------------------------------
  
  ! Apply the leapfrog step to update particle positions and velocities
  ! Save the particle accelerations for use in the next timestep
  
  call Timers_start ("move particles")
  
  do i = 1, p_count
     
     particles(VELX_PART_PROP,i) = bterm*particles(VELX_PART_PROP,i) + &
          wterm*particles(ACCX_PART_PROP,i) + &
          woldterm*particles(OACX_PART_PROP,i)
     particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + &
          aterm*particles(VELX_PART_PROP,i)
     particles(OACX_PART_PROP,i) = particles(ACCX_PART_PROP,i)
     
     if (NDIM >= 2) then
        particles(VELY_PART_PROP,i)= bterm*particles(VELY_PART_PROP,i) + &
             wterm*particles(ACCY_PART_PROP,i) + &
             woldterm*particles(OACY_PART_PROP,i)
        particles(POSY_PART_PROP,i)=  particles(POSY_PART_PROP,i) + &
             aterm*particles(VELY_PART_PROP,i)
        particles(OACY_PART_PROP,i) = particles(ACCY_PART_PROP,i)
     endif
     if (NDIM == 3) then
        particles(VELZ_PART_PROP,i)= bterm*particles(VELZ_PART_PROP,i) + &
             wterm*particles(ACCZ_PART_PROP,i) + &
             woldterm*particles(OACZ_PART_PROP,i)
        particles(POSZ_PART_PROP,i) =  particles(POSZ_PART_PROP,i) + &
             aterm*particles(VELZ_PART_PROP,i)
        particles(OACZ_PART_PROP,i) = particles(ACCZ_PART_PROP,i)
     endif
     
  enddo
  
  call Timers_stop ("move particles")
  
  !-------------------------------------------------------------------------------
  
  call Timers_stop ("particles")
  
  return
  
end subroutine pt_advanceLeapfrog_cosmo    

