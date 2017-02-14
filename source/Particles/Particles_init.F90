!!****f* source/Particles/Particles_init
!!
!! NAME
!!    Particles_init
!!
!! SYNOPSIS
!!    Particles_init( logical(in) :: restart )
!!
!! DESCRIPTION
!!
!!    General initialization routine for the particle module.
!!
!! ARGUMENTS
!!
!!    restart:   indicates if run is starting from scratch or restarting
!!               from checkpoint file
!!
!! PARAMETERS
!!
!!    useParticles   BOOLEAN [TRUE]  Should particles be used in this simulation?
!!    pt_maxPerProc  INTEGER [100]   Maximum number of particles per processor. Allocates array space
!!                                   Particles are distributed per PROCESSOR rather than per BLOCK
!!    pt_dtFactor    REAL    [0.5]   Factor to make sure that time step is small enough that particles
!!    pt_dtChangeTolerance REAL [0.4] For uncorrected Estimated Midpoint propagation scheme:
!!                                    Do Euler step if change in time step is greater than this
!!                                    percentage.  Set to 0 to always do Euler, set to a huge
!!                                    number to always use estimated midpoint velocities
!!    pt_small       REAL    [1.0E-10] Used for general comparisons of real values 
!!                                   For example, IF (abs(real1 - real2) .lt. pt_small) THEN
!!                                   don't move farther than one block in each step
!!
!!***
  
subroutine Particles_init( restart)

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  logical, INTENT(in) :: restart

  logical, save :: testUseParticles

  !! It is a failure to invoke the stub when useParticles is set TRUE.

  call RuntimeParameters_get ("useParticles", testUseParticles)
  if (testUseParticles) then
     call Driver_abortFlash("Particles unit seems not to be compiled in, and the Particles_init stub does not &
          &allow the value of useParticles to be TRUE.")
  end if

end subroutine Particles_init

