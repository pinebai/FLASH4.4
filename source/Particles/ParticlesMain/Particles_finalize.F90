!!****if* source/Particles/ParticlesMain/Particles_finalize
!!
!! NAME
!!    Particles_finalize
!!
!! SYNOPSIS
!!    Particles_finalize( )
!!
!! DESCRIPTION
!!
!!    Finalize routine for the particle unit.  Removes memory usage
!!      set in Particles_init.
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!
!!***
  
subroutine Particles_finalize()
  use Particles_data, ONLY : particles,useParticles
  use pt_interface, ONLY : pt_initFinalize
  implicit none
  if(useParticles) deallocate(particles)
  call pt_initFinalize()
end subroutine Particles_finalize
