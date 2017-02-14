!!****if* source/Particles/ParticlesMain/Particles_putLocalNum
!!
!! NAME
!!  Particles_putLocalNum
!!
!! SYNOPSIS
!!
!!  Particles_putLocalNum(integer(IN)  :: localNumParticles)
!!                        
!!               
!!  
!! DESCRIPTION 
!!
!!  Sets the number of particles local to the current MPI task
!!  
!!
!! ARGUMENTS
!!
!!  localNumParticles - input value of local number of particles
!!
!! NOTES
!!
!!
!!***

subroutine Particles_putLocalNum(localNumParticles)

  use Particles_data, ONLY : pt_numLocal

implicit none
  integer, intent(in)  :: localNumParticles

  pt_numLocal = localNumParticles


  return
end subroutine Particles_putLocalNum
