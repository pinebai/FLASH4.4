!!****f* source/Particles/Particles_putLocalNum
!!
!! NAME
!!  Particles_putLocalNum
!!
!! SYNOPSIS
!!
!!  Particles_putLocalNum(integer(IN)  :: localNumParticles)
!!               
!!  
!! DESCRIPTION 
!!
!!  Sets the local number of particles on a given proc
!!  
!!
!! ARGUMENTS
!!
!!  localNumParticles:  the local number of particles
!!
!! NOTES
!!
!!
!!***

subroutine Particles_putLocalNum(localNumParticles)

implicit none
  integer, intent(in)  :: localNumParticles

  return
end subroutine Particles_putLocalNum
