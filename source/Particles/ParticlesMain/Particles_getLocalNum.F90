!!****if* source/Particles/ParticlesMain/Particles_getLocalNum
!!
!! NAME
!!  Particles_getLocalNum
!!
!! SYNOPSIS
!!
!!  Particles_getLocalNum(integer(IN)  :: blockID,
!!                        integer(OUT) :: localNumParticles)
!!               
!!  
!! DESCRIPTION 
!!
!!  Returns the local number of particles on a given block
!!     Needed for input/output routines
!!  
!! ARGUMENTS
!!
!!  blockID - block number. !!DEV: IGNORED!
!!  localNumParticles - returned value of local particles
!!
!! NOTES
!!
!!***

subroutine Particles_getLocalNum(blockID, localNumParticles)

  use Particles_data, ONLY : pt_numLocal

  implicit none
  integer, intent(in) :: blockID
  integer, intent(out)  :: localNumParticles

  localNumParticles = pt_numLocal

  return
end subroutine Particles_getLocalNum
