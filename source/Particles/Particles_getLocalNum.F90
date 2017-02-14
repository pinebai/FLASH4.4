!!****f* source/Particles/Particles_getLocalNum
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
!!  blockID - block number !!DEV: IGNORED!
!!  localNumParticles - returned value of local particles
!!
!! NOTES
!!
!!
!!***

subroutine Particles_getLocalNum(blockID, localNumParticles)

implicit none
  integer, intent(in) :: blockID
  integer, intent(out)  :: localNumParticles

  !return 0 particles for stub level
  localNumParticles = 0

  return
end subroutine Particles_getLocalNum
