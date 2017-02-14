!!****f* source/Particles/Particles_getGlobalNum
!!
!! NAME
!!  Particles_getGlobalNum
!!
!! SYNOPSIS
!!
!!  Particles_getGlobalNum(integer(OUT)  :: globalNumParticles)
!!                
!! DESCRIPTION 
!!
!!  Returns the global (total) number of particles in the 
!!  simulation after computing it by calling an MPI reduction routine.
!!
!!  Needed for input/output across all processors.
!!
!! ARGUMENTS
!!
!!  globalNumParticles :   integer        number of particles across all processors
!!
!! SIDE EFFECTS
!!
!!  The default implementation writes to the log file a line with the current global
!!  number of particles when first called, and then whenever the number has changed
!!  since the previous time.
!!
!!***

subroutine Particles_getGlobalNum(globalNumParticles)

implicit none
  integer, intent(out)  :: globalNumParticles

  !return 0 particles for stub level
  globalNumParticles = 0

  return
end subroutine Particles_getGlobalNum
