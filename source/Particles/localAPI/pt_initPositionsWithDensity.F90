!!****if* source/Particles/localAPI/pt_initPositionsWithDensity
!!
!! NAME
!!
!!  pt_initPositionsWithDensity
!!
!! SYNOPSIS
!!
!!  call pt_initPositionsWithDensity(integer(in)  :: blockID,
!!                                   logical(OUT) :: success)
!!
!! DESCRIPTION
!!
!!    Initialize particle locations.  This version sets up passive tracer
!!      particles that are distributed randomly according to the gas density     
!!
!! ARGUMENTS
!!
!!   blockID : ID of block in current processor
!!  success:   returns .TRUE. if positions for all particles
!!             that should be assigned to this block have been
!!             successfully initialized.
!!
!! PARAMETERS
!!
!!  pt_numParticlesWanted:   integer  Approximate number of tracer particles to use
!!                                throughout domain ??
!!
!!***
!========================================================================

subroutine pt_initPositionsWithDensity (blockID,success)

  implicit none

  integer, INTENT(in) :: blockID
  logical, INTENT(out) :: success

  success = .false.
  return

end subroutine pt_initPositionsWithDensity


