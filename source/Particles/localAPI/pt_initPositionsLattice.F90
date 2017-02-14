!!****if* source/Particles/localAPI/pt_initPositionsLattice
!!
!! NAME
!!    pt_initPositionsLattice
!!
!! SYNOPSIS
!!
!!    call pt_initPositionsLattice( integer(in)  :: blockID,
!!                                  logical(out) :: success)
!!
!! DESCRIPTION
!!    Initialize particle locations.  This version sets up particles
!!      which are evenly distributed in space along the axes.  Distribution
!!      in non-Cartesian coordinate systems is not regular.
!!
!! ARGUMENTS
!!
!!  blockID:        local block ID containing particles to create
!!  success:        returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!! PARAMETERS
!!
!!    pt_numX:      number of particles along physical x-axis of domain
!!    pt_numY:      number of particles along physical y-axis of domain
!!    pt_numZ:      number of particles along physical z-axis of domain
!!    pt_initialXMin, pt_initialXMax:  physical domain to initialize with particles in X
!!    pt_initialYMin, pt_initialYMax:  physical domain to initialize with particles in Y
!!    pt_initialZMin, pt_initialZMax:  physical domain to initialize with particles in Z
!!
!!***


subroutine pt_initPositionsLattice (blockID,success)


  implicit none

  integer, INTENT(in) :: blockID
  logical,intent(OUT) :: success

  success = .FALSE.
  return

!----------------------------------------------------------------------
  
end subroutine pt_initPositionsLattice


