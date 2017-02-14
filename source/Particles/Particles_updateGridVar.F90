!!****f* source/Particles/Particles_updateGridVar
!!
!! NAME
!!    Particles_updateGridVar
!!
!! SYNOPSIS
!!    Particles_updateGridVar(integer(IN)    :: partProp,
!!                            integer(IN)    :: varGrid,
!!                            integer(IN),optional    :: mode)
!!
!! DESCRIPTION
!!
!!    Updates the given grid variable with data from the given particle
!!    property.
!!
!! ARGUMENTS
!!               partProp:  Index of particle attribute to interpolate onto 
!!                          the mesh
!!               varGrid:   Index of gridded variable to receive interpolated
!!                          quantity
!!               mode:      (Optional) If zero (default), zero varGrid first;
!!                          if present and nonzero, do not zero varGrid first
!!                          but add data from particles to existing grid data.
!!
!! PARAMETERS
!! 
!!***
  
subroutine Particles_updateGridVar(partProp, varGrid, mode)

  implicit none

  integer, INTENT(in) :: partProp, varGrid
  integer, INTENT(in), optional :: mode

  return

end subroutine Particles_updateGridVar
