!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_gridLineRegular
!!
!! NAME
!!
!!  ed_gridLineRegular
!!
!! SYNOPSIS
!!
!!  call ed_gridLineRegular (real    (in)  :: lineLength,
!!                           integer (in)  :: nGridPoints,
!!                           real    (out) :: delta)
!!
!! DESCRIPTION
!!
!!  Given a line defined by its length, the routine places a regular line grid onto it,
!!  such that the number of tics on the line matches the number of grid points wanted.
!!  The tics cover the entire line and for a regular grid the separation between consecutive
!!  tics (delta) is constant. The number of grid points can always be honored by this type
!!  of grid. Below is a picture of the grid:
!!
!!
!!                        |---|---|---|---|-+-|---|---|---|---|
!!
!!                        |   <-----   line length   ----->   |
!!
!! ARGUMENTS
!!
!!  lineLength  : the length of the line
!!  nGridPoints : the number of grid points wanted
!!  delta       : the separation distance between two consecutive grid points
!!
!! NOTES
!!
!!***

subroutine ed_gridLineRegular (lineLength, nGridPoints,    delta)

  implicit none

  real,    intent (in)  :: lineLength
  integer, intent (in)  :: nGridPoints
  real,    intent (out) :: delta

  delta = 0.0

  return
end subroutine ed_gridLineRegular
