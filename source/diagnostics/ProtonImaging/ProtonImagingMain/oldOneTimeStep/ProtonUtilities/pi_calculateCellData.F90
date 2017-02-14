!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonUtilities/pi_calculateCellData
!!
!! NAME
!!
!!  pi_calculateCellData
!!
!! SYNOPSIS
!!
!!  pi_calculateCellData (real    (in)  :: coordMin,
!!                        real    (in)  :: coordMax,
!!                        integer (in)  :: indexMin,
!!                        integer (in)  :: indexMax,
!!                        real    (out) :: delta,
!!                        real    (out) :: coords (:))
!!
!! DESCRIPTION
!!
!!  This routine calculates the cell delta and all cell coordinates based on a minimum and
!!  maximum coordinate value and the number of cells along that coordinate:
!!
!!                              cell
!!                            |-----|-----|-----|-----|-----|-----|
!!                               |     |     |
!!                              1st   2nd   3rd
!!                             index index index
!!                               |                             |
!!                           indexMin                      indexMax
!!
!!                            |     |     |
!!                           1st   2nd   3rd
!!                          coord coord coord
!!
!!
!!  Note that the number of cell coordinates equals the number of cell indices + 1. Cells
!!  are assumed to be of equal size (delta). In order to strictly match the endpoint cell
!!  coordinates to the given minimum and maximum coordinates, the endpoint cell coordinates
!!  are set equal to those minimum and maximum coordinates.
!!
!! ARGUMENTS
!!
!!  coordMin   : the minimum coordinate value
!!  coordMax   : the maximum coordinate value
!!  indexMin   : the minimum (first) index
!!  indexMax   : the maximum (last) index
!!  delta      : the cell size
!!  coords     : the calculated cell coordinates
!!
!!***

subroutine pi_calculateCellData (coordMin, coordMax, indexMin, indexMax, delta, coords)

  implicit none

  real,    intent (in)  :: coordMin, coordMax
  integer, intent (in)  :: indexMin, indexMax
  real,    intent (out) :: delta
  real,    intent (out) :: coords (indexMin:indexMax+1)

  integer :: i
!
!
!     ...Calculate the cell data.
!
!
  delta = (coordMax - coordMin) / real (indexMax - indexMin + 1)

  do i = indexMin, indexMax + 1
     coords (i) = coordMin + real (i - indexMin) * delta
  end do
!
!
!     ...Ready!
!
!
  return
end subroutine pi_calculateCellData
