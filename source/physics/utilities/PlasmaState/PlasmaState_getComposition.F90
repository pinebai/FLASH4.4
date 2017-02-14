!!****f* source/physics/utilities/PlasmaState/PlasmaState_getComposition
!!
!! NAME
!!
!!  PlasmaState_getComposition
!!
!! SYNOPSIS
!!
!!  call PlasmaState_getComposition(real(OUT), dimension(:)  :: vecz,
!!                                  real(OUT), dimension(:)  :: veca,
!!                                  real(OUT), dimension(:)  :: vecfrac,
!!                                  integer(OUT)  :: nfrac,
!!                                  real(IN), dimension(:)  :: solnstate)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   vecz :       returns list of element Z values as a one-dimensional array
!!
!!   veca :       returns list of element A values as a one-dimensional array
!!
!!   vecfrac :    returns list of number fractions values as a one-dimensional array.
!!                They should add up to 1.
!!
!!   nfrac :      Numer of valid elements in the arrays that are returned.
!!
!!   solnstate :  solution state in one cell
!!
!!
!!
!!***

subroutine PlasmaState_getComposition(vecZ, vecA, vecFrac, nFrac, solnState)
  implicit none

  real, intent(OUT), dimension(:) :: vecZ, vecA, vecFrac
  integer, intent(OUT)            :: nFrac
  real, intent(IN),  dimension(:) :: solnState

  nFrac = 0
  vecZ (:) = 0.0
  vecA (:) = 0.0
  vecFrac (:) = 0.0

end subroutine PlasmaState_getComposition
