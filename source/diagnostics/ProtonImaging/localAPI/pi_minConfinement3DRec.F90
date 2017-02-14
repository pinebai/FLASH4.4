!!****if* source/diagnostics/ProtonImaging/localAPI/pi_minConfinement3DRec
!!
!! NAME
!!
!!  pi_minConfinement3DRec
!!
!! SYNOPSIS
!!
!!  pi_minConfinement3DRec (integer, intent (in) :: nc,
!!                          real,    intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Places the minimum confinement rules for 3D rectangular ray tracing Runge Kutta runs
!!  as functions of the dependent variables array.
!!
!! ARGUMENTS
!!
!!  nc : number of confined dependent variables (must be 3 and <= size (y), not checked here!)
!!  y  : dependent variable array (not used here)
!!
!!***

function pi_minConfinement3DRec (nc,y)

  implicit none

  integer, intent (in) :: nc               ! it is absolutely mandatory to
  real,    intent (in) :: y (:)            ! declare the variables and
                                           ! array function in the way shown.
  real :: pi_minConfinement3DRec (1:nc)    ! (compatible to Runge Kutta interfaces)

  pi_minConfinement3DRec (1:nc) = 0.0

  return
end function pi_minConfinement3DRec
