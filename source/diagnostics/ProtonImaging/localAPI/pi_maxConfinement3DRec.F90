!!****if* source/diagnostics/ProtonImaging/localAPI/pi_maxConfinement3DRec
!!
!! NAME
!!
!!  pi_maxConfinement3DRec
!!
!! SYNOPSIS
!!
!!  pi_maxConfinement3DRec (integer, intent (in) :: nc,
!!                          real,    intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Places the maximum confinement rules for 3D rectangular ray tracing Runge Kutta runs
!!  as functions of the dependent variables array.
!!
!! ARGUMENTS
!!
!!  nc : number of confined dependent variables (must be 3 and <= size (y), not checked here!)
!!  y  : dependent variable array (not used here)
!!
!!***

function pi_maxConfinement3DRec (nc,y)

  implicit none

  integer, intent (in) :: nc               ! it is absolutely mandatory to
  real,    intent (in) :: y (:)            ! declare the variables and
                                           ! array function in the way shown.
  real :: pi_maxConfinement3DRec (1:nc)    ! (compatible to Runge Kutta interfaces)

  pi_maxConfinement3DRec (1:nc) = 0.0

  return
end function pi_maxConfinement3DRec
