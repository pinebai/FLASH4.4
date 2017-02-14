!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonsTrace/pi_maxConfinement3DRec
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

  use pi_traceODEfunctionData, ONLY :  xmaxCell, ymaxCell, zmaxCell

  implicit none

  integer, intent (in) :: nc               ! it is absolutely mandatory to
  real,    intent (in) :: y (:)            ! declare the variables and
                                           ! array function in the way shown.
  real :: pi_maxConfinement3DRec (1:nc)    ! (compatible to Runge Kutta interfaces)
!
!
!     ...Set up the confinement vector.
!
!
  pi_maxConfinement3DRec (1) = xmaxCell
  pi_maxConfinement3DRec (2) = ymaxCell
  pi_maxConfinement3DRec (3) = zmaxCell
!
!
!     ...Ready!
!
!
  return
end function pi_maxConfinement3DRec
