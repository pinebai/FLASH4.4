!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonsTrace/pi_minConfinement3DRec
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

  use pi_traceODEfunctionData, ONLY :  xminCell, yminCell, zminCell

  implicit none

  integer, intent (in) :: nc               ! it is absolutely mandatory to
  real,    intent (in) :: y (:)            ! declare the variables and
                                           ! array function in the way shown.
  real :: pi_minConfinement3DRec (1:nc)    ! (compatible to Runge Kutta interfaces)
!
!
!     ...Set up the confinement vector.
!
!
  pi_minConfinement3DRec (1) = xminCell
  pi_minConfinement3DRec (2) = yminCell
  pi_minConfinement3DRec (3) = zminCell
!
!
!     ...Ready!
!
!
  return
end function pi_minConfinement3DRec
