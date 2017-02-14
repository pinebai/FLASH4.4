!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/RK/ed_maxConfinement1DRec
!!
!! NAME
!!
!!  ed_maxConfinement1DRec
!!
!! SYNOPSIS
!!
!!  ed_maxConfinement1DRec (integer, intent (in) :: nc,
!!                          real,    intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Places the maximum confinement rules for 1D rectangular ray tracing Runge Kutta runs
!!  as functions of the dependent variables array.
!!
!! ARGUMENTS
!!
!!  nc : number of confined dependent variables (must be 1 and <= size (y), not checked here!)
!!  y  : dependent variable array (rx,vx,P)
!!
!!***

function ed_maxConfinement1DRec (nc,y)

  use ed_raytraceODEfunctionData, ONLY :  xmaxCell

  implicit none

  integer, intent (in) :: nc               ! it is absolutely mandatory to
  real,    intent (in) :: y (:)            ! declare the variables and
                                           ! array function in the way shown.
  real :: ed_maxConfinement1DRec (1:nc)    ! (compatible to Runge Kutta interfaces)
!
!
!     ...Set up the confinement vector.
!
!
  ed_maxConfinement1DRec (1) = xmaxCell
!
!
!     ...Ready!
!
!
  return
end function ed_maxConfinement1DRec
