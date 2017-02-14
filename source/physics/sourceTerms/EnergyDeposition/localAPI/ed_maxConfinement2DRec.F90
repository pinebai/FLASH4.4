!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_maxConfinement2DRec
!!
!! NAME
!!
!!  ed_maxConfinement2DRec
!!
!! SYNOPSIS
!!
!!  ed_maxConfinement2DRec (integer, intent (in) :: nc,
!!                          real,    intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Places the maximum confinement rules for 2D rectangular ray tracing Runge Kutta runs
!!  as functions of the dependent variables array.
!!
!! ARGUMENTS
!!
!!  nc : number of confined dependent variables (must be 2 and <= size (y), not checked here!)
!!  y  : dependent variable array (rx,ry,vx,vy,P)
!!
!!***

function ed_maxConfinement2DRec (nc,y)

  implicit none

  integer, intent (in) :: nc               ! it is absolutely mandatory to
  real,    intent (in) :: y (:)            ! declare the variables and
                                           ! array function in the way shown.
  real :: ed_maxConfinement2DRec (1:nc)    ! (compatible to Runge Kutta interfaces)

  ed_maxConfinement2DRec (1:nc) = 0.0

  return
end function ed_maxConfinement2DRec
