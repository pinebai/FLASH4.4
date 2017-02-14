!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_maxConfinement3DRec
!!
!! NAME
!!
!!  ed_maxConfinement3DRec
!!
!! SYNOPSIS
!!
!!  ed_maxConfinement3DRec (integer, intent (in) :: nc,
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
!!  y  : dependent variable array (rx,ry,rz,vx,vy,vz,P)
!!
!!***

function ed_maxConfinement3DRec (nc,y)

  implicit none

  integer, intent (in) :: nc               ! it is absolutely mandatory to
  real,    intent (in) :: y (:)            ! declare the variables and
                                           ! array function in the way shown.
  real :: ed_maxConfinement3DRec (1:nc)    ! (compatible to Runge Kutta interfaces)

  ed_maxConfinement3DRec (1:nc) = 0.0

  return
end function ed_maxConfinement3DRec
