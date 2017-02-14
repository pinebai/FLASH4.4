!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/RK/ed_minConfinement3DRec
!!
!! NAME
!!
!!  ed_minConfinement3DRec
!!
!! SYNOPSIS
!!
!!  ed_minConfinement3DRec (integer, intent (in) :: nc,
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
!!  y  : dependent variable array (rx,ry,rz,vx,vy,vz,P)
!!
!!***

function ed_minConfinement3DRec (nc,y)

  use ed_raytraceODEfunctionData, ONLY :  xminCell, yminCell, zminCell

  implicit none

  integer, intent (in) :: nc              ! it is absolutely mandatory to
  real,    intent (in) :: y (:)           ! declare the variables and
                                          ! array function in the way shown.
  real :: ed_minConfinement3DRec (1:nc)   ! (compatible to Runge Kutta interfaces)
!
!
!     ...Set up the confinement vector.
!
!
  ed_minConfinement3DRec (1) = xminCell
  ed_minConfinement3DRec (2) = yminCell
  ed_minConfinement3DRec (3) = zminCell
!
!
!     ...Ready!
!
!
  return
end function ed_minConfinement3DRec
