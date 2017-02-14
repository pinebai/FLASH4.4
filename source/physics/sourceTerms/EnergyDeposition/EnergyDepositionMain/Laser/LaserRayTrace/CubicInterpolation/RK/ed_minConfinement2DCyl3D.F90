!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/RK/ed_minConfinement2DCyl3D
!!
!! NAME
!!
!!  ed_minConfinement2DCyl3D
!!
!! SYNOPSIS
!!
!!  ed_minConfinement2DCyl3D (integer, intent (in) :: nc,
!!                            real,    intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Places the minimum confinement rules for 3D in 2D cylindrical ray tracing Runge Kutta runs
!!  as functions of the dependent variables array.
!!
!! ARGUMENTS
!!
!!  nc : number of confined dependent variables (must be 3 and <= size (y), not checked here!)
!!  y  : dependent variable array (rx,ry,rz,vx,vy,vz,P)
!!
!!***

function ed_minConfinement2DCyl3D (nc,y)

  use ed_raytraceODEfunctionData, ONLY : wedgeSlope,           &
                                         xminWedge, zminWedge

  implicit none

  integer, intent (in) :: nc               ! it is absolutely mandatory to
  real,    intent (in) :: y (:)            ! declare the variables and
                                           ! array function in the way shown.
  real :: ed_minConfinement2DCyl3D (1:nc)  ! (compatible to Runge Kutta interfaces)
!
!
!     ...Set up the confinement vector.
!
!
  ed_minConfinement2DCyl3D (1) = xminWedge
  ed_minConfinement2DCyl3D (2) = - wedgeSlope * y (1)       ! y (1) contains rayX
  ed_minConfinement2DCyl3D (3) = zminWedge
!
!
!     ...Ready!
!
!
  return
end function ed_minConfinement2DCyl3D
