!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_raytraceODEfunction2DCyl3D
!!
!! NAME
!!
!!  ed_raytraceODEfunction2DCyl3D
!!
!! SYNOPSIS
!!
!!  ed_raytraceODEfunction2DCyl3D (real, intent (in) :: t,
!!                                 real, intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Defines the ray trace ODE function for the Runge Kutta unit, to perform 3D ray tracing in
!!  2D cylindrical coordinates.
!!
!! ARGUMENTS
!!
!!  t    : independent variable (time)
!!  y    : dependent variable array (rx,ry,rz,vx,vy,vz,P)
!!
!!***

function ed_raytraceODEfunction2DCyl3D (t,y)

  implicit none

  real, intent (in) :: t                             ! it is absolutely mandatory to
  real, intent (in) :: y (:)                         ! declare the variables and
                                                     ! array function in the way shown.
  real :: ed_raytraceODEfunction2DCyl3D (1:size (y)) ! (compatible to Runge Kutta interfaces)

  ed_raytraceODEfunction2DCyl3D (1:size (y)) = 0.0

  return
end function ed_raytraceODEfunction2DCyl3D
