!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_raytraceODEfunction1DRec
!!
!! NAME
!!
!!  ed_raytraceODEfunction1DRec
!!
!! SYNOPSIS
!!
!!  ed_raytraceODEfunction1DRec (real, intent (in) :: t,
!!                               real, intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Defines the ray trace ODE function in 1D rectangular coordinates to be used for Runge Kutta
!!  ray tracing.
!!
!! ARGUMENTS
!!
!!  t    : independent variable (time)
!!  y    : dependent variable array (rx,vx,P)
!!
!!***

function ed_raytraceODEfunction1DRec (t,y)

  implicit none

  real, intent (in) :: t                             ! it is absolutely mandatory to
  real, intent (in) :: y (:)                         ! declare the variables and
                                                     ! array function in the way shown.
  real :: ed_raytraceODEfunction1DRec (1:size (y))   ! (compatible to Runge Kutta interfaces)

  ed_raytraceODEfunction1DRec (1:size (y)) = 0.0

  return
end function ed_raytraceODEfunction1DRec
