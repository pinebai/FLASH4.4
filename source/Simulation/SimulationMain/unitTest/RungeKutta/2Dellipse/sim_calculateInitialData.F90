!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/2Dellipse/sim_calculateInitialData
!!
!! NAME
!!
!!  sim_calculateInitialData
!! 
!! SYNOPSIS
!!
!!  call sim_calculateInitialData
!!
!! DESCRIPTION
!!
!!  Calculates the positive 'k' value for the elliptical 2D ODE equation from the given aspect
!!  ratio A:
!!
!!                 k = (A^2 + 1) / (A^2 - 1)
!!
!!  The routine also gives the parameters characterizing the elliptical path that will be traced.
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_calculateInitialData ()

  use  Simulation_data, ONLY: sim_ellipseMajorSemiAxis, &
                              sim_ellipseMinorSemiAxis, &
                              sim_ellipseRotationAngle, &
                              sim_ellipseCenterX,       &
                              sim_ellipseCenterY,       &
                              sim_ellipseAspectRatio,   &
                              sim_k

  implicit none

  real :: AxA
  real :: sqrt2
!
!
!   ...Calculate the 'k' value.
!
!
  AxA = sim_ellipseAspectRatio * sim_ellipseAspectRatio

  sim_k = (AxA + 1.0) / (AxA - 1.0)
!
!
!   ...Calculate the parameters that will characterize the elliptical path.
!
!
  sqrt2 = sqrt (2.0)

  sim_ellipseMinorSemiAxis = sqrt2
  sim_ellipseMajorSemiAxis = sqrt2 * sim_ellipseAspectRatio
  sim_ellipseRotationAngle = 45.0                             ! in degrees
  sim_ellipseCenterX       = 0.0
  sim_ellipseCenterY       = 0.0
!
!
!   ...Ready!
!
!
  return
end subroutine sim_calculateInitialData
