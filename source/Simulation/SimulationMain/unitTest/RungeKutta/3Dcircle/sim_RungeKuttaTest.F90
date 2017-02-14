!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/3Dcircle/sim_RungeKuttaTest
!!
!! NAME
!!
!!  sim_RungeKuttaTest
!! 
!! SYNOPSIS
!!
!!  call sim_RungeKuttaTest
!!
!! DESCRIPTION
!!
!!  This function tests the Runge Kutta stepper of the Runge Kutta utility unit.
!!  A ODE function is defined for this test and the stepper is called many times
!!  recording advancement of the function + error.
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_RungeKuttaTest ()

  use  Simulation_data,      ONLY: sim_stepSize,                &
                                   sim_errorFraction,           &
                                   sim_numberOfCircles,         &
                                   sim_numberOfRungeKuttaSteps, &
                                   sim_RungeKuttaMethod,        &
                                   sim_rx0,                     &
                                   sim_ry0,                     &
                                   sim_rz0,                     &
                                   sim_vx0,                     &
                                   sim_vy0,                     &
                                   sim_vz0

  use  RungeKutta_interface, ONLY: RungeKutta_Step,             &
                                   RungeKutta_StepSizeEstimate

  use  sim_interface,        ONLY: sim_ODEfunction

  implicit none

  integer :: i,n
  integer :: nHalfCircle


  integer :: nstep
  real    :: hsum


  real    :: ep, ev
  real    :: htry, hused, hnext
  real    :: t
  real    :: r,v
  real    :: rx, ry, rz
  real    :: rzprev
  real    :: vx, vy, vz

  real    :: y     (1:6)
  real    :: ynew  (1:6)
  real    :: error (1:6)
  real    :: eBase (1:6)
!
!
!   ...Set error base for all dependent variables.
!
!  
!  eBase (1) = 100.0
!  eBase (2) = 100.0
!  eBase (3) = 100.0
!  eBase (4) = 100.0
!  eBase (5) = 100.0
!  eBase (6) = 100.0

  eBase (1) = 1.0
  eBase (2) = 1.0
  eBase (3) = 1.0
  eBase (4) = 1.0
  eBase (5) = 1.0
  eBase (6) = 1.0
!
!
!   ...Initialize dependent variables.
!
!  
  y (1) = sim_rx0
  y (2) = sim_ry0
  y (3) = sim_rz0
  y (4) = sim_vx0
  y (5) = sim_vy0
  y (6) = sim_vz0

  rx = y (1)
  ry = y (2)
  rz = y (3)
  vx = y (4)
  vy = y (5)
  vz = y (6)

  r = sqrt (rx * rx + ry * ry + rz * rz)
  v = sqrt (vx * vx + vy * vy + vz * vz)

  ep = 0.0
  ev = 0.0

  write (*,'(a)'    ) '     radius         speed          z      error r     error v     h used   Half-Circle'
  write (*,'(a)'    ) ' -------------------------------------------------------------------------------------'
  write (*,'(f14.10,f14.10,2x,f9.4,2x,es10.3,2x,es10.3)') r,v,rz,ep,ev
!
!
!   ...Call the Runge Kutta stepper.
!
!  
  t = 0.0
!  htry = sim_stepSize

  htry = RungeKutta_StepSizeEstimate (sim_RungeKuttaMethod, &
                                      sim_ODEfunction,      &
                                      t,                    &
                                      y,                    &
                                      sim_errorFraction,    &
                                      eBase,                &
                                      hmax = 1.0            )

  write (*,*) ' h estimate = ',htry

  nstep = 0
  nHalfCircle = 0

  do while (nHalfCircle < 2 * sim_numberOfCircles)

     call RungeKutta_Step (sim_RungeKuttaMethod,   &
                           sim_ODEfunction,        &
                           t,                      &
                           y,                      &
                           sim_errorFraction,      &
                           eBase,                  &
                                            htry,  &
                                            hused, &
                                            hnext, &
                                            ynew,  &
                                            error  )


     t = t + hused
     nstep = nstep + 1

     y (:) = ynew (:)

     rzprev = rz

     rx = y (1)
     ry = y (2)
     rz = y (3)
     vx = y (4)
     vy = y (5)
     vz = y (6)

     r = sqrt (rx * rx + ry * ry + rz * rz)
     v = sqrt (vx * vx + vy * vy + vz * vz)

     ep = maxval (abs (error (1:3)))
     ev = maxval (abs (error (4:6)))

     if (abs (rzprev - rz) > max (abs (rzprev), abs (rz))) then
         nHalfCircle = nHalfCircle + 1
     end if

     write (*,'(f14.10,f14.10,2x,f9.4,2x,es10.3,2x,es10.3,2x,es10.3,3x,i3)') r,v,rz,ep,ev,hused,nHalfCircle

     htry = hnext

  end do

  write (*,*) ' h avg = ',t / real (nstep)

!
!
!   ...Ready!
!
!
  return
end subroutine sim_RungeKuttaTest
