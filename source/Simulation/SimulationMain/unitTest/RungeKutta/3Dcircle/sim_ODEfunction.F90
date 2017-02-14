!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/3Dcircle/sim_ODEfunction
!!
!! NAME
!!
!!  sim_ODEfunction
!! 
!! SYNOPSIS
!!
!!  sim_ODEfunction (real, intent (in) :: t,
!!                   real, intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Defines the ODE function to be used for the Runge Kutta unit test. The function
!!  corresponds to a classical path evaluation function (input x,y,z positions and
!!  x,y,z velocities) given a well defined acceleration pattern.
!!
!! ARGUMENTS 
!!
!! NOTES
!!
!!  Even though this particular problem does not depend on time explicitly, a time
!!  variable must be defined as possible input argument in order to be compatible
!!  with the Runge Kutta unit routines. 
!!
!!***

function sim_ODEfunction (t,y)

  use  Simulation_data, ONLY: sim_speed

  implicit none

  real, intent (in) :: t                    ! it is absolutely mandatory to
  real, intent (in) :: y (:)                ! declare the variables and
                                            ! array function in the way shown.
  real :: sim_ODEfunction (1:size (y))      ! (compatible to Runge Kutta interfaces)

  real :: a, r, v
  real :: ax, ay, az
  real :: factor
  real :: rx, ry, rz
  real :: vx, vy, vz
!
!
!   ...Define the ODE function.
!
!  
  v = sim_speed

  rx = y (1)
  ry = y (2)
  rz = y (3)
  vx = y (4)
  vy = y (5)
  vz = y (6)

  r = sqrt (rx * rx + ry * ry + rz * rz)
  a = v * v / r                                ! acceleration magnitude

  factor = a / r

  ax = - rx * factor                           ! acceleration x component
  ay = - ry * factor                           ! acceleration y component
  az = - rz * factor                           ! acceleration z component

  sim_ODEfunction (1) = vx
  sim_ODEfunction (2) = vy
  sim_ODEfunction (3) = vz
  sim_ODEfunction (4) = ax
  sim_ODEfunction (5) = ay
  sim_ODEfunction (6) = az
!
!
!   ...Ready!
!
!
  return
end function sim_ODEfunction
