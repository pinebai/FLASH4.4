!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/3Dcircle/sim_calculateInitialData
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
!!  Calculates the radius and the initial velocity components according to the starting
!!  position. The velocity vector must be at right angles to the position vector (scalar
!!  product = 0) and the velocity vector must be contained in the plane formed by the
!!  position vector and the z axis.
!!
!!  The radius is:  sqrt (rx^2 + ry^2 + rz^2)
!!
!!  The velocity components are:
!!
!!                                rx * rz * |v|
!!              vx  =  -/+  ------------------------
!!                          |r| * sqrt (rx^2 + ry^2)
!!
!!
!!                                ry * rz * |v|
!!              vy  =  -/+  ------------------------
!!                          |r| * sqrt (rx^2 + ry^2)
!!
!!
!!                          |v| * sqrt (rx^2 + ry^2)
!!              vz  =  +/-  ------------------------
!!                                    |r|
!!
!!
!!  where |r| and |v| are the values of the radius and the speed, respectively. The
!!  two possible sign choices come from the two possible orientations of the velocity
!!  vector on the plane. We take the (+) choice for vx/vy and the (-) choice for vz.
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_calculateInitialData ()

  use  Simulation_data, ONLY: sim_radius, &
                              sim_rx0,    &
                              sim_ry0,    &
                              sim_rz0,    &
                              sim_speed,  &
                              sim_vx0,    &
                              sim_vy0,    &
                              sim_vz0

  implicit none

  real :: r, v
  real :: rx, ry, rz
  real :: rxrx, ryry, rzrz
  real :: rxy
  real :: vx, vy, vz
!
!
!   ...Calculate radius.
!
!
  v  = sim_speed
  rx = sim_rx0
  ry = sim_ry0
  rz = sim_rz0

  rxrx = rx * rx
  ryry = ry * ry
  rzrz = rz * rz

  r = sqrt (rxrx + ryry + rzrz)
  rxy = sqrt (rxrx + ryry)
!
!
!   ...Calculate velocity components.
!
!
  vx = (rx * rz * v) / (r * rxy)
  vy = (ry * rz * v) / (r * rxy)
  vz =   - (rxy * v) /  r
!
!
!   ...Store data.
!
!
  sim_radius = r
  sim_vx0    = vx
  sim_vy0    = vy
  sim_vz0    = vz
!
!
!   ...Ready!
!
!
  return
end subroutine sim_calculateInitialData
