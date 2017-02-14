!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/2Dellipse/sim_ODEfunction
!!
!! NAME
!!
!!  sim_ODEfunction
!! 
!! SYNOPSIS
!!
!!  sim_ODEfunction (real, intent (in) :: t,
!!                   real, intent (in) :: p (:))
!!
!! DESCRIPTION
!!
!!  Defines the ODE function to be used for the 2D ellipse Runge Kutta unit test. The function
!!  is defined through the 2D ODE:
!!
!!                   dp     |  1   k |   | x |
!!                   --  =  |        | * |   |
!!                   dt     | -k  -1 |   | y |
!!
!!  or, writing each component explicitly:
!!
!!                   x  =   x + k * y
!!                   y  = - y - k * x
!!
!!  where p is the 2D position vector and x/y are its cartesian 2D components. The constant k
!!  defines the ellipticity of the problem, with k > 1. For k close to 1, the ellipse has a high
!!  aspect ratio, almost a line with long extension in space. For very large k, the path taken
!!  by the particle becomes almost a circle.
!!
!! ARGUMENTS 
!!
!!  t : time variable
!!  p : the x/y position components: p (1) -> x and p (2) -> y
!!
!! NOTES
!!
!!  Even though this particular problem does not depend on time explicitly, a time
!!  variable must be defined as possible input argument in order to be compatible
!!  with the Runge Kutta unit routines. 
!!
!!***

function sim_ODEfunction (t,p)

  use  Simulation_data, ONLY: sim_k

  implicit none

  real, intent (in) :: t                    ! it is absolutely mandatory to
  real, intent (in) :: p (:)                ! declare the variables and
                                            ! array function in the way shown.
  real :: sim_ODEfunction (1:size (p))      ! (compatible to Runge Kutta interfaces)
!
!
!   ...Define the ODE function.
!
!  
  sim_ODEfunction (1) =   p (1) + sim_k * p (2)
  sim_ODEfunction (2) = - p (2) - sim_k * p (1)
!
!
!   ...Ready!
!
!
  return
end function sim_ODEfunction
