!!****if* source/Simulation/SimulationMain/unitTest/ProtonImaging/Validation1/pi_statisticalXYcircle
!!
!! NAME
!!
!!  pi_statisticalXYcircle
!!
!! SYNOPSIS
!!
!!  call pi_statisticalXYcircle (real,    intent (in)  :: radius,
!!                               integer, intent (in)  :: arraySize,
!!                               integer, intent (out) :: nCircle,
!!                               real,    intent (out) :: xCircle (1:arraySize),
!!                               real,    intent (out) :: yCircle (1:arraySize))
!!
!! DESCRIPTION
!!
!!  Overriding version of the original. This particular version sets the (x,y) pairs
!!  equal to the x,y-coordinates of the sphere triples (x,y,z), which have been stored
!!  temporarily in 'sim_' variables.
!!
!! ARGUMENTS
!!
!!  radius    : radius of the circle
!!  arraySize : maximum number of circumference (x,y) pairs that can be returned per call
!!  nCircle   : the actual number of circumference (x,y) circle pairs returned
!!  xCircle   : the x-coordinates of the circumference (x,y) circle pairs
!!  yCircle   : the y-coordinates of the circumference (x,y) circle pairs
!!
!! NOTES
!!
!!  none
!!
!!***

subroutine pi_statisticalXYcircle (radius,             &
                                   arraySize,          &
                                              nCircle, &
                                              xCircle, &
                                              yCircle  )

  use Driver_interface,   ONLY : Driver_abortFlash

  use Simulation_data,    ONLY : sim_nCircle, &
                                 sim_xCircle, &
                                 sim_yCircle

  implicit none

  real,    intent (in)  :: radius
  integer, intent (in)  :: arraySize
  integer, intent (out) :: nCircle
  real,    intent (out) :: xCircle (1:arraySize)
  real,    intent (out) :: yCircle (1:arraySize)
!
!
!     ...Retrieve the info from the stored sphere (x,y,z) triples.
!
!
  nCircle = sim_nCircle
  xCircle (1:nCircle) = sim_xCircle (1:sim_nCircle)
  yCircle (1:nCircle) = sim_yCircle (1:sim_nCircle)
!
!
!     ...Ready!
!
!
  return
end subroutine pi_statisticalXYcircle
