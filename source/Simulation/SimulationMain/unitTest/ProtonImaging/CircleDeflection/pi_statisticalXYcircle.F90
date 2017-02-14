!!****if* source/Simulation/SimulationMain/unitTest/ProtonImaging/CircleDeflection/pi_statisticalXYcircle
!!
!! NAME
!!
!!  pi_statisticalXYcircle
!!
!! SYNOPSIS
!!
!!  call pi_statisticalXYcircle (real,    intent (in)  :: radius,
!!                               integer, intent (in)  :: arraySize,
!!                               integer, intent (out) :: nXY,
!!                               real,    intent (out) :: xCircle (1:arraySize),
!!                               real,    intent (out) :: yCircle (1:arraySize))
!!
!! DESCRIPTION
!!
!!  Overriding version of the original. This particular version returns chuncks of (x,y)
!!  pairs sitting on the circumference of a circle of certain radius. There is no
!!  statistical distribution of the (x,y) pairs on the circle. The returned (x,y) pairs
!!  lay sequentially on the circumference of the circle.
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

  use Simulation_data,    ONLY : sim_nBeamPoints, &
                                 sim_beamPointsX, &
                                 sim_beamPointsY

  implicit none

  real,    intent (in)  :: radius
  integer, intent (in)  :: arraySize
  integer, intent (out) :: nCircle
  real,    intent (out) :: xCircle (1:arraySize)
  real,    intent (out) :: yCircle (1:arraySize)

  integer, save :: collectedXYpairs = 0    ! save for next call, initial value is zero

  integer :: n
!
!
!     ...Retrieve the current array of circumference (x,y) circle pairs.
!
!
  nCircle = 0

  do n = 1, arraySize

     collectedXYpairs = collectedXYpairs + 1

     if (collectedXYpairs <= sim_nBeamPoints) then
         nCircle = nCircle + 1
         xCircle (nCircle) = radius * sim_beamPointsX (collectedXYpairs)
         yCircle (nCircle) = radius * sim_beamPointsY (collectedXYpairs)
     else
         return
     end if

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine pi_statisticalXYcircle
