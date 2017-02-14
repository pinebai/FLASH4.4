!!****if* source/Simulation/SimulationMain/unitTest/ProtonImaging/CircleDeflection/pi_statisticalXYZsphere
!!
!! NAME
!!
!!  pi_statisticalXYZsphere
!!
!! SYNOPSIS
!!
!!  call pi_statisticalXYZsphere (real,    intent (in)    :: radius,
!!                                integer, intent (in)    :: arraySize,
!!                                integer, intent (out)   :: nSphere,
!!                                real,    intent (out)   :: xSphere (1:arraySize),
!!                                real,    intent (out)   :: ySphere (1:arraySize),
!!                                real,    intent (out)   :: zSphere (1:arraySize))
!!
!! DESCRIPTION
!!
!!  Overriding version of the original. This particular version returns chuncks of (x,y,z)
!!  triples with z=0, sitting on the circumference of a sphere of certain radius. There is
!!  no statistical distribution of the (x,y,0) triples on the circumference. The returned
!!  (x,y,0) triples lay sequentially on the circumference of the sphere.
!!
!! ARGUMENTS
!!
!!  radius    : radius of the sphere
!!  arraySize : maximum number of circumference (x,y,0) triples that can be returned per call
!!  nSphere   : the actual number of circumference (x,y,0) sphere triples returned
!!  xSphere   : the x-coordinates of the circumference (x,y,0) sphere triples
!!  ySphere   : the y-coordinates of the circumference (x,y,0) sphere triples
!!  zSphere   : the z-coordinates of the circumference (x,y,0) sphere triples (all zeros)
!!
!! NOTES
!!
!!  none
!!
!!***

subroutine pi_statisticalXYZsphere (radius,             &
                                    arraySize,          &
                                               nSphere, &
                                               xSphere, &
                                               ySphere, &
                                               zSphere  )

  use Driver_interface,   ONLY : Driver_abortFlash

  use Simulation_data,    ONLY : sim_nBeamPoints, &
                                 sim_beamPointsX, &
                                 sim_beamPointsY

  implicit none

  real,    intent (in)    :: radius
  integer, intent (in)    :: arraySize
  integer, intent (out)   :: nSphere
  real,    intent (out)   :: xSphere (1:arraySize)
  real,    intent (out)   :: ySphere (1:arraySize)
  real,    intent (out)   :: zSphere (1:arraySize)

  integer, save :: collectedXYZtriples = 0    ! save for next call, initial value is zero

  integer :: n
!
!
!     ...Retrieve the current array of circumference (x,y,0) sphere triples.
!
!
  nSphere = 0

  do n = 1, arraySize

     collectedXYZtriples = collectedXYZtriples + 1

     if (collectedXYZtriples <= sim_nBeamPoints) then
         nSphere = nSphere + 1
         xSphere (nSphere) = radius * sim_beamPointsX (collectedXYZtriples)
         ySphere (nSphere) = radius * sim_beamPointsY (collectedXYZtriples)
         zSphere (nSphere) = 0.0
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
end subroutine pi_statisticalXYZsphere
