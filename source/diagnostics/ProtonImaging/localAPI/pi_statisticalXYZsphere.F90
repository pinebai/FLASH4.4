!!****if* source/diagnostics/ProtonImaging/localAPI/pi_statisticalXYZsphere
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
!!  Returns a collection of statistical (x,y,z) triples within a sphere of certain radius.
!!  The random number generator has to be inititalized (seeded) in order to be able to use
!!  this routine. The supplied (x,y,z) arrays will be used first to collect a sequence of
!!  random number triples. The (x,y,z) triples in the circle are placed at the beginning of
!!  the (x,y,z) arrays.
!!
!! ARGUMENTS
!!
!!  radius    : radius of the sphere
!!  arraySize : maximum number of statistical (x,y,z) triples that can be returned per call
!!  nSphere   : the actual number of statistical (x,y,z) sphere triples returned
!!  xSphere   : the x-coordinates of the statistical (x,y,z) sphere triples
!!  ySphere   : the y-coordinates of the statistical (x,y,z) sphere triples
!!  zSphere   : the z-coordinates of the statistical (x,y,z) sphere triples
!!
!! NOTES
!!
!!  The current algorithm uses the sphere (x,y,z) arrays for harvesting the random numbers needed.
!!  Since some of these numbers lay not within the sphere, the number of accepted random numbers
!!  is approximately pi/6 or 1/2 less than the maximum number of random numbers harvested.
!!  In order for this procedure to be efficient in collecting all the sphere triples, the size of
!!  the (x,y,z) arrays should not be too small.
!!
!!  Needs the (seeded) random number generator as provided in the FLASH untilities directory.
!!
!!***

subroutine pi_statisticalXYZsphere (radius,             &
                                    arraySize,          &
                                               nSphere, &
                                               xSphere, &
                                               ySphere, &
                                               zSphere  )
  implicit none

  real,    intent (in)    :: radius
  integer, intent (in)    :: arraySize
  integer, intent (out)   :: nSphere
  real,    intent (out)   :: xSphere (1:arraySize)
  real,    intent (out)   :: ySphere (1:arraySize)
  real,    intent (out)   :: zSphere (1:arraySize)

  nSphere               = 0
  xSphere (1:arraySize) = 0.0
  ySphere (1:arraySize) = 0.0
  zSphere (1:arraySize) = 0.0

  return
end subroutine pi_statisticalXYZsphere
