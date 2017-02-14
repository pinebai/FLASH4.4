!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonUtilities/pi_statisticalXYZsphere
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

  use ProtonImaging_data, ONLY : pi_largestPositiveInteger,    &
                                 pi_randomNumberSeedArray,     &
                                 pi_randomNumberSeedIncrement, &
                                 pi_randomNumberSeedInitial

  use Driver_interface,   ONLY : Driver_abortFlash

  use ut_randomInterface, ONLY : ut_randomNumberArray, &
                                 ut_randomSeed

  implicit none

  real,    intent (in)    :: radius
  integer, intent (in)    :: arraySize
  integer, intent (out)   :: nSphere
  real,    intent (out)   :: xSphere (1:arraySize)
  real,    intent (out)   :: ySphere (1:arraySize)
  real,    intent (out)   :: zSphere (1:arraySize)

  logical :: inSphere
  integer :: maxSeed
  integer :: n
  real    :: x,y,z
!
!
!     ...Change and commit the seed first (not needed at the moment).
!
!
!  maxSeed = maxval (pi_randomNumberSeedArray)
!  if (maxSeed > pi_largestPositiveInteger - pi_randomNumberSeedIncrement) then
!      pi_randomNumberSeedArray (:) = pi_randomNumberSeedInitial
!  else
!      pi_randomNumberSeedArray (:) = pi_randomNumberSeedArray (:) + pi_randomNumberSeedIncrement
!  end if
!  call ut_randomSeed (ut_put = pi_randomNumberSeedArray)
!
!
!     ...Retrieve the current array of (x,y,z) sphere triples.
!
!
  call ut_randomNumberArray (xSphere)     ! fill up x-coordinate array with random numbers
  call ut_randomNumberArray (ySphere)     ! fill up y-coordinate array with random numbers
  call ut_randomNumberArray (zSphere)     ! fill up z-coordinate array with random numbers

  nSphere = 0

  do n = 1, arraySize

     x = xSphere (n)                             ! random number between [0,1]
     y = ySphere (n)                             ! random number between [0,1]
     z = zSphere (n)                             ! random number between [0,1]

     x = x + x - 1.0                             ! shift to random number between [-1,1]
     y = y + y - 1.0                             ! shift to random number between [-1,1]
     z = z + z - 1.0                             ! shift to random number between [-1,1]

     inSphere = (x * x + y * y + z * z <= 1.0)   ! check, if inside unit sphere

     if (inSphere) then
         nSphere = nSphere + 1
         xSphere (nSphere) = x
         ySphere (nSphere) = y
         zSphere (nSphere) = z
     end if

  end do

  xSphere (1:nSphere) = radius * xSphere (1:nSphere)   ! unit sphere -> real sphere
  ySphere (1:nSphere) = radius * ySphere (1:nSphere)   ! unit sphere -> real sphere
  zSphere (1:nSphere) = radius * zSphere (1:nSphere)   ! unit sphere -> real sphere
!
!
!     ...Ready!
!
!
  return
end subroutine pi_statisticalXYZsphere
