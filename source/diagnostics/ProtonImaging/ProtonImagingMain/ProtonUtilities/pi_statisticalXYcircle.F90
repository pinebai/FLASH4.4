!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonUtilities/pi_statisticalXYcircle
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
!!  Returns a collection of statistical (x,y) pairs within a circle of certain radius.
!!  The random number generator has to be inititalized (seeded) in order to be able to use
!!  this routine. The supplied (x,y) arrays will be used first to collect a sequence of
!!  random number pairs. The (x,y) pairs in the circle are placed at the beginning of
!!  the (x,y) arrays.
!!
!! ARGUMENTS
!!
!!  radius    : radius of the circle
!!  arraySize : maximum number of statistical (x,y) pairs that can be returned per call
!!  nCircle   : the actual number of statistical (x,y) circle pairs returned
!!  xCircle   : the x-coordinates of the statistical (x,y) circle pairs
!!  yCircle   : the y-coordinates of the statistical (x,y) circle pairs
!!
!! NOTES
!!
!!  The current algorithm uses the circle (x,y) arrays for harvesting the random numbers needed.
!!  Since some of these numbers lay not within the circle, the number of accepted random numbers
!!  is approximately pi/4 or 3/4 less than the maximum number of random numbers harvested.
!!  In order for this procedure to be efficient in collecting all the circle pairs, the size of
!!  the (x,y) arrays should not be too small.
!!
!!  Needs the (seeded) random number generator as provided in the FLASH untilities directory.
!!
!!***

subroutine pi_statisticalXYcircle (radius,             &
                                   arraySize,          &
                                              nCircle, &
                                              xCircle, &
                                              yCircle  )

  use ProtonImaging_data, ONLY : pi_largestPositiveInteger,    &
                                 pi_randomNumberSeedArray,     &
                                 pi_randomNumberSeedIncrement, &
                                 pi_randomNumberSeedInitial

  use Driver_interface,   ONLY : Driver_abortFlash

  use ut_randomInterface, ONLY : ut_randomNumberArray, &
                                 ut_randomSeed

  implicit none

  real,    intent (in)  :: radius
  integer, intent (in)  :: arraySize
  integer, intent (out) :: nCircle
  real,    intent (out) :: xCircle (1:arraySize)
  real,    intent (out) :: yCircle (1:arraySize)

  logical :: inCircle
  integer :: maxSeed
  integer :: n
  real    :: x,y
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
!     ...Retrieve the current array of (x,y) circle pairs.
!
!
  call ut_randomNumberArray (xCircle)    ! get as many random numbers for x-coordinate as possible
  call ut_randomNumberArray (yCircle)    ! get as many random numbers for y-coordinate as possible

  nCircle = 0

  do n = 1, arraySize

     x = xCircle (n)                     ! random number between [0,1]
     y = yCircle (n)                     ! random number between [0,1]

     x = x + x - 1.0                     ! shift to random number between [-1,1]
     y = y + y - 1.0                     ! shift to random number between [-1,1]

     inCircle = (x * x + y * y <= 1.0)   ! check, if inside unit circle

     if (inCircle) then
         nCircle = nCircle + 1
         xCircle (nCircle) = x
         yCircle (nCircle) = y
     end if

  end do

  if (nCircle > 0) then
      xCircle (1:nCircle) = radius * xCircle (1:nCircle)   ! unit circle -> real circle
      yCircle (1:nCircle) = radius * yCircle (1:nCircle)   ! unit circle -> real circle
  end if
!
!
!     ...Ready!
!
!
  return
end subroutine pi_statisticalXYcircle
