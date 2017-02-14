!!****if* source/diagnostics/ProtonImaging/localAPI/pi_statisticalXYcircle
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

  implicit none

  real,    intent (in)  :: radius
  integer, intent (in)  :: arraySize
  integer, intent (out) :: nCircle
  real,    intent (out) :: xCircle (1:arraySize)
  real,    intent (out) :: yCircle (1:arraySize)

  nCircle               = 0
  xCircle (1:arraySize) = 0.0
  yCircle (1:arraySize) = 0.0

  return
end subroutine pi_statisticalXYcircle
