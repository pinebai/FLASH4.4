!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_gridEllipseSquare
!!
!! NAME
!!
!!  ed_gridEllipseSquare
!!
!! SYNOPSIS
!!
!!  call ed_gridEllipseSquare (real    (in)    :: aspectRatio,
!!                             integer (inout) :: nGridPoints,
!!                             integer (out)   :: nTicsSemiaxisMajor,
!!                             integer (out)   :: nTicsSemiaxisMinor,
!!                             real    (out)   :: deltaNormalized)
!!
!! DESCRIPTION
!!
!!  Given an ellipse defined by its aspect ratio (a/b) of its two major (a) and minor (b) semiaxis,
!!  the routine tries to place a square grid onto it, such that the number of grid points inside
!!  the ellipse is close to the number of such grid points wanted. Shown below is a picture in which
!!  one quadrant of the ellipse has been filled with the square grid.
!!
!!
!!                                  *  |__*_
!!                              *      |__|__ *
!!                           *       b |__|__|__ *
!!                         *           |__|__|__|__*
!!                        *            |  |  |  |  |*
!!                        *------------+------------*
!!                        *            |     a      *
!!                         *           |           *
!!                           *         |         *
!!                              *      |      *
!!                                  *  |  *
!!
!!
!!  The routine deals only with the normalized ellipse, which is obtained from the original
!!  ellipse by downscaling with the major semiaxis (a). The normalized ellipse has then
!!  a major semiaxis value of 1 and a minor semiaxis value of (aspect ratio)^(-1).
!!  The square grid is completely defined by the number of tics along both semiaxes as
!!  well as the normalized separation between consecutive tics (in units of the major semiaxis).
!!  The routine only defines a grid in case the number of requested grid points is > 1.
!!  Nothing is done, if this number is =< 1.
!!
!!
!!  Algorithm:
!!
!!      If N is the number of grid points wanted inside an ellipse with semiaxis (a,b;a>=b),
!!      then the number of grid points inside the circumscribing rectangle is approximately:
!!
!!           # of grid points in circumscribing rectangle = N * (area rectangle) / (area ellipse)
!!
!!                                                        = N * ([4ab]/[a*b*pi])
!!
!!                                                        = N * (4/pi)
!!
!!      Denote the two possible ratios of the two semiaxes by:
!!
!!                                         f = b / a
!!                                         g = a / b
!!
!!      In the rectangle, if A is the number of tics along the 'a' semiaxis starting with 0, then
!!      there will be:
!!
!!                                    grid points along the axes =  2A + 2fA + 1
!!             grid points in the 4 quadrants excluding the axes =  4f * A^2
!!
!!      Hence we have to solve the quadratic equation in A:
!!
!!                     4f * A^2  +  (2f + 2) * A  -  (4N/pi - 1)  =  0
!!
!!      whose solution is:
!!
!!                      A = (1/4) * [ - (1 + g) + sqrt [(1 - g)^2 + 16Ng/pi] ]
!!
!!      Since A has to be integer, after some tests we decided to take the ceiling value,
!!      since this usually gives a closer # of grid points when compared to N. After A has been 
!!      determined, there is one more refinement step concerning the delta value (separation
!!      between two consecutive tics). The delta values are examined in between the range:
!!
!!                      # of tics = A - 1  ---->  # of tics = A + 1
!!
!!      over a predefined number of steps and the average delta value leading to the number of
!!      grid points closest to N is chosen. Note, that it is not guaranteed that this way we obtain
!!      the truly optimum delta value leading to the optimum number of grid points, but in the
!!      vast majority of cases we are very close.
!!
!!
!! ARGUMENTS
!!
!!  aspectRatio        : the aspect ratio (major semiaxis divided by minor semiaxis) of the ellipse
!!  nGridPoints        : on input  -> the number of grid points wanted
!!                       on output -> the optimum number of grid points close to the input value
!!  nTicsSemiaxisMajor : the number of tics on the major semiaxis for the square grid 
!!  nTicsSemiaxisMinor : the number of tics on the minor semiaxis for the square grid 
!!  deltaNormalized    : the normalized separation distance between two consecutive tics
!!                       (in units of major semiaxis)
!!
!! NOTES
!!
!!  Only the aspect ratio is needed for this routine. Hence the absolute size of the ellipse
!!  is not important. The returned tic separation value is in units of the major semiaxis, thus
!!  the grid is valid for the whole class of ellipses belonging to the specified aspect ratio.
!!  The code can be thought off as getting the square grid information from a major semiaxis
!!  normalized (major semiaxis = 1) ellipse.
!!
!!***

subroutine ed_gridEllipseSquare (aspectRatio,                     &
                                              nGridPoints,        &
                                              nTicsSemiaxisMajor, &
                                              nTicsSemiaxisMinor, &
                                              deltaNormalized     )

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"

  real,    intent (in)    :: aspectRatio
  integer, intent (inout) :: nGridPoints
  integer, intent (out)   :: nTicsSemiaxisMajor
  integer, intent (out)   :: nTicsSemiaxisMinor
  real,    intent (out)   :: deltaNormalized

  integer :: aTics, bTics
  integer :: deltaCount
  integer :: minimum
  integer :: n
  integer :: nAxisPoints
  integer :: nEllipsePoints
  integer :: nOptimumPoints
  integer :: nPointsDifference
  integer :: nQuadrantPoints
  integer :: nSteps
  integer :: step
  integer :: x

  real    :: aspectRatioInv
  real    :: ax, by
  real    :: delta
  real    :: deltaInc, deltaInv, deltaSum
  real    :: deltaMin, deltaMax
  real    :: r,s,t

  integer, parameter :: nStepsDefault = 1000       ! controls optimicity quality of final delta result
!
!
!     ...Return, if the number of grid points wanted is =< 1. Catch a bad aspect ratio.
!
!
  if (nGridPoints <= 1) then
      return
  end if

  if (aspectRatio < 1.0) then
      call Driver_abortFlash ("ed_gridEllipseSquare: Bad ellipse aspect ratio!")
  end if
!
!
!     ...Calculate: 1) initial number of grid points
!                   2) delta (normalized) range to be tested (maximum and minimum)
!                   3) delta (normalized) increments
!
!        Note, that from the construction formula for the initial number 'n' of grid points
!        we should always get n > 0. This is easily seen by the lowest possible number of
!        grid points (= 1) and the lowest possible aspect ratio (= 1). In this case we get for
!        the variables below: aspectRatio = 1.0, r = 2.0, s = 0.0, t = 5.0929... and therefore
!        n = ceiling (0.06418...) = 1. If n = 1, then the search for an optimum delta value
!        is bypassed.
!
!
  r = 1.0 + aspectRatio
  s = 1.0 - aspectRatio
  t = 16.0 * real (nGridPoints) * aspectRatio / PI

  n = ceiling ( 0.25 * ( - r + sqrt (s * s + t)) )

  if (n == 1) then
      nSteps   = 0            ! means no refinement
      deltaMax = 1.0          ! normalized, in units of major semiaxis
      deltaInc = 0.0
  else
      nSteps   = nStepsDefault
      deltaMax = 1.0 / (real (n - 1))
      deltaMin = 1.0 / (real (n + 1))
      deltaInc = (deltaMax - deltaMin) / real (nSteps)
  end if
!
!
!     ...Perform the one-stage refinement for the delta value. 
!
!
  delta = deltaMax
  minimum = huge (1)
  aspectRatioInv = 1.0 / aspectRatio

  do step = 0,nSteps

     deltaInv = 1.0 / delta

     aTics = int (deltaInv)
     bTics = int (deltaInv * aspectRatioInv)

     nAxisPoints     = aTics + aTics + bTics + bTics + 1     ! includes the (0,0) point
     nQuadrantPoints = 0

     do x = 1, aTics
        ax = real (x) * delta                                ! real x-coordinate of normalized ellipse
        by = aspectRatioInv * sqrt (1.0 - ax * ax)           ! real y-coordinate of normalized ellipse
        bTics = int (by * deltaInv)
        nQuadrantPoints = nQuadrantPoints + bTics
     end do

     nEllipsePoints    = 4 * nQuadrantPoints + nAxisPoints
     nPointsDifference = abs (nGridPoints - nEllipsePoints)

     if (nPointsDifference < minimum) then

         minimum = nPointsDifference
         deltaSum = delta
         deltaCount = 1
         nOptimumPoints = nEllipsePoints

     else if (nPointsDifference == minimum) then

         deltaSum = deltaSum + delta
         deltaCount = deltaCount + 1

     else
         exit
     end if

     delta = delta - deltaInc

  end do
!
!
!     ...Calculate final normalized delta value and corresponding number of tics on both semiaxes.
!
!
  delta = deltaSum / real (deltaCount)
  deltaInv = 1.0 / delta

  nGridPoints        = nOptimumPoints
  nTicsSemiaxisMajor = int (deltaInv)
  nTicsSemiaxisMinor = int (deltaInv * aspectRatioInv)
  deltaNormalized    = delta
!
!
!     ...Ready!
!
!
  return
end subroutine ed_gridEllipseSquare
