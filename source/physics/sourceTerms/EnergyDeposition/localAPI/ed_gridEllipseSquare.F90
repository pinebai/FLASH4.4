!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_gridEllipseSquare
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

  implicit none

  real,    intent (in)    :: aspectRatio
  integer, intent (inout) :: nGridPoints
  integer, intent (out)   :: nTicsSemiaxisMajor
  integer, intent (out)   :: nTicsSemiaxisMinor
  real,    intent (out)   :: deltaNormalized

  nGridPoints        = 0
  nTicsSemiaxisMajor = 0
  nTicsSemiaxisMinor = 0
  deltaNormalized    = 0.0

  return
end subroutine ed_gridEllipseSquare
