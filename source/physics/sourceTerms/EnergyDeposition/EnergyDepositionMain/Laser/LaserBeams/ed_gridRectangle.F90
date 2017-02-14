!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_gridRectangle
!!
!! NAME
!!
!!  ed_gridRectangle
!!
!! SYNOPSIS
!!
!!  call ed_gridRectangle (real    (in)    :: aspectRatio,
!!                         integer (inout) :: nGridPoints,
!!                         integer (out)   :: nTicsSemiaxisMajor,
!!                         integer (out)   :: nTicsSemiaxisMinor,
!!                         real    (out)   :: deltaNormalizedMajor,
!!                         real    (out)   :: deltaNormalizedMinor)
!!
!! DESCRIPTION
!!
!!  Given a rectangle defined by its aspect ratio (a/b) of its two major (a) and minor (b) semiaxis,
!!  the routine tries to place a rectangular grid onto it (covering the whole rectangle area), such
!!  that the number of grid points inside the rectangle is close to the number of such grid points
!!  wanted. Shown below is a picture in which one quadrant of the rectangle has been filled with the
!!  rectangular grid.
!!
!!
!!                            -----------------------------
!!                           |              |__|__|__|__|__|
!!                           |            b |__|__|__|__|__|
!!                           |              |  |  |  |  |  |
!!                            --------------+--------------
!!                           |              |      a       |
!!                           |              |              |
!!                           |              |              |
!!                            -----------------------------
!!
!!  The routine deals only with the normalized rectangle, which is obtained from the original
!!  rectangle by downscaling with the major semiaxis (a). The normalized rectangle has then
!!  a major semiaxis value of 1 and a minor semiaxis value of (aspect ratio)^(-1).
!!  The rectangular grid is completely defined by the number of tics along both semiaxes as
!!  well as the normalized separation between consecutive tics (in units of the major semiaxis).
!!  The routine only defines a grid in case the number of requested grid points is > 1.
!!  Nothing is done, if this number is =< 1.
!!
!!
!!  Algorithm:
!!
!!      Let N be the number of grid points wanted inside a rectanle with semiaxis (a,b;a>=b).
!!      Denote the two possible ratios of the two semiaxes by:
!!
!!                                         f = b / a
!!                                         g = a / b (aspect ratio)
!!
!!      In the rectangle, if A is the number of tics along the 'a' semiaxis starting with 0, then
!!      there will be:
!!
!!                                    grid points along the axes =  2A + 2fA + 1
!!             grid points in the 4 quadrants excluding the axes =  4f * A^2
!!
!!      Hence we have to solve the quadratic equation in A:
!!
!!                     4f * A^2  +  (2f + 2) * A  + (1 - N)  =  0
!!
!!      whose solution is:
!!
!!                      A = (1/4) * [ - (1 + g) + sqrt [(1 + g)^2 + 4(N - 1)g] ]
!!
!!      Since A has to be integer >=0, we take the ceiling value. Note that A is equal to zero
!!      only when N = 1. Next we find B, the number of tics along the 'b' semiaxis. Since the
!!      total number of tics is 4AB + 2(A+B) + 1 and this number has to be close to N, we can
!!      set:
!!
!!                      B = (N - 2A - 1) / (4A + 2)
!!
!!      and take again the ceiling of it. We need to make sure however that B >= 1. This has
!!      given us the initial A and B values. Since taking the ceiling for B and enforcing B >= 1,
!!      we know that the number of grid points will be >= N, i.e.:
!!
!!                    # grid points = 4AB + 2(A+B) + 1 >= N  
!!
!!      To get better (# grid points closer to N) values for A and B, we test the following
!!      three combinations (this is why we needed to enforce B >= 1):
!!
!!                         A and B , A and B - 1 , A + 1 and B - 1
!!
!!      and choose the combination that leads to the # of grid points closest to N.
!!      Once the final A and B have been chosen, the normalized tic spacing corresponding to
!!      the major and minor semiaxes are:
!!
!!               tic spacing (delta) for major semiaxis = 1 / A     (set to 0, if A = 0)
!!               tic spacing (delta) for minor semiaxis = 1 / (Bg)  (set to 0, if B = 0)
!!
!!
!! ARGUMENTS
!!
!!  aspectRatio          : the aspect ratio (major semiaxis divided by minor semiaxis) of the rectangle
!!  nGridPoints          : on input  -> the number of grid points wanted
!!                         on output -> the optimum number of grid points close to the input value
!!  nTicsSemiaxisMajor   : the number of tics on the major rectangle semiaxis for the rectangular grid 
!!  nTicsSemiaxisMinor   : the number of tics on the minor rectangle semiaxis for the rectangular grid 
!!  deltaNormalizedMajor : the normalized tic spacing on the major rectangle semiaxis
!!  deltaNormalizedMinor : the normalized tic spacing on the minor rectangle semiaxis
!!
!! NOTES
!!
!!  Only the aspect ratio is needed for this routine. Hence the absolute size of the rectangle
!!  is not important. The returned tic separation value is in units of the major semiaxis, thus
!!  the grid is valid for the whole class of rectangles belonging to the specified aspect ratio.
!!
!!***

subroutine ed_gridRectangle (aspectRatio,                       &
                                          nGridPoints,          &
                                          nTicsSemiaxisMajor,   &
                                          nTicsSemiaxisMinor,   &
                                          deltaNormalizedMajor, &
                                          deltaNormalizedMinor  )

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  real,    intent (in)    :: aspectRatio
  integer, intent (inout) :: nGridPoints
  integer, intent (out)   :: nTicsSemiaxisMajor
  integer, intent (out)   :: nTicsSemiaxisMinor
  real,    intent (out)   :: deltaNormalizedMajor
  real,    intent (out)   :: deltaNormalizedMinor

  integer :: A,B,N
  integer :: diff1, diff2, diff3, diffmin
  integer :: ngp1, ngp2, ngp3

  real    :: r,t
!
!
!     ...Return, if the number of grid points wanted is =< 1. Catch a bad aspect ratio.
!
!
  if (nGridPoints <= 1) then
      return
  end if

  if (aspectRatio < 1.0) then
      call Driver_abortFlash ("ed_gridRectangle: Bad rectangle aspect ratio!")
  end if
!
!
!     ...Calculate initial number of grid points A and B
!
!
  N = nGridPoints

  r = 1.0 + aspectRatio
  t = 4.0 * real (N - 1) * aspectRatio

  A = ceiling ( 0.25 * ( - r + sqrt (r * r + t)) )               ! A is >= 1, because N > 1 -> t > 0
  B = ceiling ( real (N - 2 * A - 1) / real (4 * A + 2) )
  B = max (B,1)                                                  ! B is >= 1
!
!
!     ...Readjust A and B to optimum values and set the number
!        of tics for both major and minor semiaxes. 
!
!
  ngp1 = 4 * A * B + 2 * (A + B) + 1                             ! # of grid points for A and B
  ngp2 = ngp1 - (4 * A + 2)                                      ! # of grid points for A and B - 1
  ngp3 = ngp1 - 4 * (A - B + 1)                                  ! # of grid points for A + 1 and B - 1

  diff1 = abs (ngp1 - N)
  diff2 = abs (ngp2 - N)
  diff3 = abs (ngp3 - N)

  diffmin = min (diff1,diff2,diff3)

  if (diffmin == diff2) then
      B = B - 1
  else if (diffmin == diff3) then
      A = A + 1
      B = B - 1
  end if

  nTicsSemiaxisMajor = A
  nTicsSemiaxisMinor = B
!
!
!     ...Calculate the nomalized tic spacings. 
!
!
  deltaNormalizedMajor = 1.0 / real (A)
  deltaNormalizedMinor = 1.0 / (real (B) * aspectRatio)
!
!
!     ...Ready!
!
!
  return
end subroutine ed_gridRectangle
