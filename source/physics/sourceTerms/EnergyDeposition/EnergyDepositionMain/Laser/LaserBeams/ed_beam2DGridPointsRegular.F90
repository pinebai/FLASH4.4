!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beam2DGridPointsRegular
!!
!! NAME
!!
!!  ed_beam2DGridPointsRegular
!!
!! SYNOPSIS
!!
!!  call ed_beam2DGridPointsRegular (real,    intent (in)    :: semiAxis,
!!                                   integer, intent (in)    :: nTics,
!!                                   real,    intent (in)    :: delta,
!!                                   real,    intent (in)    :: firstTic,
!!                                   logical, intent (inout) :: startGrid,
!!                                   integer, intent (in)    :: maxGridPoints,
!!                                   logical, intent (out)   :: moreGridPoints,
!!                                   integer, intent (out)   :: nGridPoints,
!!                                   real,    intent (out)   :: xGrid (1:maxGridPoints))
!!
!! DESCRIPTION
!!
!!  Returns a collection of regular grid points for a 2D beam. The regular grid has to be
!!  set up beforehand to be able to use this routine. The total number of regular grid
!!  points of the regular grid is cut into several equally sized arrays, each of which
!!  is returned consecutively by a call to this routine.
!!
!! ARGUMENTS
!!
!!  semiAxis       : the semiaxis (1/2 length) of the 2D beam cross sectional area
!!  nTics          : the total number of tics on the grid (equal to the number of grid points)
!!  delta          : the tic spacing (for regular grids)
!!  firstTic       : the position of the 1st tic from the lowest boundary
!!  startGrid      : if true, the grid points will start from the beginning 
!!  maxGridPoints  : the maximum number of grid points that can be returned
!!  moreGridPoints : if true, more grid points are expected
!!  nGridPoints    : the actual number of grid points returned
!!  xGrid          : the semiaxis based coordinates of the grid points
!!
!! NOTES
!!
!!***

subroutine ed_beam2DGridPointsRegular (semiAxis,                      &
                                       nTics,                         &
                                       delta,                         &
                                       firstTic,                      &
                                       startGrid,                     &
                                       maxGridPoints,                 &
                                                      moreGridPoints, &
                                                      nGridPoints,    &
                                                      xGrid           )

  use Driver_interface,   ONLY : Driver_abortFlash

  implicit none

  real,    intent (in)    :: semiAxis
  integer, intent (in)    :: nTics
  real,    intent (in)    :: delta
  real,    intent (in)    :: firstTic
  logical, intent (inout) :: startGrid
  integer, intent (in)    :: maxGridPoints
  logical, intent (out)   :: moreGridPoints
  integer, intent (out)   :: nGridPoints
  real,    intent (out)   :: xGrid (1:maxGridPoints)

  integer       :: i
  integer, save :: istart           ! save semiaxis grid starting point for next call

  real          :: x
!
!
!     ...If requested, start the grid points.
!
!
  if (startGrid) then
      istart    = 1
      startGrid = .false.
  end if
!
!
!     ...Retrieve the current array of grid points.
!
!
  nGridPoints = 0

  do i = istart, nTics
     x = - semiAxis + firstTic + real (i - 1) * delta

     nGridPoints = nGridPoints + 1

     if (nGridPoints > maxGridPoints) then
         istart = i                     ! save current semiaxis grid starting point for next call
         nGridPoints = nGridPoints - 1  ! adjust to the true # of grid points returned
         moreGridPoints = .true.        ! this is always the case
         return
     else
         xGrid (nGridPoints) = x
         moreGridPoints = .not. (i == nTics)
     end if
  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beam2DGridPointsRegular
