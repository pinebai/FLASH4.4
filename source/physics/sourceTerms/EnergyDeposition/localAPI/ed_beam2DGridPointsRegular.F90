!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beam2DGridPointsRegular
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

  startGrid               = .false.
  moreGridPoints          = .false.
  nGridPoints             =  0
  xGrid (1:maxGridPoints) =  0.0

  return
end subroutine ed_beam2DGridPointsRegular
