!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beam2DGridPointsStatistical
!!
!! NAME
!!
!!  ed_beam2DGridPointsStatistical
!!
!! SYNOPSIS
!!
!!  call ed_beam2DGridPointsStatistical (real,    intent (in)    :: semiAxis,
!!                                       real,    intent (in)    :: seed,
!!                                       integer, intent (in)    :: targetTotalGridPoints,
!!                                       logical, intent (inout) :: startGrid,
!!                                       integer, intent (in)    :: maxGridPoints,
!!                                       logical, intent (out)   :: moreGridPoints,
!!                                       integer, intent (out)   :: nGridPoints,
!!                                       real,    intent (out)   :: xGrid (1:maxGridPoints))
!!
!! DESCRIPTION
!!
!!  Returns a collection of statistical grid points for a 2D beam. The statistical grid has
!!  to be set up beforehand to be able to use this routine. The total number of statistical
!!  grid points of the statistical grid is cut into several equally sized arrays, each of which
!!  is returned consecutively by a call to this routine.
!!
!! ARGUMENTS
!!
!!  semiAxis              : the semiaxis of the linear 2D beam cross section
!!  seed                  : the seed value defining the grid (random number sequence)
!!  targetTotalGridPoints : the complete number of grid points wanted
!!  startGrid             : if true, the grid points will start from the beginning 
!!  maxGridPoints         : the maximum number of grid points that can be returned
!!  moreGridPoints        : if true, more grid points are expected
!!  nGridPoints           : the actual number of grid points returned
!!  xGrid                 : the semiaxis based coordinates of the grid points
!!
!! NOTES
!!
!!  All random numbers harvested by the random number generator correspond to a grid point.
!!
!!***

subroutine ed_beam2DGridPointsStatistical (semiAxis,                           &
                                           seed,                               &
                                           targetTotalGridPoints,              &
                                           startGrid,                          &
                                           maxGridPoints,                      &
                                                               moreGridPoints, &
                                                               nGridPoints,    &
                                                               xGrid           )

  implicit none

  real,    intent (in)    :: semiAxis
  integer, intent (in)    :: seed
  integer, intent (in)    :: targetTotalGridPoints
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
end subroutine ed_beam2DGridPointsStatistical
