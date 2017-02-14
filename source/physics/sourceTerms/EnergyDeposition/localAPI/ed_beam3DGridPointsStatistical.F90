!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beam3DGridPointsStatistical
!!
!! NAME
!!
!!  ed_beam3DGridPointsStatistical
!!
!! SYNOPSIS
!!
!!  call ed_beam3DGridPointsStatistical (real,    intent (in)    :: semiAxisMajor,
!!                                       real,    intent (in)    :: semiAxisMinor,
!!                                       real,    intent (in)    :: seed,
!!                                       integer, intent (in)    :: targetTotalGridPoints,
!!                                       logical, intent (inout) :: startGrid,
!!                                       integer, intent (in)    :: maxGridPoints,
!!                                       logical, intent (out)   :: moreGridPoints,
!!                                       integer, intent (out)   :: nGridPoints,
!!                                       real,    intent (out)   :: xGrid (1:maxGridPoints),
!!                                       real,    intent (out)   :: yGrid (1:maxGridPoints))
!!
!! DESCRIPTION
!!
!!  Returns a collection of statistical grid points for a 3D beam. The statistical grid has
!!  to be set up beforehand to be able to use this routine. The total number of statistical
!!  grid points of the statistical grid is cut into several equally sized arrays, each of which
!!  is returned consecutively by a call to this routine.
!!
!! ARGUMENTS
!!
!!  semiAxisMajor         : the elliptical major semiaxis of the 3D beam
!!  semiAxisMinor         : the elliptical minor semiaxis of the 3D beam
!!  seed                  : the seed value defining the statistical grid (random number sequence)
!!  targetTotalGridPoints : the complete number of grid points wanted
!!  startGrid             : if true, the grid points will start from the beginning 
!!  maxGridPoints         : the maximum number of grid points that can be returned
!!  moreGridPoints        : if true, more grid points are expected
!!  nGridPoints           : the actual number of grid points returned
!!  xGrid                 : the major semiaxis based coordinates of the grid points
!!  yGrid                 : the minor semiaxis based coordinates of the grid points
!!
!! NOTES
!!
!!  The current algorithm uses the arrays reserved for the final grid points for harvesting
!!  the random numbers needed. Since some of these numbers lay not within the ellipsoidal grid,
!!  the number of accepted grid points per random number harvesting step is always approximately
!!  pi/4 or 3/4 less than the maximum number of random points harvested. In order for this
!!  procedure to be efficient in collecting all the grid points, the size of the grid arrays
!!  should not be too small.
!!
!!***

subroutine ed_beam3DGridPointsStatistical (semiAxisMajor,                      &
                                           semiAxisMinor,                      &
                                           seed,                               &
                                           targetTotalGridPoints,              &
                                           startGrid,                          &
                                           maxGridPoints,                      &
                                                               moreGridPoints, &
                                                               nGridPoints,    &
                                                               xGrid,          &
                                                               yGrid           )

  implicit none

  real,    intent (in)    :: semiAxisMajor
  real,    intent (in)    :: semiAxisMinor
  integer, intent (in)    :: seed
  integer, intent (in)    :: targetTotalGridPoints
  logical, intent (inout) :: startGrid
  integer, intent (in)    :: maxGridPoints
  logical, intent (out)   :: moreGridPoints
  integer, intent (out)   :: nGridPoints
  real,    intent (out)   :: xGrid (1:maxGridPoints)
  real,    intent (out)   :: yGrid (1:maxGridPoints)

  startGrid               = .false.
  moreGridPoints          = .false.
  nGridPoints             =  0
  xGrid (1:maxGridPoints) =  0.0
  yGrid (1:maxGridPoints) =  0.0

  return
end subroutine ed_beam3DGridPointsStatistical
