!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beam3DGridPointsDelta
!!
!! NAME
!!
!!  ed_beam3DGridPointsDelta
!!
!! SYNOPSIS
!!
!!  call ed_beam3DGridPointsDelta (real,    intent (in)    :: semiAxisMajor,
!!                                 real,    intent (in)    :: semiAxisMinor,
!!                                 integer, intent (in)    :: nTicsSemiAxisMajor,
!!                                 integer, intent (in)    :: nTicsSemiAxisMinor,
!!                                 real,    intent (in)    :: deltaSemiAxisMajor,
!!                                 real,    intent (in)    :: deltaSemiAxisMinor,
!!                                 real,    intent (in)    :: firstTicSemiAxisMajor,
!!                                 real,    intent (in)    :: firstTicSemiAxisMinor,
!!                                 logical, intent (inout) :: startGrid,
!!                                 integer, intent (in)    :: maxGridPoints,
!!                                 logical, intent (out)   :: moreGridPoints,
!!                                 integer, intent (out)   :: nGridPoints,
!!                                 real,    intent (out)   :: xGrid (1:maxGridPoints),
!!                                 real,    intent (out)   :: yGrid (1:maxGridPoints))
!!
!! DESCRIPTION
!!
!!  Returns a collection of delta grid points for a 3D beam. The delta grid has to be
!!  set up beforehand to be able to use this routine. The total number of delta grid
!!  points of the delta grid is cut into several equally sized arrays, each of which
!!  is returned consecutively by a call to this routine.
!!
!! ARGUMENTS
!!
!!  semiAxisMajor         : the elliptical major semiaxis of the 3D beam
!!  semiAxisMinor         : the elliptical minor semiaxis of the 3D beam
!!  nTicsSemiAxisMajor    : # of grid positions along the major semiaxis
!!  nTicsSemiAxisMinor    : # of grid positions along the minor semiaxis
!!  deltaSemiAxisMajor    : the tic separation for the elliptical major semiaxis of the 3D beam
!!  deltaSemiAxisMinor    : the tic separation for the elliptical minor semiaxis of the 3D beam
!!  firstTicSemiAxisMajor : position of the 1st tic along the major semiaxis (in both directions)
!!  firstTicSemiAxisMinor : position of the 1st tic along the minor semiaxis (in both directions)
!!  startGrid             : if true, the grid points will start from the beginning 
!!  maxGridPoints         : the maximum number of grid points that can be returned
!!  moreGridPoints        : if true, more grid points are expected
!!  nGridPoints           : the actual number of grid points returned
!!  xGrid                 : the major semiaxis based coordinates of the grid points
!!  yGrid                 : the minor semiaxis based coordinates of the grid points
!!
!! NOTES
!!
!!***

subroutine ed_beam3DGridPointsDelta (semiAxisMajor,                      &
                                     semiAxisMinor,                      &
                                     nTicsSemiAxisMajor,                 &
                                     nTicsSemiAxisMinor,                 &
                                     deltaSemiAxisMajor,                 &
                                     deltaSemiAxisMinor,                 &
                                     firstTicSemiAxisMajor,              &
                                     firstTicSemiAxisMinor,              &
                                     startGrid,                          &
                                     maxGridPoints,                      &
                                                         moreGridPoints, &
                                                         nGridPoints,    &
                                                         xGrid,          &
                                                         yGrid           )

  implicit none

  real,    intent (in)    :: semiAxisMajor
  real,    intent (in)    :: semiAxisMinor
  integer, intent (in)    :: nTicsSemiAxisMajor
  integer, intent (in)    :: nTicsSemiAxisMinor
  real,    intent (in)    :: deltaSemiAxisMajor
  real,    intent (in)    :: deltaSemiAxisMinor
  real,    intent (in)    :: firstTicSemiAxisMajor
  real,    intent (in)    :: firstTicSemiAxisMinor
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
end subroutine ed_beam3DGridPointsDelta
