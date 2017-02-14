!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beam3DGridPointsRecBeam
!!
!! NAME
!!
!!  ed_beam3DGridPointsRecBeam
!!
!! SYNOPSIS
!!
!!  call ed_beam3DGridPointsRecBeam (integer, intent (in)    :: nTicsSemiAxisMajor,
!!                                   integer, intent (in)    :: nTicsSemiAxisMinor,
!!                                   real,    intent (in)    :: deltaSemiAxisMajor,
!!                                   real,    intent (in)    :: deltaSemiAxisMinor,
!!                                   logical, intent (inout) :: startGrid,
!!                                   integer, intent (in)    :: maxGridPoints,
!!                                   logical, intent (out)   :: moreGridPoints,
!!                                   integer, intent (out)   :: nGridPoints,
!!                                   real,    intent (out)   :: xGrid (1:maxGridPoints),
!!                                   real,    intent (out)   :: yGrid (1:maxGridPoints))
!!
!! DESCRIPTION
!!
!!  Returns a collection of rectangular grid points for a rectangular 3D beam. The rectangular
!!  grid has to be set up beforehand to be able to use this routine. The total number of rectangular
!!  grid points of the rectangular grid is cut into several equally sized arrays, each of which
!!  is returned consecutively by a call to this routine.
!!
!! ARGUMENTS
!!
!!  nTicsSemiAxisMajor : # of grid positions along the major semiaxis
!!  nTicsSemiAxisMinor : # of grid positions along the minor semiaxis
!!  deltaSemiAxisMajor : the separation distance between two consecutive tics on the major semiaxis
!!  deltaSemiAxisMinor : the separation distance between two consecutive tics on the minor semiaxis
!!  startGrid          : if true, the grid points will start from the beginning 
!!  maxGridPoints      : the maximum number of grid points that can be returned
!!  moreGridPoints     : if true, more grid points are expected
!!  nGridPoints        : the actual number of grid points returned
!!  xGrid              : the major semiaxis based coordinates of the grid points
!!  yGrid              : the minor semiaxis based coordinates of the grid points
!!
!! NOTES
!!
!!***

subroutine ed_beam3DGridPointsRecBeam (nTicsSemiAxisMajor,                 &
                                       nTicsSemiAxisMinor,                 &
                                       deltaSemiAxisMajor,                 &
                                       deltaSemiAxisMinor,                 &
                                       startGrid,                          &
                                       maxGridPoints,                      &
                                                           moreGridPoints, &
                                                           nGridPoints,    &
                                                           xGrid,          &
                                                           yGrid           )

  implicit none

  integer, intent (in)    :: nTicsSemiAxisMajor
  integer, intent (in)    :: nTicsSemiAxisMinor
  real,    intent (in)    :: deltaSemiAxisMajor
  real,    intent (in)    :: deltaSemiAxisMinor
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
end subroutine ed_beam3DGridPointsRecBeam
