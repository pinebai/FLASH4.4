!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beam3DGridPointsRadial
!!
!! NAME
!!
!!  ed_beam3DGridPointsRadial
!!
!! SYNOPSIS
!!
!!  call ed_beam3DGridPointsRadial (real,    intent (in)    :: semiAxisMajor,
!!                                  real,    intent (in)    :: semiAxisMinor,
!!                                  integer, intent (in)    :: nTicsRadial,
!!                                  integer, intent (in)    :: nTicsAngular,
!!                                  real,    intent (in)    :: deltaRadial,
!!                                  real,    intent (in)    :: deltaAngular,
!!                                  real,    intent (in)    :: firstTicRadial,
!!                                  real,    intent (in)    :: firstTicAngular,
!!                                  logical, intent (inout) :: startGrid,
!!                                  integer, intent (in)    :: maxGridPoints,
!!                                  logical, intent (out)   :: moreGridPoints,
!!                                  integer, intent (out)   :: nGridPoints,
!!                                  real,    intent (out)   :: xGrid (1:maxGridPoints),
!!                                  real,    intent (out)   :: yGrid (1:maxGridPoints))
!!
!! DESCRIPTION
!!
!!  Returns a collection of radial grid points for a 3D beam. The radial grid has to be
!!  set up beforehand to be able to use this routine. The total number of radial grid
!!  points of the radial grid is cut into several equally sized arrays, each of which
!!  is returned consecutively by a call to this routine.
!!
!! ARGUMENTS
!!
!!  semiAxisMajor   : the elliptical major semiaxis of the 3D beam
!!  semiAxisMinor   : the elliptical minor semiaxis of the 3D beam
!!  nTicsRadial     : # of grid positions along each radial spike
!!  nTicsAngular    : # of angular slices along the angular dimension
!!  deltaRadial     : the tic spacing along each radial spike (as fraction of radius)
!!  deltaAngular    : the angular spacing along the angular dimension (as fraction of 2pi)
!!  firstTicRadial  : position of 1st tic along each radial spike (as fraction of radius)
!!  firstTicAngular : position of 1st tic along the angular dimension (as fraction of 2pi)
!!  startGrid       : if true, the grid points will start from the beginning 
!!  maxGridPoints   : the maximum number of grid points that can be returned
!!  moreGridPoints  : if true, more grid points are expected
!!  nGridPoints     : the actual number of grid points returned
!!  xGrid           : the major semiaxis based coordinates of the grid points
!!  yGrid           : the minor semiaxis based coordinates of the grid points
!!
!! NOTES
!!
!!  The angles are measured from the major semiaxis in counterclockwise direction. The major
!!  semiaxis constitutes thus the angular basis from which all angles will be measured.
!!  It is in this routine where the abstract definition of the radial grid is given a
!!  definite angular orientation inside the beam.
!!
!!***

subroutine ed_beam3DGridPointsRadial (semiAxisMajor,                   &
                                      semiAxisMinor,                   &
                                      nTicsRadial,                     &
                                      nTicsAngular,                    &
                                      deltaRadial,                     &
                                      deltaAngular,                    &
                                      firstTicRadial,                  &
                                      firstTicAngular,                 &
                                      startGrid,                       &
                                      maxGridPoints,                   &
                                                       moreGridPoints, &
                                                       nGridPoints,    &
                                                       xGrid,          &
                                                       yGrid           )

  implicit none

  real,    intent (in)    :: semiAxisMajor
  real,    intent (in)    :: semiAxisMinor
  integer, intent (in)    :: nTicsRadial
  integer, intent (in)    :: nTicsAngular
  real,    intent (in)    :: deltaRadial
  real,    intent (in)    :: deltaAngular
  real,    intent (in)    :: firstTicRadial
  real,    intent (in)    :: firstTicAngular
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
end subroutine ed_beam3DGridPointsRadial
