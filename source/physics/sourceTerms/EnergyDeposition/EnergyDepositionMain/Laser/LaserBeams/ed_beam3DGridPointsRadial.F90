!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beam3DGridPointsRadial
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

  use Driver_interface,   ONLY : Driver_abortFlash

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

  logical       :: noGrid

  integer       :: i, j
  integer, save :: istart                  ! save radial grid starting point for next call
  integer, save :: jstart                  ! save angular grid starting point for next call

  real,    save :: alpha, beta
  real,    save :: angleSine, angleCosine
  real          :: angleSineSave
  real          :: halfSine
  real          :: radialFraction
  real          :: x, y
!
!
!     ...Check, if a grid is present at all. If not, then there is only 1 grid point at
!        the center of the beam.
!
!
  noGrid = (nTicsRadial == 0) .and. (nTicsAngular == 0)

  if (noGrid) then
      xGrid (1)      = 0.0
      yGrid (1)      = 0.0
      nGridPoints    = 1
      moreGridPoints = .false.
      return
  end if
!
!
!     ...If requested, start the grid points.
!
!
  if (startGrid) then
      istart      = 1
      jstart      = 1
      halfSine    = sin (0.5 * deltaAngular)
      alpha       = 2.0 * halfSine * halfSine
      beta        = sin (deltaAngular)
      angleSine   = sin (firstTicAngular)
      angleCosine = cos (firstTicAngular)
      startGrid   = .false.
  end if
!
!
!     ...Retrieve the current array of grid points.
!
!
  nGridPoints = 0

  do j = jstart, nTicsAngular

     radialFraction = firstTicRadial + (istart - 1) * deltaRadial

     do i = istart, nTicsRadial

        x = semiAxisMajor * radialFraction * angleCosine        ! angle measured from major semiaxis
        y = semiAxisMinor * radialFraction * angleSine          ! in counterclockwise direction

        nGridPoints = nGridPoints + 1

        if (nGridPoints > maxGridPoints) then
            istart = i                     ! save current radial grid starting point for next call
            jstart = j                     ! save current angular grid starting point for next call
            nGridPoints = nGridPoints - 1  ! adjust to the true # of grid points returned
            moreGridPoints = .true.        ! this is always the case
            return
        else
            xGrid (nGridPoints) = x
            yGrid (nGridPoints) = y
            moreGridPoints = .not. ((i == nTicsRadial) .and. (j == nTicsAngular))
        end if

        radialFraction = radialFraction + deltaRadial

     end do

     angleSineSave = angleSine
     angleSine     = angleSine   - (alpha * angleSine   - beta * angleCosine  )
     angleCosine   = angleCosine - (alpha * angleCosine + beta * angleSineSave)

     istart = 1

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beam3DGridPointsRadial
