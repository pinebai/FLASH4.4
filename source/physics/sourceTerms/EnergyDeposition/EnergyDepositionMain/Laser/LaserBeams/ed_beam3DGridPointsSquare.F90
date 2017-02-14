!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beam3DGridPointsSquare
!!
!! NAME
!!
!!  ed_beam3DGridPointsSquare
!!
!! SYNOPSIS
!!
!!  call ed_beam3DGridPointsSquare (real,    intent (in)    :: semiAxisMajor,
!!                                  real,    intent (in)    :: semiAxisMinor,
!!                                  integer, intent (in)    :: nTicsSemiAxisMajor,
!!                                  integer, intent (in)    :: nTicsSemiAxisMinor,
!!                                  real,    intent (in)    :: delta,
!!                                  logical, intent (inout) :: startGrid,
!!                                  integer, intent (in)    :: maxGridPoints,
!!                                  logical, intent (out)   :: moreGridPoints,
!!                                  integer, intent (out)   :: nGridPoints,
!!                                  real,    intent (out)   :: xGrid (1:maxGridPoints),
!!                                  real,    intent (out)   :: yGrid (1:maxGridPoints))
!!
!! DESCRIPTION
!!
!!  Returns a collection of square grid points for a 3D beam. The square grid has to be
!!  set up beforehand to be able to use this routine. The total number of square grid
!!  points of the square grid is cut into several equally sized arrays, each of which
!!  is returned consecutively by a call to this routine.
!!
!! ARGUMENTS
!!
!!  semiAxisMajor      : the elliptical major semiaxis of the 3D beam
!!  semiAxisMinor      : the elliptical minor semiaxis of the 3D beam
!!  nTicsSemiAxisMajor : # of grid positions along the major semiaxis
!!  nTicsSemiAxisMinor : # of grid positions along the minor semiaxis
!!  delta              : the separation distance between two consecutive tics
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

subroutine ed_beam3DGridPointsSquare (semiAxisMajor,                      &
                                      semiAxisMinor,                      &
                                      nTicsSemiAxisMajor,                 &
                                      nTicsSemiAxisMinor,                 &
                                      delta,                              &
                                      startGrid,                          &
                                      maxGridPoints,                      &
                                                          moreGridPoints, &
                                                          nGridPoints,    &
                                                          xGrid,          &
                                                          yGrid           )

  use Driver_interface,   ONLY : Driver_abortFlash

  implicit none

  real,    intent (in)    :: semiAxisMajor
  real,    intent (in)    :: semiAxisMinor
  integer, intent (in)    :: nTicsSemiAxisMajor
  integer, intent (in)    :: nTicsSemiAxisMinor
  real,    intent (in)    :: delta
  logical, intent (inout) :: startGrid
  integer, intent (in)    :: maxGridPoints
  logical, intent (out)   :: moreGridPoints
  integer, intent (out)   :: nGridPoints
  real,    intent (out)   :: xGrid (1:maxGridPoints)
  real,    intent (out)   :: yGrid (1:maxGridPoints)

  logical       :: noGrid

  integer       :: i, j
  integer, save :: istart                  ! save major semiaxis grid starting point for next call
  integer, save :: jstart                  ! save minor semiaxis grid starting point for next call
  integer       :: nTicsX

  real,    save :: aspectRatio
  real,    save :: deltaInv
  real          :: x, y
!
!
!     ...Check, if a grid is present at all. If not, then there is only 1 grid point at
!        the center of the beam.
!
!
  noGrid = (nTicsSemiAxisMajor == 0) .and. (nTicsSemiAxisMinor == 0)

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
      istart      = - nTicsSemiAxisMajor
      jstart      = - nTicsSemiAxisMinor
      deltaInv    = 1.0 / delta
      aspectRatio = semiAxisMajor / semiAxisMinor
      startGrid   = .false.
  end if
!
!
!     ...Retrieve the current array of grid points.
!
!
  nGridPoints = 0

  do j = jstart, nTicsSemiAxisMinor
     y = real (j) * delta
     x = aspectRatio * sqrt ((semiAxisMinor - y) * (semiAxisMinor + y))
     nTicsX = min (nTicsSemiAxisMajor , int (x * deltaInv))
     istart = max (istart , - nTicsX)

     do i = istart, nTicsX
        x = real (i) * delta
        nGridPoints = nGridPoints + 1

        if (nGridPoints > maxGridPoints) then
            istart = i                     ! save current major semiaxis grid starting point for next call
            jstart = j                     ! save current minor semiaxis grid starting point for next call
            nGridPoints = nGridPoints - 1  ! adjust to the true # of grid points returned
            moreGridPoints = .true.        ! this is always the case
            return
        else
            xGrid (nGridPoints) = x
            yGrid (nGridPoints) = y
            moreGridPoints = .not. ((i == nTicsX) .and. (j == nTicsSemiAxisMinor))
        end if

     end do
     istart = - nTicsSemiAxisMajor
  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beam3DGridPointsSquare
