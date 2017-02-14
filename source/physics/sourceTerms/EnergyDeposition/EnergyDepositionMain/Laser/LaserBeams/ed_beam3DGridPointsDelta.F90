!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beam3DGridPointsDelta
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

  use Driver_interface,   ONLY : Driver_abortFlash

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
!     ...Check, if a grid is present at all. If not, then there are no grid points
!        and the program must stop.
!
!
  noGrid = (nTicsSemiAxisMajor == 0) .and. (nTicsSemiAxisMinor == 0)

  if (noGrid) then
      call Driver_abortFlash ("ed_beam3DGridPointsDelta: No grid points present!")
  end if
!
!
!     ...If requested, start the grid points.
!
!
  if (startGrid) then
      istart      = - nTicsSemiAxisMajor
      jstart      = - nTicsSemiAxisMinor
      deltaInv    = 1.0 / deltaSemiAxisMajor
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

     if (j < 0) then
         y = real (j+1) * deltaSemiAxisMinor - firstTicSemiAxisMinor
     else if (j > 0) then
         y = real (j-1) * deltaSemiAxisMinor + firstTicSemiAxisMinor
     else
         cycle
     end if

     x = aspectRatio * sqrt ((semiAxisMinor - y) * (semiAxisMinor + y))
     nTicsX = int ((x - firstTicSemiAxisMajor) * deltaInv) + 1
     istart = max (istart , - nTicsX)

     do i = istart, nTicsX

        if (i < 0) then
            x = real (i+1) * deltaSemiAxisMajor - firstTicSemiAxisMajor
        else if (i > 0) then
            x = real (i-1) * deltaSemiAxisMajor + firstTicSemiAxisMajor
        else
            cycle
        end if

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
end subroutine ed_beam3DGridPointsDelta
