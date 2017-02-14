!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beam2DGridPointsStatistical
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
!!  Needs the random number generator as provided in the FLASH untilities directory.
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

  use Driver_interface,   ONLY : Driver_abortFlash
  use ut_randomInterface, ONLY : ut_randomNumber, &
                                 ut_randomSeed

  implicit none

  real,    intent (in)    :: semiAxis
  integer, intent (in)    :: seed
  integer, intent (in)    :: targetTotalGridPoints
  logical, intent (inout) :: startGrid
  integer, intent (in)    :: maxGridPoints
  logical, intent (out)   :: moreGridPoints
  integer, intent (out)   :: nGridPoints
  real,    intent (out)   :: xGrid (1:maxGridPoints)

  integer, save :: collectedGridPoints     ! save for next call
  integer       :: remainingGridPoints
  integer       :: seedArraySize

  real          :: notUsed

  integer, allocatable :: seedArray (:)
!
!
!     ...If requested, start the statistical grid points using the seed value.
!
!
  if (startGrid) then

      call ut_randomSeed    (ut_size = seedArraySize) ! get compiler dependent seed array size
      allocate (seedArray   (1:seedArraySize))        ! allocate properly sized seed array
      seedArray (:) = seed                            ! set the seed array
      call ut_randomSeed    (ut_put  = seedArray)     ! same seed values ensure the same sequence of random #'s
      call ut_randomNumber  (notUsed)                 ! for consecutive seeds -> 1st random # almost identical 
      deallocate (seedArray)                          ! deallocate the seed array

      collectedGridPoints = 0 
      startGrid = .false.
  end if
!
!
!     ...Retrieve the current array of grid points.
!
!
  remainingGridPoints = targetTotalGridPoints - collectedGridPoints

  if (remainingGridPoints > 0) then

      nGridPoints = min (remainingGridPoints , maxGridPoints)

      call random_number (xGrid (1:nGridPoints))                 ! random numbers between [0,1]

      xGrid (1:nGridPoints) = 2 * xGrid (1:nGridPoints) - 1.0    ! shift to random numbers between [-1,1]
      xGrid (1:nGridPoints) = semiAxis * xGrid (1:nGridPoints)   ! semiaxis based coordinates

      collectedGridPoints = collectedGridPoints + nGridPoints    ! update # of grid points counter

      moreGridPoints = (collectedGridPoints < targetTotalGridPoints)
  else
      moreGridPoints = .false.
  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beam2DGridPointsStatistical
