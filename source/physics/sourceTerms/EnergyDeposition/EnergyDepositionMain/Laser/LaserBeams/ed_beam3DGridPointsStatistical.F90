!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beam3DGridPointsStatistical
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
!!  Needs the random number generator as provided in the FLASH untilities directory.
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

  use Driver_interface,   ONLY : Driver_abortFlash
  use ut_randomInterface, ONLY : ut_randomNumber, &
                                 ut_randomSeed

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

  logical       :: fullGrid
  logical       :: inGrid

  integer, save :: collectedGridPoints     ! save for next call
  integer       :: n
  integer       :: seedArraySize

  real          :: notUsed
  real          :: x, y

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
  call random_number (xGrid)           ! get as many random numbers for x-coordinate as possible
  call random_number (yGrid)           ! get as many random numbers for y-coordinate as possible

  nGridPoints = 0

  do n = 1,maxGridPoints

     x = xGrid (n)                     ! random number between [0,1]
     y = yGrid (n)                     ! random number between [0,1]

     x = x + x - 1.0                   ! shift to random number between [-1,1]
     y = y + y - 1.0                   ! shift to random number between [-1,1]

     inGrid = (x * x + y * y <= 1.0)   ! check, if in elliptical grid

     if (inGrid) then

         nGridPoints = nGridPoints + 1

         xGrid (nGridPoints) = x
         yGrid (nGridPoints) = y

         collectedGridPoints = collectedGridPoints + 1

         fullGrid = (collectedGridPoints == targetTotalGridPoints)

         if (fullGrid) then

             xGrid (1:nGridPoints) = semiAxisMajor * xGrid (1:nGridPoints)
             yGrid (1:nGridPoints) = semiAxisMinor * yGrid (1:nGridPoints)

             moreGridPoints = .false.
             return
         end if

     end if

  end do
!
!
!     ...Return the collected grid points. If this part of the code is reached, there
!        are more grid points to be collected.
!
!
  xGrid (1:nGridPoints) = semiAxisMajor * xGrid (1:nGridPoints)
  yGrid (1:nGridPoints) = semiAxisMinor * yGrid (1:nGridPoints)

  moreGridPoints = .true.
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beam3DGridPointsStatistical
