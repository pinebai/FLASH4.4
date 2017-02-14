!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beam2DGridSetupRegular
!!
!! NAME
!!
!!  ed_beam2DGridSetupRegular
!!
!! SYNOPSIS
!!
!!  call ed_beam2DGridSetupRegular (real,    intent (in)  :: semiAxis,
!!                                  integer, intent (in)  :: nGridPoints,
!!                                  integer, intent (out) :: nTics,
!!                                  real,    intent (out) :: delta,
!!                                  real,    intent (out) :: firstTic)
!!
!! DESCRIPTION
!!
!!  Generates the regular 1-dimensional grid for a 2D beam. The regular grid is defined
!!  through a delta value (separation of consecutive tics) and the position of the first
!!  tic on the beams cross sectional line (shown for 9 grid points):
!!
!!
!!                         1st tic
!!                            |
!!                            1   2   3   4   5   6   7   8   9
!!                          |-|---|---|---|---|---|---|---|---|-|
!!
!!       lowest boundary -> |  <--- cross section length ---->  |
!!                                            | <- semiaxis --> |
!!
!!                        
!!  The number of grid points can always be honored exactly for 1-dimensional linear grids.
!!
!! ARGUMENTS
!!
!!  semiAxis    : the semiaxis (1/2 length) of the 2D beam cross sectional area
!!  nGridPoints : the number of grid points wanted
!!  nTics       : the total number of tics on the grid (equal to the number of grid points)
!!  delta       : the tic spacing (for regular grids)
!!  firstTic    : the position of the 1st tic from the lowest boundary
!!
!! NOTES
!!
!!  The regular grid is currently being defined as a 1-dimensional grid having its first
!!  and last tics 1/2 the delta value from the grid boundary.
!!
!!***

subroutine ed_beam2DGridSetupRegular (semiAxis,             &
                                      nGridPoints,          &
                                                   nTics,   &
                                                   delta,   &
                                                   firstTic )

  use ed_interface, ONLY : ed_gridLineRegular

  implicit none

  real,    intent (in)   :: semiAxis
  integer, intent (in)   :: nGridPoints
  integer, intent (out)  :: nTics
  real,    intent (out)  :: delta
  real,    intent (out)  :: firstTic

  real :: crossSectionLength
!
!
!     ...Set up the regular grid. The routine below sets the tics evenly between
!        the grid boundaries, including two tics on the boundary. To get the specific
!        regular grid mentioned above, we call this routine with one extra grid point
!        (to get the right delta value) and set the frist tic to 1/2 delta. Note, that
!        this covers the 1 grid point case as well.
!
!
  crossSectionLength = semiAxis + semiAxis

!  call ed_gridLineRegular (crossSectionLength,       &
!                           nGridPoints,              &
!                                               delta )
!
!  nTics    = nGridPoints
!  firstTic = 0.0

  call ed_gridLineRegular (crossSectionLength,       &
                           nGridPoints + 1,          &
                                               delta )

  nTics    = nGridPoints
  firstTic = 0.5 * delta
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beam2DGridSetupRegular
