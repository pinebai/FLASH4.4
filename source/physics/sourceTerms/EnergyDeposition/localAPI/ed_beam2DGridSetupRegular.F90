!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beam2DGridSetupRegular
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

  implicit none

  real,    intent (in)   :: semiAxis
  integer, intent (in)   :: nGridPoints
  integer, intent (out)  :: nTics
  real,    intent (out)  :: delta
  real,    intent (out)  :: firstTic

  nTics    = 0
  delta    = 0.0
  firstTic = 0.0

  return
end subroutine ed_beam2DGridSetupRegular
