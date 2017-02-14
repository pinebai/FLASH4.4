!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beam3DGridSetupRadial
!!
!! NAME
!!
!!  ed_beam3DGridSetupRadial
!!
!! SYNOPSIS
!!
!!  call ed_beam3DGridSetupRadial (integer, intent (inout) :: nGridPoints,
!!                                 integer, intent (inout) :: nTicsRadial,
!!                                 integer, intent (inout) :: nTicsAngular,
!!                                 real,    intent (out)   :: deltaRadial,
!!                                 real,    intent (out)   :: deltaAngular,
!!                                 real,    intent (out)   :: firstTicRadial,
!!                                 real,    intent (out)   :: firstTicAngular)
!!
!! DESCRIPTION
!!
!!  Sets up information about the radial cross sectional grid for a particular 3D beam.
!!  On input, the number of grid points is the wanted number of grid points. On output, the
!!  number of grid points might have changed and it overrides the input value.
!!
!!  The radial grid consists of a collection of same sized angular spaced radial spikes and
!!  equidistant concentrical ellipses inside the boundary beam cross section ellipse.
!!  The collection of the concentrical ellipses gives rise to equidistant grid positions
!!  (radial tics) on each radial spike.
!!
!!  There is the option of passing predefined values for the number of concentrical
!!  ellipses (nTicsRadial) and the number of radial spikes (nTicsAngular). The following
!!  actions are being taken, depending on the input values of these two grid parameters:
!!
!!                 nTicsRadial  =< 0   ->   automatic evaluation
!!                 nTicsAngular =< 0   ->   automatic evaluation
!!
!!                 nTicsRadial   > 0   ->   take this value
!!                 nTicsAngular =< 0   ->   adjust to number of grid points wanted
!!
!!                 nTicsRadial  =< 0   ->   adjust to number of grid points wanted
!!                 nTicsAngular  > 0   ->   take this value
!!
!!                 nTicsRadial   > 0   ->   take this value
!!                 nTicsAngular  > 0   ->   take this value
!!
!!
!! ARGUMENTS
!!
!!  nGridPoints     : on input  -> the number of grid points wanted
!!                    on output -> the optimum number of grid points close to the input value
!!  nTicsRadial     : # of grid positions along each radial spike (if > 0 -> assumed user given)
!!  nTicsAngular    : # of angular slices along the angular dimension (if > 0 -> assumed user given)
!!  deltaRadial     : the tic spacing along each radial spike (as fraction of radius)
!!  deltaAngular    : the angular spacing along the angular dimension (as fraction of 2pi)
!!  firstTicRadial  : position of 1st tic along each radial spike (as fraction of radius)
!!  firstTicAngular : position of 1st tic along the angular dimension (as fraction of 2pi)
!!
!! NOTES
!!
!!  Since the radial grid is defined in terms of concentric inner ellipses, there is no
!!  need to specify the boundary ellipses detailed structure, i.e. the length of both its
!!  semiaxes. This contrasts the case of the square grid, which needs at least the ratio
!!  between both semiaxes in order to establish the square grid.
!!
!!  The structure of the radial grid is currently set in such a way that the innermost
!!  ellipse is 1/2 radial tic spacing away from the center and there is no grid point
!!  at the center. The first angular tic is currently set to 0.
!!
!!***

subroutine ed_beam3DGridSetupRadial (nGridPoints,    &
                                     nTicsRadial,    &
                                     nTicsAngular,   &
                                     deltaRadial,    &
                                     deltaAngular,   &
                                     firstTicRadial, &
                                     firstTicAngular )

  implicit none

  integer, intent (inout) :: nGridPoints
  integer, intent (inout) :: nTicsRadial
  integer, intent (inout) :: nTicsAngular
  real,    intent (out)   :: deltaRadial
  real,    intent (out)   :: deltaAngular
  real,    intent (out)   :: firstTicRadial
  real,    intent (out)   :: firstTicAngular

  nGridPoints     = 0
  nTicsRadial     = 0
  nTicsAngular    = 0
  deltaRadial     = 0.0
  deltaAngular    = 0.0
  firstTicRadial  = 0.0
  firstTicAngular = 0.0

  return
end subroutine ed_beam3DGridSetupRadial
