!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beam2DGridWeightRegular
!!
!! NAME
!!
!!  ed_beam2DGridWeightRegular
!!
!! SYNOPSIS
!!
!!  call ed_beam2DGridWeightRegular (real,              intent (in)  :: gridRadialOrigin,
!!                                   real,              intent (in)  :: ux,
!!                                   real,              intent (in)  :: target2LensMagnify,
!!                                   real,              intent (in)  :: semiAxis,
!!                                   character (len=*), intent (in)  :: crossSectionFunctionType,
!!                                   real,              intent (in)  :: gaussianExponent,
!!                                   real,              intent (in)  :: gaussianRadius,
!!                                   real,              intent (in)  :: gaussianCenter,
!!                                   integer,           intent (in)  :: nTics,
!!                                   real,              intent (in)  :: delta,
!!                                   real,              intent (in)  :: firstTic,
!!                                   integer,           intent (in)  :: totalGridPoints,
!!                                   real,              intent (out) :: gridWeight)
!!
!! DESCRIPTION
!!
!!  Calculates the weight for the linear 2D beam regular grid, which is defined as the sum of the
!!  individual regular grid point weights. The regular grid must have been set up for the grid points
!!  to be ready for retrieval. If the number of grid points processed mismatches the total number of
!!  grid points of the beam grid, the program is aborted with a message.
!!
!!  For 2D cylindrical geometries there is an additional radial weight to ensure proper beam
!!  intensities (power / area). The radial weight is multiplied with the beam's cross sectional
!!  weight and consists of two components: 1) absolute radial position of the beam's grid origin
!!  plus 2) the local grid point locations on the grid.
!!
!! ARGUMENTS
!!
!!  gridRadialOrigin         : the (absolute) radial location of the grid's origin
!!                             (only relevant for 2D cylindrical geometries)
!!  ux                       : radial component of a unit vector that ponts along the grid
!!                             (only relevant for 2D cylindrical geometries)
!!  target2LensMagnify       : factor that scales target size to "lens" size 
!!                             (only used for 2D cylindrical geometries)
!!  semiAxis                 : the semiaxis defining the linear cross section of the 2D beam
!!  crossSectionFunctionType : the cross section function type defining the grid point weigths
!!  gaussianExponent         : the gaussian exponent
!!  gaussianRadius           : the gaussian radius along the semiaxis
!!  gaussianCenter           : the gaussian center location along the semiaxis 
!!  nTics                    : the total number of tics on the grid
!!  delta                    : the tic spacing (for regular grids)
!!  firstTic                 : the position of the 1st tic from the lowest boundary
!!  totalGridPoints          : the total number of grid points of the beam grid
!!  gridWeight               : the value of the grid weight
!!
!! NOTES
!!
!!***

subroutine ed_beam2DGridWeightRegular (gridRadialOrigin,                  &
                                       ux,                                &
                                       target2LensMagnify,                &
                                       semiAxis,                          &
                                       crossSectionFunctionType,          &
                                       gaussianExponent,                  &
                                       gaussianRadius,                    &
                                       gaussianCenter,                    &
                                       nTics,                             &
                                       delta,                             &
                                       firstTic,                          &
                                       totalGridPoints,                   &
                                                               gridWeight )

  implicit none

  real,              intent (in)  :: gridRadialOrigin
  real,              intent (in)  :: ux
  real,              intent (in)  :: target2LensMagnify
  real,              intent (in)  :: semiAxis
  character (len=*), intent (in)  :: crossSectionFunctionType
  real,              intent (in)  :: gaussianExponent
  real,              intent (in)  :: gaussianRadius
  real,              intent (in)  :: gaussianCenter
  integer,           intent (in)  :: nTics
  real,              intent (in)  :: delta
  real,              intent (in)  :: firstTic
  integer,           intent (in)  :: totalGridPoints
  real,              intent (out) :: gridWeight

  gridWeight = 0.0

  return
end subroutine ed_beam2DGridWeightRegular
