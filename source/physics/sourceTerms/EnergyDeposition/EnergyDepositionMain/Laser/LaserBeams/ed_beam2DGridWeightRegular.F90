!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beam2DGridWeightRegular
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

  use Driver_interface,             ONLY : Driver_abortFlash
  use ed_interface,                 ONLY : ed_beam2DGridPointsRegular
  use ed_beamCrossSectionFunctions, ONLY : ed_beamCrossSectionWeight
  use EnergyDeposition_data,        ONLY : ed_gridGeometry

  implicit none

#include "EnergyDeposition.h"

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

  logical :: grid2Dcylindrical
  logical :: moreGridPoints
  logical :: startGrid

  integer :: maxGridPoints
  integer :: n
  integer :: nGridPoints
  integer :: usedGridPoints

  real :: weightRadial, weightCrossSection
  real :: x

  real, allocatable :: xGrid (:)
!
!
!     ...Set geometrical flag.
!
!
  grid2Dcylindrical = (ed_gridGeometry == GRID_2DCYLINDRICAL)
!
!
!     ...Calculate the grid weight.
!
!
  gridWeight     =  0.0

  startGrid      = .true.
  moreGridPoints = .true.
  usedGridPoints =  0
  maxGridPoints  =  BEAM_GRID_ARRAYSIZE

  allocate (xGrid (1:maxGridPoints))

  do while (moreGridPoints)

     call ed_beam2DGridPointsRegular (semiAxis,                      &
                                      nTics,                         &
                                      delta,                         &
                                      firstTic,                      &
                                      startGrid,                     &
                                      maxGridPoints,                 &
                                                     moreGridPoints, &
                                                     nGridPoints,    &
                                                     xGrid           )
!
!
!     ...For radial grids, adjust the weight such that the resulting beam intensity (power/area)
!        profile matches the cross sectional function ID given. Since for radial grids the
!        density of grid points (# of grid points / area) reduces proportional to the radius, each
!        grid position has to be weighted by the radius. 
!
!
     if (grid2Dcylindrical .AND. target2LensMagnify*ux.NE.0.0) then

         do n = 1, nGridPoints

            x = xGrid (n)

            weightRadial = abs (gridRadialOrigin + x*target2LensMagnify*ux)

            weightCrossSection = ed_beamCrossSectionWeight (crossSectionFunctionType,   &
                                                            x        = x,               &
                                                            Cx       = gaussianCenter,  &
                                                            Rx       = gaussianRadius,  &
                                                            Exponent = gaussianExponent )

            gridWeight = gridWeight + weightRadial * weightCrossSection

         end do

     else

         do n = 1, nGridPoints

            x = xGrid (n)

            weightCrossSection = ed_beamCrossSectionWeight (crossSectionFunctionType,   &
                                                            x        = x,               &
                                                            Cx       = gaussianCenter,  &
                                                            Rx       = gaussianRadius,  &
                                                            Exponent = gaussianExponent )

            gridWeight = gridWeight + weightCrossSection

         end do

     end if

     usedGridPoints = usedGridPoints + nGridPoints

  end do

  if (usedGridPoints /= totalGridPoints) then
      call Driver_abortFlash ("ed_beam2DGridWeightRegular: # of used and total grid points mismatch !")
  end if

  deallocate (xGrid)
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beam2DGridWeightRegular
