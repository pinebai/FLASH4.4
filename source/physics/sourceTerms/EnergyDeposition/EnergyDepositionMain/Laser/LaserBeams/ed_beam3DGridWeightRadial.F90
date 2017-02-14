!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beam3DGridWeightRadial
!!
!! NAME
!!
!!  ed_beam3DGridWeightRadial
!!
!! SYNOPSIS
!!
!!  call ed_beam3DGridWeightRadial (real,              intent (in)  :: semiAxisMajor,
!!                                  real,              intent (in)  :: semiAxisMinor,
!!                                  character (len=*), intent (in)  :: crossSectionFunctionType,
!!                                  real,              intent (in)  :: gaussianExponent,
!!                                  real,              intent (in)  :: gaussianRadiusMajor,
!!                                  real,              intent (in)  :: gaussianRadiusMinor,
!!                                  real,              intent (in)  :: gaussianCenterMajor,
!!                                  real,              intent (in)  :: gaussianCenterMinor,
!!                                  integer,           intent (in)  :: nTicsRadial,
!!                                  integer,           intent (in)  :: nTicsAngular,
!!                                  real,              intent (in)  :: deltaRadial,
!!                                  real,              intent (in)  :: deltaAngular,
!!                                  real,              intent (in)  :: firstTicRadial,
!!                                  real,              intent (in)  :: firstTicAngular,
!!                                  integer,           intent (in)  :: totalGridPoints,
!!                                  real,              intent (out) :: gridWeight)
!!
!! DESCRIPTION
!!
!!  Calculates the weight for the 3D beam radial grid, which is defined as the sum of the individual
!!  radial grid point weights. The radial grid must have been set up for the grid points to be ready
!!  for retrieval. If the number of grid points processed mismatches the total number of grid points
!!  of the beam grid, the program is aborted with a message.
!!
!! ARGUMENTS
!!
!!  semiAxisMajor            : the elliptical major semiaxis of the 3D beam
!!  semiAxisMinor            : the elliptical minor semiaxis of the 3D beam 
!!  crossSectionFunctionType : the cross section function type defining the grid point weigths
!!  gaussianExponent         : the gaussian exponent
!!  gaussianRadiusMajor      : the gaussian radius along the major semiaxis
!!  gaussianRadiusMinor      : the gaussian radius along the minor semiaxis
!!  gaussianCenterMajor      : the gaussian center location along the major semiaxis
!!  gaussianCenterMinor      : the gaussian center location along the minor semiaxis
!!  nTicsRadial              : # of grid positions along each radial spike
!!  nTicsAngular             : # of angular slices along the angular dimension
!!  deltaRadial              : the tic spacing along each radial spike (as fraction of radius)
!!  deltaAngular             : the angular spacing along the angular dimension (as fraction of 2pi)
!!  firstTicRadial           : position of 1st tic along each radial spike (as fraction of radius)
!!  firstTicAngular          : position of 1st tic along the angular dimension (as fraction of 2pi)
!!  totalGridPoints          : the total number of grid points of the beam grid
!!  gridWeight               : the value of the grid weight
!!
!! NOTES
!!
!!***

subroutine ed_beam3DGridWeightRadial (semiAxisMajor,                    &
                                      semiAxisMinor,                    &
                                      crossSectionFunctionType,         &
                                      gaussianExponent,                 &
                                      gaussianRadiusMajor,              &
                                      gaussianRadiusMinor,              &
                                      gaussianCenterMajor,              &
                                      gaussianCenterMinor,              &
                                      nTicsRadial,                      &
                                      nTicsAngular,                     &
                                      deltaRadial,                      &
                                      deltaAngular,                     &
                                      firstTicRadial,                   &
                                      firstTicAngular,                  &
                                      totalGridPoints,                  &
                                                             gridWeight )

  use Driver_interface,             ONLY : Driver_abortFlash
  use ed_interface,                 ONLY : ed_beam3DGridPointsRadial
  use ed_beamCrossSectionFunctions, ONLY : ed_beamCrossSectionWeight

  implicit none

#include "EnergyDeposition.h"

  real,              intent (in)  :: semiAxisMajor
  real,              intent (in)  :: semiAxisMinor
  character (len=*), intent (in)  :: crossSectionFunctionType
  real,              intent (in)  :: gaussianExponent
  real,              intent (in)  :: gaussianRadiusMajor
  real,              intent (in)  :: gaussianRadiusMinor
  real,              intent (in)  :: gaussianCenterMajor
  real,              intent (in)  :: gaussianCenterMinor
  integer,           intent (in)  :: nTicsRadial
  integer,           intent (in)  :: nTicsAngular
  real,              intent (in)  :: deltaRadial
  real,              intent (in)  :: deltaAngular
  real,              intent (in)  :: firstTicRadial
  real,              intent (in)  :: firstTicAngular
  integer,           intent (in)  :: totalGridPoints
  real,              intent (out) :: gridWeight

  logical :: moreGridPoints
  logical :: startGrid

  integer :: maxGridPoints
  integer :: n
  integer :: nGridPoints
  integer :: usedGridPoints

  real :: aspectRatio
  real :: weightRadial, weightCrossSection
  real :: x,y

  real, allocatable :: xGrid (:)
  real, allocatable :: yGrid (:)
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
  allocate (yGrid (1:maxGridPoints))

  do while (moreGridPoints)

     call ed_beam3DGridPointsRadial (semiAxisMajor,                 &
                                     semiAxisMinor,                 &
                                     nTicsRadial,                   &
                                     nTicsAngular,                  &
                                     deltaRadial,                   &
                                     deltaAngular,                  &
                                     firstTicRadial,                &
                                     firstTicAngular,               &
                                     startGrid,                     &
                                     maxGridPoints,                 &
                                                    moreGridPoints, &
                                                    nGridPoints,    &
                                                    xGrid,          &
                                                    yGrid           )
!
!
!     ...For radial grids, adjust the weight such that the resulting beam intensity (power/area)
!        profile matches the cross sectional function type given. Since for radial grids the
!        density of grid points (# of grid points / area) reduces proportional to the ellipse radius,
!        each grid position has to be weighted by this radius. 
!
!
     aspectRatio = semiAxisMajor / semiAxisMinor

     do n = 1, nGridPoints

        x = xGrid (n)
        y = yGrid (n)

        weightCrossSection = ed_beamCrossSectionWeight (crossSectionFunctionType,       &
                                                        x        = x,                   &
                                                        y        = y,                   &
                                                        Cx       = gaussianCenterMajor, &
                                                        Cy       = gaussianCenterMinor, &
                                                        Rx       = gaussianRadiusMajor, &
                                                        Ry       = gaussianRadiusMinor, &
                                                        Exponent = gaussianExponent     )

        y = y * aspectRatio

        weightRadial = sqrt (x * x + y * y)

        gridWeight = gridWeight + weightRadial * weightCrossSection

     end do

     usedGridPoints = usedGridPoints + nGridPoints

  end do

  if (usedGridPoints /= totalGridPoints) then
      call Driver_abortFlash ("ed_beam3DGridWeightRadial: # of used and total grid points mismatch !")
  end if

  deallocate (xGrid)
  deallocate (yGrid)
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beam3DGridWeightRadial
