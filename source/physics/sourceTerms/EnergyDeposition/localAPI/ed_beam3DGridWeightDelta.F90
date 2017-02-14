!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beam3DGridWeightDelta
!!
!! NAME
!!
!!  ed_beam3DGridWeightDelta
!!
!! SYNOPSIS
!!
!!  call ed_beam3DGridWeightDelta (real,              intent (in)  :: semiAxisMajor,
!!                                 real,              intent (in)  :: semiAxisMinor,
!!                                 character (len=*), intent (in)  :: crossSectionFunctionType,
!!                                 real,              intent (in)  :: gaussianExponent,
!!                                 real,              intent (in)  :: gaussianRadiusMajor,
!!                                 real,              intent (in)  :: gaussianRadiusMinor,
!!                                 real,              intent (in)  :: gaussianCenterMajor,
!!                                 real,              intent (in)  :: gaussianCenterMinor,
!!                                 integer,           intent (in)  :: nTicsSemiAxisMajor,
!!                                 integer,           intent (in)  :: nTicsSemiAxisMinor,
!!                                 real,              intent (in)  :: deltaSemiAxisMajor,
!!                                 real,              intent (in)  :: deltaSemiAxisMinor,
!!                                 real,              intent (in)  :: firstTicSemiAxisMajor,
!!                                 real,              intent (in)  :: firstTicSemiAxisMinor,
!!                                 integer,           intent (in)  :: totalGridPoints,
!!                                 real,              intent (out) :: gridWeight)
!!
!! DESCRIPTION
!!
!!  Calculates the weight for the 3D beam delta grid, which is defined as the sum of the individual
!!  delta grid point weights. The delta grid must have been set up for the grid points to be ready
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
!!  nTicsSemiAxisMajor       : # of grid positions along the major semiaxis
!!  nTicsSemiAxisMinor       : # of grid positions along the minor semiaxis
!!  deltaSemiAxisMajor       : the tic separation for the elliptical major semiaxis of the 3D beam
!!  deltaSemiAxisMinor       : the tic separation for the elliptical minor semiaxis of the 3D beam
!!  firstTicSemiAxisMajor    : position of the 1st tic along the major semiaxis (in both directions)
!!  firstTicSemiAxisMinor    : position of the 1st tic along the minor semiaxis (in both directions)
!!  totalGridPoints          : the total number of grid points of the beam grid
!!  gridWeight               : the value of the grid weight
!!
!! NOTES
!!
!!***

subroutine ed_beam3DGridWeightDelta (semiAxisMajor,                    &
                                     semiAxisMinor,                    &
                                     crossSectionFunctionType,         &
                                     gaussianExponent,                 &
                                     gaussianRadiusMajor,              &
                                     gaussianRadiusMinor,              &
                                     gaussianCenterMajor,              &
                                     gaussianCenterMinor,              &
                                     nTicsSemiAxisMajor,               &
                                     nTicsSemiAxisMinor,               &
                                     deltaSemiAxisMajor,               &
                                     deltaSemiAxisMinor,               &
                                     firstTicSemiAxisMajor,            &
                                     firstTicSemiAxisMinor,            &
                                     totalGridPoints,                  &
                                                            gridWeight )

  implicit none

  real,              intent (in)  :: semiAxisMajor
  real,              intent (in)  :: semiAxisMinor
  character (len=*), intent (in)  :: crossSectionFunctionType
  real,              intent (in)  :: gaussianExponent
  real,              intent (in)  :: gaussianRadiusMajor
  real,              intent (in)  :: gaussianRadiusMinor
  real,              intent (in)  :: gaussianCenterMajor
  real,              intent (in)  :: gaussianCenterMinor
  integer,           intent (in)  :: nTicsSemiAxisMajor
  integer,           intent (in)  :: nTicsSemiAxisMinor
  real,              intent (in)  :: deltaSemiAxisMajor
  real,              intent (in)  :: deltaSemiAxisMinor
  real,              intent (in)  :: firstTicSemiAxisMajor
  real,              intent (in)  :: firstTicSemiAxisMinor
  integer,           intent (in)  :: totalGridPoints
  real,              intent (out) :: gridWeight

  gridWeight = 0.0

  return
end subroutine ed_beam3DGridWeightDelta
