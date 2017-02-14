!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beam3DGridWeightRecBeam
!!
!! NAME
!!
!!  ed_beam3DGridWeightRecBeam
!!
!! SYNOPSIS
!!
!!  call ed_beam3DGridWeightRecBeam (character (len=*), intent (in)  :: crossSectionFunctionType,
!!                                   real,              intent (in)  :: gaussianExponent,
!!                                   real,              intent (in)  :: gaussianRadiusMajor,
!!                                   real,              intent (in)  :: gaussianRadiusMinor,
!!                                   real,              intent (in)  :: gaussianCenterMajor,
!!                                   real,              intent (in)  :: gaussianCenterMinor,
!!                                   integer,           intent (in)  :: nTicsSemiAxisMajor,
!!                                   integer,           intent (in)  :: nTicsSemiAxisMinor,
!!                                   real,              intent (in)  :: deltaSemiAxisMajor,
!!                                   real,              intent (in)  :: deltaSemiAxisMinor,
!!                                   integer,           intent (in)  :: totalGridPoints,
!!                                   real,              intent (out) :: gridWeight)
!!
!! DESCRIPTION
!!
!!  Calculates the weight for the rectangular grid of a 3D rectangular beam, which is defined as
!!  the sum of the individual rectangular grid point weights. The rectangular grid must have been
!!  set up for the grid points to be ready for retrieval. If the number of grid points processed
!!  mismatches the total number of grid points of the beam grid, the program is aborted with a
!!  message.
!!
!! ARGUMENTS
!!
!!  crossSectionFunctionType : the cross section function type defining the grid point weigths
!!  gaussianExponent         : the gaussian exponent
!!  gaussianRadiusMajor      : the gaussian radius along the major semiaxis
!!  gaussianRadiusMinor      : the gaussian radius along the minor semiaxis
!!  gaussianCenterMajor      : the gaussian center location along the major semiaxis
!!  gaussianCenterMinor      : the gaussian center location along the minor semiaxis
!!  nTicsSemiAxisMajor       : # of grid positions along the major semiaxis
!!  nTicsSemiAxisMinor       : # of grid positions along the minor semiaxis
!!  deltaSemiAxisMajor       : the tic separation along the major semiaxis
!!  deltaSemiAxisMinor       : the tic separation along the minor semiaxis
!!  totalGridPoints          : the total number of grid points of the beam grid
!!  gridWeight               : the value of the grid weight
!!
!! NOTES
!!
!!***

subroutine ed_beam3DGridWeightRecBeam (crossSectionFunctionType,         &
                                       gaussianExponent,                 &
                                       gaussianRadiusMajor,              &
                                       gaussianRadiusMinor,              &
                                       gaussianCenterMajor,              &
                                       gaussianCenterMinor,              &
                                       nTicsSemiAxisMajor,               &
                                       nTicsSemiAxisMinor,               &
                                       deltaSemiAxisMajor,               &
                                       deltaSemiAxisMinor,               &
                                       totalGridPoints,                  &
                                                              gridWeight )

  implicit none

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
  integer,           intent (in)  :: totalGridPoints
  real,              intent (out) :: gridWeight

  gridWeight = 0.0

  return
end subroutine ed_beam3DGridWeightRecBeam
