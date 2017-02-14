!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beam3DGridWeightSquare
!!
!! NAME
!!
!!  ed_beam3DGridWeightSquare
!!
!! SYNOPSIS
!!
!!  call ed_beam3DGridWeightSquare (real,              intent (in)  :: semiAxisMajor,
!!                                  real,              intent (in)  :: semiAxisMinor,
!!                                  character (len=*), intent (in)  :: crossSectionFunctionType,
!!                                  real,              intent (in)  :: gaussianExponent,
!!                                  real,              intent (in)  :: gaussianRadiusMajor,
!!                                  real,              intent (in)  :: gaussianRadiusMinor,
!!                                  real,              intent (in)  :: gaussianCenterMajor,
!!                                  real,              intent (in)  :: gaussianCenterMinor,
!!                                  integer,           intent (in)  :: nTicsSemiAxisMajor,
!!                                  integer,           intent (in)  :: nTicsSemiAxisMinor,
!!                                  real,              intent (in)  :: delta,
!!                                  integer,           intent (in)  :: totalGridPoints,
!!                                  real,              intent (out) :: gridWeight)
!!
!! DESCRIPTION
!!
!!  Calculates the weight for the 3D beam square grid, which is defined as the sum of the individual
!!  square grid point weights. The square grid must have been set up for the grid points to be ready
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
!!  delta                    : the separation distance between two consecutive tics
!!  totalGridPoints          : the total number of grid points of the beam grid
!!  gridWeight               : the value of the grid weight
!!
!! NOTES
!!
!!***

subroutine ed_beam3DGridWeightSquare (semiAxisMajor,                    &
                                      semiAxisMinor,                    &
                                      crossSectionFunctionType,         &
                                      gaussianExponent,                 &
                                      gaussianRadiusMajor,              &
                                      gaussianRadiusMinor,              &
                                      gaussianCenterMajor,              &
                                      gaussianCenterMinor,              &
                                      nTicsSemiAxisMajor,               &
                                      nTicsSemiAxisMinor,               &
                                      delta,                            &
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
  real,              intent (in)  :: delta
  integer,           intent (in)  :: totalGridPoints
  real,              intent (out) :: gridWeight

  gridWeight = 0.0

  return
end subroutine ed_beam3DGridWeightSquare
