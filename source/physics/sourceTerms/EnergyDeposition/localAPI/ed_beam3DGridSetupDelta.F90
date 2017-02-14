!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beam3DGridSetupDelta
!!
!! NAME
!!
!!  ed_beam3DGridSetupDelta
!!
!! SYNOPSIS
!!
!!  call ed_beam3DGridSetupDelta (real,    intent (in)    :: semiAxisMajor,
!!                                real,    intent (in)    :: semiAxisMinor,
!!                                real,    intent (in)    :: deltaSemiAxisMajor,
!!                                real,    intent (in)    :: deltaSemiAxisMinor,
!!                                integer, intent (inout) :: nGridPoints,
!!                                integer, intent (out)   :: nTicsSemiAxisMajor,
!!                                integer, intent (out)   :: nTicsSemiAxisMinor,
!!                                real,    intent (out)   :: firstTicSemiAxisMajor,
!!                                real,    intent (out)   :: firstTicSemiAxisMinor)
!!
!! DESCRIPTION
!!
!!  Sets up information about a rectangular cross sectional grid based on given tic
!!  separations along both major and minor semiaxes for a particular ellipsoidal
!!  3D beam. On output, the number of grid points overrides the input value, which
!!  is meaningless here for this routine.
!!
!! ARGUMENTS
!!
!!  semiAxisMajor         : the elliptical major semiaxis of the 3D beam
!!  semiAxisMinor         : the elliptical minor semiaxis of the 3D beam
!!  deltaSemiAxisMajor    : the tic separation for the elliptical major semiaxis of the 3D beam
!!  deltaSemiAxisMinor    : the tic separation for the elliptical minor semiaxis of the 3D beam
!!  nGridPoints           : on input  -> irrelevant, not needed here
!!                          on output -> the number of grid points determined (abort, if 0)
!!  nTicsSemiAxisMajor    : # of grid positions along the major semiaxis
!!  nTicsSemiAxisMinor    : # of grid positions along the minor semiaxis
!!  firstTicSemiAxisMajor : position of the 1st tic along the major semiaxis (in both directions)
!!  firstTicSemiAxisMinor : position of the 1st tic along the minor semiaxis (in both directions)
!!
!! NOTES
!!
!!  Currently the 1st tics start at 1/2 the tic separation along both axes in both
!!  directions (+ve and -ve axes).
!!
!!***

subroutine ed_beam3DGridSetupDelta (semiAxisMajor,                             &
                                    semiAxisMinor,                             &
                                    deltaSemiAxisMajor,                        &
                                    deltaSemiAxisMinor,                        &
                                                        nGridPoints,           &
                                                        nTicsSemiAxisMajor,    &
                                                        nTicsSemiAxisMinor,    &
                                                        firstTicSemiAxisMajor, &
                                                        firstTicSemiAxisMinor  )

  implicit none

  real,    intent (in)    :: semiAxisMajor
  real,    intent (in)    :: semiAxisMinor
  real,    intent (in)    :: deltaSemiAxisMajor
  real,    intent (in)    :: deltaSemiAxisMinor
  integer, intent (inout) :: nGridPoints
  integer, intent (out)   :: nTicsSemiAxisMajor
  integer, intent (out)   :: nTicsSemiAxisMinor
  real,    intent (out)   :: firstTicSemiAxisMajor
  real,    intent (out)   :: firstTicSemiAxisMinor

  nGridPoints           = 0
  nTicsSemiAxisMajor    = 0
  nTicsSemiAxisMinor    = 0
  firstTicSemiAxisMajor = 0.0
  firstTicSemiAxisMinor = 0.0

  return
end subroutine ed_beam3DGridSetupDelta
