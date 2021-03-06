!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beam3DGridSetupRecBeam
!!
!! NAME
!!
!!  ed_beam3DGridSetupRecBeam
!!
!! SYNOPSIS
!!
!!  call ed_beam3DGridSetupRecBeam (real,    intent (in)    :: semiAxisMajor,
!!                                  real,    intent (in)    :: semiAxisMinor,
!!                                  integer, intent (inout) :: nGridPoints,
!!                                  integer, intent (inout) :: nTicsSemiAxisMajor,
!!                                  integer, intent (inout) :: nTicsSemiAxisMinor,
!!                                  real,    intent (out)   :: deltaSemiAxisMajor,
!!                                  real,    intent (out)   :: deltaSemiAxisMinor)
!!
!! DESCRIPTION
!!
!!  Sets up information about the rectangular beam cross sectional grid for a particular 3D beam.
!!  The rectagular beam is defined as having a rectangular crossection:
!!
!!                            -----------------------------
!!                           |              |              |
!!                           |            b |              |
!!                           |              |              |   a = major semi axis
!!                            -----------------------------
!!                           |              |      a       |   b = minor semi axis
!!                           |              |              |
!!                           |              |              |
!!                            -----------------------------
!!
!!  On input, the number of grid points is the wanted number of grid points. On output, the
!!  number of grid points might have changed and it overrides the input value.
!!
!!  There is the option of passing predefined values for the number of tics along the major
!!  and minor semi axes. The following actions are being taken, depending on the input values
!!  of these two grid parameters:
!!
!!            nTicsSemiAxisMajor =< 0   ->   automatic evaluation to rectangular tiling
!!            nTicsSemiAxisMinor =< 0   ->   automatic evaluation to rectangular tiling
!!
!!            nTicsSemiAxisMajor  > 0   ->   take this value
!!            nTicsSemiAxisMinor =< 0   ->   adjust to number of grid points wanted
!!
!!            nTicsSemiAxisMajor =< 0   ->   adjust to number of grid points wanted
!!            nTicsSemiAxisMinor  > 0   ->   take this value
!!
!!            nTicsSemiAxisMajor  > 0   ->   take this value
!!            nTicsSemiAxisMinor  > 0   ->   take this value
!!
!!
!! ARGUMENTS
!!
!!  semiAxisMajor      : the rectangular major semiaxis of the 3D beam
!!  semiAxisMinor      : the rectangular minor semiaxis of the 3D beam
!!  nGridPoints        : on input  -> the number of grid points wanted
!!                       on output -> the optimum number of grid points close to the input value
!!  nTicsSemiAxisMajor : # of grid positions along the major semi axis (if > 0 -> assumed user given)
!!  nTicsSemiAxisMinor : # of grid positions along the minor semi axis (if > 0 -> assumed user given)
!!  deltaSemiAxisMajor : the tic spacing along the major semi axis
!!  deltaSemiAxisMinor : the tic spacing along the minor semi axis
!!
!! NOTES
!!
!!  1) The automatic rectangular tiling is performed in such a way as to preserve the closest
!!     way possible the rectangular structure of the beam crossection. The shape of each rectangular
!!     tile inside the beam would correspond roughly the shape of the beam crossection.
!!
!!  2) The outer grid points will lay on the outer surface of the rectangular beam, i.e. the
!!     grid will cover the entire rectangular crossection.
!!
!!***

subroutine ed_beam3DGridSetupRecBeam (semiAxisMajor,      &
                                      semiAxisMinor,      &
                                      nGridPoints,        &
                                      nTicsSemiAxisMajor, &
                                      nTicsSemiAxisMinor, &
                                      deltaSemiAxisMajor, &
                                      deltaSemiAxisMinor  )

  implicit none

  real,    intent (in)    :: semiAxisMajor
  real,    intent (in)    :: semiAxisMinor
  integer, intent (inout) :: nGridPoints
  integer, intent (inout) :: nTicsSemiAxisMajor
  integer, intent (inout) :: nTicsSemiAxisMinor
  real,    intent (out)   :: deltaSemiAxisMajor
  real,    intent (out)   :: deltaSemiAxisMinor

  nGridPoints        = 0
  nTicsSemiAxisMajor = 0
  nTicsSemiAxisMinor = 0
  deltaSemiAxisMajor = 0.0
  deltaSemiAxisMinor = 0.0

  return
end subroutine ed_beam3DGridSetupRecBeam
