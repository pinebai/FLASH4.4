!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beam3DGridSetupSquare
!!
!! NAME
!!
!!  ed_beam3DGridSetupSquare
!!
!! SYNOPSIS
!!
!!  call ed_beam3DGridSetupSquare (real,    intent (in)    :: semiAxisMajor,
!!                                 real,    intent (in)    :: semiAxisMinor,
!!                                 integer, intent (inout) :: nGridPoints,
!!                                 integer, intent (out)   :: nTicsSemiAxisMajor,
!!                                 integer, intent (out)   :: nTicsSemiAxisMinor,
!!                                 real,    intent (out)   :: delta)
!!
!! DESCRIPTION
!!
!!  Sets up information about the square cross sectional grid for a particular ellipsoidal
!!  3D beam. On input, the number of grid points is the wanted number of grid points.
!!  On output, the number of grid points might have changed and it overrides the input value.
!!
!! ARGUMENTS
!!
!!  semiAxisMajor      : the elliptical major semiaxis of the 3D beam
!!  semiAxisMinor      : the elliptical minor semiaxis of the 3D beam
!!  nGridPoints        : on input  -> the number of grid points wanted
!!                       on output -> the optimum number of grid points close to the input value
!!  nTicsSemiAxisMajor : # of grid positions along the major semiaxis
!!  nTicsSemiAxisMinor : # of grid positions along the minor semiaxis
!!  delta              : the separation distance between two consecutive tics
!!
!! NOTES
!!
!!***

subroutine ed_beam3DGridSetupSquare (semiAxisMajor,                     &
                                     semiAxisMinor,                     &
                                                    nGridPoints,        &
                                                    nTicsSemiAxisMajor, &
                                                    nTicsSemiAxisMinor, &
                                                    delta               )

  use ed_interface,  ONLY : ed_gridEllipseSquare

  implicit none

  real,    intent (in)    :: semiAxisMajor
  real,    intent (in)    :: semiAxisMinor
  integer, intent (inout) :: nGridPoints
  integer, intent (out)   :: nTicsSemiAxisMajor
  integer, intent (out)   :: nTicsSemiAxisMinor
  real,    intent (out)   :: delta

  real :: aspectRatio
  real :: deltaNormalized
!
!
!     ...Set up the square grid inside the elliptical beam cross section area. The special
!        situation of having only 1 grid point at the center is handled separately.
!        This 'no grid' situation is characterized by having no grid tics at all.
!
!
  if (nGridPoints == 1) then

      nTicsSemiAxisMajor = 0
      nTicsSemiAxisMinor = 0

  else

      aspectRatio = semiAxisMajor / semiAxisMinor

      call ed_gridEllipseSquare (aspectRatio,                     &
                                              nGridPoints,        &
                                              nTicsSemiAxisMajor, &
                                              nTicsSemiAxisMinor, &
                                              deltaNormalized     )

      delta = deltaNormalized * semiAxisMajor

  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beam3DGridSetupSquare
