!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonUtilities/pi_printBeamsData
!!
!! NAME
!!
!!  pi_printBeamsData
!!
!! SYNOPSIS
!!
!!  call pi_printBeamsData ()
!!
!! DESCRIPTION
!!
!!  Utility routine, which prints detailed info regarding the generated data
!!  for all beams to a text file. The information is written out to a file named
!!  <basenm>ProtonBeamsPrint.txt, where <basenm> is the runtime parameter for
!!  output file names.
!!
!!***

subroutine pi_printBeamsData ()

  use Driver_interface,    ONLY : Driver_abortFlash

  use ProtonImaging_data,  ONLY : pi_baseName,      &
                                  pi_beams,         &
                                  pi_globalMe,      &
                                  pi_numberOfBeams, &
                                  pi_speedOfLight

  implicit none
   
#include "constants.h"
#include "ProtonImaging.h"
#include "Flash.h"

  character (len =  MAX_STRING_LENGTH) :: fileName

  logical :: ignoreBnd

  integer :: beam
  integer :: capsuleGL, capsuleGN
  integer :: detector
  integer :: dimensions
  integer :: fileUnit
  integer :: nProtons, nProtonsGrain
  integer :: seed
  integer :: ut_getFreeFileUnit

  real    :: apertureAngle
  real    :: axisUnit1X, axisUnit1Y, axisUnit1Z
  real    :: axisUnit2X, axisUnit2Y, axisUnit2Z
  real    :: axisUnit3X, axisUnit3Y, axisUnit3Z
  real    :: c
  real    :: capsuleGS
  real    :: capsuleRadius
  real    :: capsuleX, capsuleY, capsuleZ
  real    :: distC2T
  real    :: protonEnergy, protonSpeed
  real    :: targetRadius
  real    :: targetX, targetY, targetZ
  real    :: time2Launch
!
!
!   ...Do the printout only on the master processor.
!
!
  if (pi_globalMe /= MASTER_PE) then
      return
  end if
!
!
!   ...Open the printout file.
!
!
  fileUnit = ut_getFreeFileUnit ()
  fileName = trim (pi_baseName) // "ProtonBeamsPrint.txt"

  open (fileUnit, file = fileName)
!
!
!     ...Loop over all beams.
!
!
  c = pi_speedOfLight

  do beam = 1, pi_numberOfBeams

     apertureAngle = pi_beams (beam) % apertureAngle
     axisUnit1X    = pi_beams (beam) % axisUnit1X
     axisUnit1Y    = pi_beams (beam) % axisUnit1Y
     axisUnit1Z    = pi_beams (beam) % axisUnit1Z
     axisUnit2X    = pi_beams (beam) % axisUnit2X
     axisUnit2Y    = pi_beams (beam) % axisUnit2Y
     axisUnit2Z    = pi_beams (beam) % axisUnit2Z
     axisUnit3X    = pi_beams (beam) % axisUnit3X
     axisUnit3Y    = pi_beams (beam) % axisUnit3Y
     axisUnit3Z    = pi_beams (beam) % axisUnit3Z
     capsuleGL     = pi_beams (beam) % capsuleGrainLevel
     capsuleGN     = pi_beams (beam) % capsuleNumberOfGrains
     capsuleGS     = pi_beams (beam) % capsuleGrainSize
     capsuleRadius = pi_beams (beam) % capsuleRadius
     capsuleX      = pi_beams (beam) % capsuleX
     capsuleY      = pi_beams (beam) % capsuleY
     capsuleZ      = pi_beams (beam) % capsuleZ
     detector      = pi_beams (beam) % detector
     dimensions    = pi_beams (beam) % dimensionality
     distC2T       = pi_beams (beam) % distanceCapsule2Target
     protonSpeed   = pi_beams (beam) % initialProtonSpeed
     ignoreBnd     = pi_beams (beam) % noBoundaryCondition
     nProtons      = pi_beams (beam) % numberOfProtons
     nProtonsGrain = pi_beams (beam) % numberOfProtonsPerGrain
     protonEnergy  = pi_beams (beam) % protonEnergy
     seed          = pi_beams (beam) % randomNumberSeed
     targetRadius  = pi_beams (beam) % targetRadius
     targetX       = pi_beams (beam) % targetX
     targetY       = pi_beams (beam) % targetY
     targetZ       = pi_beams (beam) % targetZ
     time2Launch   = pi_beams (beam) % time2Launch

     write (fileUnit,'(/)'        )
     write (fileUnit,'(a,i2)'     ) "               PROTON BEAM NR ",beam
     write (fileUnit,'(/)'        )
     write (fileUnit,'(a,es20.12)') " Beam aperture angle (rad)                      = ", apertureAngle
     write (fileUnit,'(a,es20.12)') " Distance capsule center --> target center      = ", distC2T
     write (fileUnit,'(a,es20.12)') " Proton energy (in MeV)                         = ", protonEnergy
     write (fileUnit,'(a,es20.12)') " Initial proton speed (in units of cm/s)        = ", protonSpeed
     write (fileUnit,'(a,es20.12)') " Initial proton speed (in units of light speed) = ", protonSpeed / c
     write (fileUnit,'(a,es20.12)') " Beam launches protons, if simulation time is  >= ", time2Launch
     write (fileUnit,'(a,L1)'     ) " The domain boundary conditions are ignored     = ", ignoreBnd
     write (fileUnit,'(a,i20)'    ) " Number of protons in beam                      = ", nProtons
     write (fileUnit,'(a,i20)'    ) " Beam protons are recorded by detector number   = ", detector

     select case (dimensions)

     case (3)

     write (fileUnit,'(a,i20)'    ) " Random number generator seed                   = ", seed
     write (fileUnit,'(a,i20)'    ) " Capsule grain level                            = ", capsuleGL
     write (fileUnit,'(a,i20)'    ) " Capsule number of grains                       = ", capsuleGN
     write (fileUnit,'(a,i20)'    ) " Capsule number of protons per grain            = ", nProtonsGrain
     write (fileUnit,'(a,es20.12)') " Capsule radius                                 = ", capsuleRadius
     write (fileUnit,'(a,es20.12)') " Capsule grain size (cube side length)          = ", capsuleGS
     write (fileUnit,'(a,es20.12)') " Capsule center global x-coordinate             = ", capsuleX
     write (fileUnit,'(a,es20.12)') " Capsule center global y-coordinate             = ", capsuleY
     write (fileUnit,'(a,es20.12)') " Capsule center global z-coordinate             = ", capsuleZ
     write (fileUnit,'(a,es20.12)') " 1st axis unit vector x-coordinate              = ", axisUnit1X
     write (fileUnit,'(a,es20.12)') " 1st axis unit vector y-coordinate              = ", axisUnit1Y
     write (fileUnit,'(a,es20.12)') " 1st axis unit vector z-coordinate              = ", axisUnit1Z
     write (fileUnit,'(a,es20.12)') " 2nd axis unit vector x-coordinate              = ", axisUnit2X
     write (fileUnit,'(a,es20.12)') " 2nd axis unit vector y-coordinate              = ", axisUnit2Y
     write (fileUnit,'(a,es20.12)') " 2nd axis unit vector z-coordinate              = ", axisUnit2Z
     write (fileUnit,'(a,es20.12)') " 3rd axis unit vector x-coordinate              = ", axisUnit3X
     write (fileUnit,'(a,es20.12)') " 3rd axis unit vector y-coordinate              = ", axisUnit3Y
     write (fileUnit,'(a,es20.12)') " 3rd axis unit vector z-coordinate              = ", axisUnit3Z
     write (fileUnit,'(a,es20.12)') " Target radius length                           = ", targetRadius
     write (fileUnit,'(a,es20.12)') " Target global x-coordinate                     = ", targetX
     write (fileUnit,'(a,es20.12)') " Target global y-coordinate                     = ", targetY
     write (fileUnit,'(a,es20.12)') " Target global z-coordinate                     = ", targetZ
     write (fileUnit,'(a)'        ) " ------------------------------------------------ "

     case default

     call Driver_abortFlash ("pi_printBeamsData: Bad dimensionality of beam (programmers fault)!")

     end select

  end do
!
!
!   ...Close the printout file.
!
!
  close (fileUnit)
!
!
!    ...Ready!
!
!  
  return
end subroutine pi_printBeamsData
