!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beamsInfo1DRec
!!
!! NAME
!!
!!  ed_beamsInfo1DRec
!!
!! SYNOPSIS
!!
!!  call ed_beamsInfo1DRec ()
!!
!! DESCRIPTION
!!
!!  Generates information about the beams for those geometries consisting formally of 1D
!!  rectangular grids (cartesian + spherical). In here all beam information is generated
!!  that can be generated at initialization. Currently it contains the following:
!!
!!         1) Calculate the lens to target distance and incorporate some
!!            additional info (time duration) regarding the pulses associated
!!            with each beam.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsInfo1DRec ()

  use Driver_interface,       ONLY : Driver_abortFlash

  use Logfile_interface,      ONLY : Logfile_stamp

  use EnergyDeposition_data,  ONLY : ed_beams,                 &
                                     ed_numberOfBeams,         &
                                     ed_pulseNumberOfSections, &
                                     ed_pulses

  implicit none

#include "EnergyDeposition.h"
#include "Flash.h"
#include "constants.h"

  integer :: beam
  integer :: lastSection
  integer :: pulseNumber

  real    :: beamLength
  real    :: beginTime, endTime
  real    :: lensX
  real    :: targetX
!
!
!     ...Loop over all beams.
!
!
  do beam = 1, ed_numberOfBeams

     lensX       = ed_beams (beam) % lensX
     pulseNumber = ed_beams (beam) % pulseNumber
     targetX     = ed_beams (beam) % targetX
!
!
!     ...Set the dimensionality indicator of the beam.
!
!
     ed_beams (beam) % dimensionality = 1
!
!
!     ...Calculate distance from lens to target and pulse time bounds for each beam.
!
!
     beamLength = abs (targetX - lensX)

     lastSection = ed_pulseNumberOfSections (pulseNumber)
     beginTime   = ed_pulses (1,             pulseNumber) % pulseTime
     endTime     = ed_pulses (lastSection,   pulseNumber) % pulseTime

     ed_beams (beam) % distanceLens2Target = beamLength
     ed_beams (beam) % pulseStartingTime   = beginTime
     ed_beams (beam) % pulseEndingTime     = endTime

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beamsInfo1DRec
