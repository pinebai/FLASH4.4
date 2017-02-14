!!****if* source/diagnostics/ProtonImaging/localAPI/pi_setupDetectorFileNames
!!
!! NAME
!!
!!  pi_setupDetectorFileNames
!!
!! SYNOPSIS
!!
!!  call pi_setupDetectorFileNames (real, intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Sets the detector file names for further use during a proton imaging step. The
!!  names of the detector files are as follows:
!!
!!              <basenm> ProtonDetectorFile <detectorID> _ <timeSimulation>
!!
!!  where <basenm> is the simulation base name, <detectorID> contains the detector number
!!  and (optionally) <timeSimulation> is the current simulation time.
!!
!! ARGUMENTS
!!
!!  timeSimulation  : current simulation time
!!
!! NOTES
!!          
!!  1) As many detector file names are established as there were number of detectors
!!     specified. The detector number range is currently limited to 01-99. If more
!!     detectors were specified, the routine aborts.
!!
!!  2) A runtime parameter exists which enables suppression of the _<timeSimulation> part
!!     of the detector file names.
!!
!!  3) Although only the master processor will ever need these names, all processors
!!     get a copy of the detector names.
!!
!!***

subroutine pi_setupDetectorFileNames (timeSimulation)

  implicit none

  real, intent (in) :: timeSimulation

  return
end subroutine pi_setupDetectorFileNames
