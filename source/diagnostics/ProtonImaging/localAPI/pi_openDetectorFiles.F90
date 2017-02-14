!!****if* source/diagnostics/ProtonImaging/localAPI/pi_openDetectorFiles
!!
!! NAME
!!
!!  pi_openDetectorFiles
!!
!! SYNOPSIS
!!
!!  call pi_openDetectorFiles (real, intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Opens the detector files for recording screen protons for active proton beams only.
!!  To each detector there corresponds one file, characterized by the detector number and
!!  the current simulation time. Only the master processor opens the detector files. The
!!  names of the detector files is as follows:
!!
!!              <basenm> ProtonDetectorFile <detectorID> % <timeSimulation>
!!
!!  where <basenm> is the simulation base name, <detectorID> contains the detector number
!!  and <timeSimulation> is the current simulation time.
!!
!! ARGUMENTS
!!
!!  timeSimulation  : current simulation time
!!
!! NOTES
!!          
!!  Only detector files are opened for active beams. The detector number range is
!!  currently limited to 01-99. If this routine is called, then there is at least one
!!  active proton beam.
!!
!!***

subroutine pi_openDetectorFiles (timeSimulation)

  implicit none

  real, intent (in) :: timeSimulation

  return
end subroutine pi_openDetectorFiles
