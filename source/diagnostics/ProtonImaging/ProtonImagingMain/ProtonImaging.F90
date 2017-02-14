!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonImaging
!!
!! NAME
!!
!!  ProtonImaging
!!
!! SYNOPSIS
!!
!!  call ProtonImaging (integer, intent (in) :: blockCount, 
!!                      integer, intent (in) :: blockList (:), 
!!                      real,    intent (in) :: timeStep,
!!                      real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Main routine controlling the actions to be performed during the current
!!  time step:
!!
!!         1) transport old disk protons throught the domain
!!         2) create and transport beam protons through the domain
!!
!!  The routine calls two subroutines, each handling the above 2 steps.
!!  It is important that first old remaining disk protons from previous time
!!  steps must be processed before considering new beam protons, hence the
!!  order of these two routines matter.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  timeStep       : The current time step duration
!!  timeSimulation : current simulation time
!!
!! NOTES
!!
!!  If no time resolved proton imaging is performed, the first step handling
!!  disk protons is omitted.
!!          
!!***

subroutine ProtonImaging (blockCount, blockList, timeStep, timeSimulation)

  use ProtonImaging_data, ONLY : pi_currentStepNumber,         &
                                 pi_globalMe,                  &
                                 pi_IOplotProtons,             &
                                 pi_monitorFileUnit,           &
                                 pi_timeResolvedProtonImaging, &
                                 pi_useIOprotonPlot,           &
                                 pi_useProtonImaging

  use Logfile_interface,  ONLY : Logfile_stamp

  use Driver_interface,   ONLY : Driver_abortFlash, &
                                 Driver_getNStep

  use Timers_interface,   ONLY : Timers_start, &
                                 Timers_stop

  use IO_interface,       ONLY : IO_checkForPlot,     &
                                 IO_endProtonWrite,   &
                                 IO_startProtonWrite

  use pi_interface,       ONLY : pi_closeMonitorFile,       &
                                 pi_openMonitorFile,        &
                                 pi_setupDetectorFileNames, &
                                 pi_transportBeamProtons,   &
                                 pi_transportDiskProtons

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeStep
  real,    intent (in) :: timeSimulation
!
!
!     ...Check, if proton imaging is needed at all. If not, return at once.
!
!
  if (.not. pi_useProtonImaging) then
       return
  end if
!
!
!     ...Get the current step number and set up the detector file names to be
!        used during the current time step. Open the monitor file and record
!        the step number into the monitor file.
!
!
  call Driver_getNStep (pi_currentStepNumber)
  call pi_setupDetectorFileNames (timeSimulation)
  call pi_openMonitorFile ()

  if (pi_globalMe == MASTER_PE) then
      write (pi_monitorFileUnit,'(/)')
      write (pi_monitorFileUnit,'(a,i5)') ' CURRENT STEP NUMBER = ',pi_currentStepNumber
      write (pi_monitorFileUnit,'(/)')
  end if
!
!
!     ...Do the logistics for plotting IO protons (if requested) to the plotfile
!        for the current time step.
!
!
  pi_IOplotProtons = .false.

  if (pi_useIOprotonPlot) then
      call IO_checkForPlot (pi_IOplotProtons)
  end if

  if (pi_IOplotProtons) then
      call IO_startProtonWrite ()
  end if
!
!
!     ...Call the 2 transportation routines.
!
!
  call Logfile_stamp ('Starting Proton Imaging ...','[Diagnostics]')
  call Timers_start  ("ProtonImaging")

  if (pi_timeResolvedProtonImaging) then
      call pi_transportDiskProtons (blockCount, blockList, timeStep, timeSimulation)
  end if

  call pi_transportBeamProtons (blockCount, blockList, timeStep, timeSimulation)

  call Timers_stop   ("ProtonImaging")
  call Logfile_stamp ('Finished Proton Imaging','[Diagnostics]')
!
!
!     ...End IO proton plotting (if requested).
!
!
  if (pi_IOplotProtons) then
      call IO_endProtonWrite ()
  end if
!
!
!     ...Close the monitor file.
!
!
  call pi_closeMonitorFile ()
!
!
!     ...Done for now!
!
!
  return
end subroutine ProtonImaging
