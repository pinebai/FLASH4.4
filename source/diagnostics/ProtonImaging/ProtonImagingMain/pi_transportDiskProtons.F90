!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/pi_transportDiskProtons
!!
!! NAME
!!
!!  pi_transportDiskProtons
!!
!! SYNOPSIS
!!
!!  call pi_transportDiskProtons (integer, intent (in) :: blockCount, 
!!                                integer, intent (in) :: blockList (:), 
!!                                real,    intent (in) :: timeStep,
!!                                real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Processes old batches of disk protons residing on disk from a previous time step.
!!  All protons that make it through the domain during the current time step are
!!  recorded on the detector screen(s). Any protons that stay in the domain are
!!  dumped back to disk as disk protons. Currently the protons do not interact with
!!  the domain, hence no domain update is necessary.
!!
!!  The code consists of the following basic steps:
!!
!!         1) Read the protons from disk back into the domain
!!         2) Follow the protons motion through the domain
!!         3) Record each proton leaving the domain on the detector screen(s)
!!         4) Store each proton remaining in the domain back to disk
!!         5) Write screen protons to the detector file(s)
!!
!!  The code can handle a larger number of protons than there is memory available
!!  by sending batches of protons and reusing available memory.
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
!!  The proton paths are calculated using classical Newton mechanics and are deflected
!!  due to average electrical and magnetic fields in each cell.
!!
!!***

subroutine pi_transportDiskProtons (blockCount, blockList, timeStep, timeSimulation)

  use ProtonImaging_data, ONLY : pi_diskProtonCount,            &
                                 pi_diskProtonOldFileName,      &
                                 pi_diskProtonOldFileNbatches,  &
                                 pi_globalComm,                 &
                                 pi_globalMe,                   &
                                 pi_IOaddDetectorScreens,       &
                                 pi_IOaddProtonsCapsule2Domain, &
                                 pi_IOplotProtons,              &
                                 pi_IOprotonCount,              &
                                 pi_IOprotonPointCount,         &
                                 pi_IOprotonPoints,             &
                                 pi_IOprotonTags,               &
                                 pi_printProtons,               &
                                 pi_protonCount,                &
                                 pi_protons,                    &
                                 pi_screenProtonCount,          &
                                 pi_timeStep

  use Logfile_interface,  ONLY : Logfile_stamp

  use Driver_interface,   ONLY : Driver_abortFlash

  use Timers_interface,   ONLY : Timers_start, &
                                 Timers_stop

  use IO_interface,       ONLY : IO_writeProtons

  use pi_interface,       ONLY : pi_closeDetectorFiles,         &
                                 pi_closeDiskProtonFile,        &
                                 pi_createProtons,              &
                                 pi_createProtonTags,           &
                                 pi_flushDiskProtons2Disk,      &
                                 pi_flushScreenProtons2Disk,    &
                                 pi_IOdetectorScreens,          &
                                 pi_IOprotonsCapsule2Domain,    &
                                 pi_openDetectorFiles,          &
                                 pi_openDiskProtonFile,         &
                                 pi_printProtonsData,           &
                                 pi_readDiskProtons,            &
                                 pi_traceProtons,               &
                                 pi_updateProtons,              &
                                 pi_xferDiskProtonsNew2OldFile

  implicit none

#include "constants.h"
#include "Flash.h"
#include "ProtonImaging.h"

  include "Flash_mpi.h"

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeStep
  real,    intent (in) :: timeSimulation

  character (len = 4 ) :: charCycleNr
  character (len = 32) :: fileLabel

  logical :: activeDisk
  logical :: activeProtons
  logical :: moreProtons
  logical :: moveProtons

  integer :: cycleNrDisk, cycleNrProtons
  integer :: error
!
!
!     ...Check, if there are disk protons on disk. If none are found, return.
!
!
  if (pi_diskProtonOldFileNbatches == 0) then
      return
  end if

  inquire (file = pi_diskProtonOldFileName, exist = activeDisk)

  if (activeDisk) then
      call pi_openDetectorFiles  (timeSimulation)
      call pi_openDiskProtonFile ('old',saveOldRecords = .false.)
      call pi_openDiskProtonFile ('new')
  else
      return
  end if
!
!
!     ...Loop over the active disk protons. Read the old disk protons from disk and
!        transport them through the domain for the current time step.
!
!
  call Logfile_stamp ('Transporting Disk Protons ...','[Diagnostics]')
  call Timers_start  ("Transport Disk Protons")

  pi_timeStep = timeStep
  cycleNrDisk = 0

  do while (activeDisk)

     cycleNrDisk = cycleNrDisk + 1
     moveProtons  = .false.        ! because proton block ID's are valid at read-in.

     pi_protonCount = 0
     pi_diskProtonCount = 0
     pi_screenProtonCount = 0

     call Timers_start        ("Read Disk Protons")
     call pi_readDiskProtons  (blockCount, blockList, activeDisk)
     call pi_updateProtons    (moveProtons) 
     call Timers_stop         ("Read Disk Protons")

     if (pi_printProtons) then
         write (charCycleNr,'(I4.4)') cycleNrDisk
         fileLabel = "printProtonsDiskCyleNr" // charCycleNr // "ProcNr"
         call pi_printProtonsData (fileLabel, pi_globalMe)
     end if
!
!
!     ...Transport the current batch of disk protons through the domain and record
!        them on the detector screens as screen protons, if they leave the domain.
!        Protons that stay in the domain are written out back to disk to a new disk
!        proton file. Plot also the IO protons (if requested).
!
!
     moveProtons = .true.
     activeProtons = .true.
     cycleNrProtons = 0

     do while (activeProtons)

        cycleNrProtons = cycleNrProtons + 1

        call Timers_start     ("Trace Disk Protons")
        call pi_traceProtons  ()
        call pi_updateProtons (moveProtons)
        call Timers_stop      ("Trace Disk Protons")

        if (pi_IOplotProtons) then
            call Timers_start    ("IO Disk Protons -> Disk")
            call IO_writeProtons (pi_IOprotonCount,     &
                                  pi_IOprotonTags,      &
                                  pi_IOprotonPoints,    &
                                  pi_IOprotonPointCount )
            call Timers_stop     ("IO Disk Protons -> Disk")
        end if

        moreProtons = (pi_protonCount > 0)

        if (moreProtons) then
            moreProtons = any (pi_protons (PROTON_BLCK,1:pi_protonCount) /= real (NONEXISTENT))
        end if

        call MPI_Allreduce (moreProtons,   &
                            activeProtons, &
                            1,             &
                            MPI_LOGICAL,   &
                            MPI_LOR,       &
                            pi_globalComm, &
                            error          )

     end do
!
!
!     ...Store the collected set of screen protons to disk by writing them out
!        to detector-specific files. Store the collected set of disk protons to
!        disk by writing (appending) them out to the new disk proton file.
!
!
     call Timers_start ("Screen Protons -> Disk")
     call pi_flushScreenProtons2Disk ()
     call Timers_stop  ("Screen Protons -> Disk")

     call Timers_start ("Disk Protons -> Disk")
     call pi_flushDiskProtons2Disk ()
     call Timers_stop  ("Disk Protons -> Disk")

  end do
!
!
!     ...Transfer the new disk protons on the new disk proton file to the beginning
!        of the old disk proton file (only the master processor will do this).
!        Since all old disk protons have been processed, the old disk proton file is
!        started from the beginning.
!
!
  call pi_xferDiskProtonsNew2OldFile (rewindOldFile = .true.)
!
!
!     ...If the user requested a plot of all detector screens, do it at this stage.
!        The screens are plotted by having a proton travel around the edges each
!        screen with a unique proton tag assignment.
!
!
  if (pi_IOplotProtons .and. pi_IOaddDetectorScreens) then
      call pi_IOdetectorScreens ()
      call Timers_start    ("IO Disk Protons -> Disk")
      call IO_writeProtons (pi_IOprotonCount,     &
                            pi_IOprotonTags,      &
                            pi_IOprotonPoints,    &
                            pi_IOprotonPointCount )
      call Timers_stop     ("IO Disk Protons -> Disk")
  end if
!
!
!     ...Do some final chores.
!
!
  call pi_closeDetectorFiles  ()
  call pi_closeDiskProtonFile ('old')
  call pi_closeDiskProtonFile ('new')

  call Timers_stop   ("Transport Disk Protons")
  call Logfile_stamp ('Finished Transporting Disk Protons','[Diagnostics]')
!
!
!     ...Done for now!
!
!
  return
end subroutine pi_transportDiskProtons
