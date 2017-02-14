!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/pi_transportBeamProtons
!!
!! NAME
!!
!!  pi_transportBeamProtons
!!
!! SYNOPSIS
!!
!!  call pi_transportBeamProtons (integer, intent (in) :: blockCount, 
!!                                integer, intent (in) :: blockList (:), 
!!                                real,    intent (in) :: timeStep,
!!                                real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Launches batches of beam protons onto the domain and transports them through
!!  the domain. Two things can happen to each proton: 1) the proton does not
!!  leave the domain during this time step or 2) the proton leaves the domain
!!  and is recorded on the screen. The domain structure and properties are assumed
!!  to remain frozen during the current time step transportation of the protons.
!!  Currently the protons do not interact with the domain, hence no domain update
!!  is necessary.
!!
!!  The code consists of the following basic steps:
!!
!!         1) Create the beam protons on the domain surface
!!         2) Follow the protons motion through the domain
!!         3) Record each proton leaving the domain on the detector screen(s)
!!         4) Store each proton remaining in the domain back to disk
!!         5) Write each screen proton to the detector file(s)
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
!!  The current implementation can handle more than one proton beam and more
!!  than one detector screen. The proton paths are calculated using classical
!!  Newton mechanics and are deflected due to average electrical and magnetic
!!  fields in each cell.
!!
!!***

subroutine pi_transportBeamProtons (blockCount, blockList, timeStep, timeSimulation)

  use ProtonImaging_data, ONLY : pi_beams,                      &
                                 pi_diskProtonCount,            &
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
                                 pi_timeResolvedProtonImaging,  &
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

  logical :: activeBeams
  logical :: activeProtons
  logical :: moreProtons
  logical :: moveProtons

  integer :: cycleNrBeams, cycleNrProtons
  integer :: error
!
!
!     ...Check, if there are active proton beams. A proton beam is considered
!        inactive, if the launching time is still less than the simulation time.
!        Once the simulation time gets passed the launching time, the beam is
!        activated and remains so until all its protons have been processed.
!        When the number of protons in a beam has dropped to 0, the beam is
!        also considered inactive.
!
!
  activeBeams = any (      pi_beams (:) % time2Launch     <= timeSimulation  &
                     .and. pi_beams (:) % numberOfProtons >  0               )

  if (activeBeams) then
      call pi_openDetectorFiles  (timeSimulation)
      if (pi_timeResolvedProtonImaging) then
          call pi_openDiskProtonFile ('old',saveOldRecords = .true.) ! to xfer new disk protons -> old file
          call pi_openDiskProtonFile ('new')                         ! scratch file to store new disk protons
      end if
  else
      return
  end if
!
!
!     ...Loop over the active beams.
!
!
  call Logfile_stamp ('Transporting Beam Protons ...','[Diagnostics]')
  call Timers_start  ("Transport Beam Protons")

  pi_timeStep = timeStep
  cycleNrBeams = 0

  do while (activeBeams)

     cycleNrBeams = cycleNrBeams + 1
     moveProtons  = .false.          ! because proton block ID's are valid at creation

     pi_protonCount = 0              !
     pi_diskProtonCount = 0          ! local processor counts
     pi_screenProtonCount = 0        !

     call Timers_start        ("Create Beam Protons")
     call pi_createProtons    (blockCount, blockList, timeSimulation)
     call pi_updateProtons    (moveProtons) 
     call pi_createProtonTags ()
     call Timers_stop         ("Create Beam Protons")

     if (pi_printProtons) then
         write (charCycleNr,'(I4.4)') cycleNrBeams
         fileLabel = "printProtonsBeamCyleNr" // charCycleNr // "ProcNr"
         call pi_printProtonsData (fileLabel, pi_globalMe)
     end if

     if (pi_IOplotProtons .and. pi_IOaddProtonsCapsule2Domain) then
         call pi_IOprotonsCapsule2Domain ()
         call Timers_start    ("IO Beam Protons -> Disk")
         call IO_writeProtons (pi_IOprotonCount,     &
                               pi_IOprotonTags,      &
                               pi_IOprotonPoints,    &
                               pi_IOprotonPointCount )
         call Timers_stop     ("IO Beam Protons -> Disk")
     end if
!
!
!     ...Transport the current batch of protons created through the domain and
!        either store them as disk protons or record them on the detector screens
!        as screen protons. Plot also the IO protons (if requested).
!
!
     moveProtons = .true.
     activeProtons = .true.
     cycleNrProtons = 0

     do while (activeProtons)

        cycleNrProtons = cycleNrProtons + 1

        call Timers_start     ("Trace Beam Protons")
        call pi_traceProtons  ()
        call pi_updateProtons (moveProtons)
        call Timers_stop      ("Trace Beam Protons")

        if (pi_IOplotProtons) then
            call Timers_start    ("IO Beam Protons -> Disk")
            call IO_writeProtons (pi_IOprotonCount,     &
                                  pi_IOprotonTags,      &
                                  pi_IOprotonPoints,    &
                                  pi_IOprotonPointCount )
            call Timers_stop     ("IO Beam Protons -> Disk")
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

     if (pi_timeResolvedProtonImaging) then
         call Timers_start ("Disk Protons -> Disk")
         call pi_flushDiskProtons2Disk ()
         call Timers_stop  ("Disk Protons -> Disk")
     end if
!
!
!     ...Check, if some proton beams are still active, in which case we have to
!        create the next set of protons, trace them through the domain and record
!        them on the detector screens.
!
!
     activeBeams = any (      pi_beams (:) % time2Launch     <= timeSimulation  &
                        .and. pi_beams (:) % numberOfProtons >  0               )

  end do
!
!
!     ...Transfer the new disk protons on the new disk proton file to the beginning
!        of the old disk proton file (only the master processor will do this).
!        Since beam protons can add disk protons onto already existing disk protons,
!        the old disk proton file will be appended.
!
!
  if (pi_timeResolvedProtonImaging) then
      call pi_xferDiskProtonsNew2OldFile (rewindOldFile = .false.)
  end if
!
!
!     ...If the user requested a plot of all detector screens, do it at this stage.
!        The screens are plotted by having a proton travel around the edges each
!        screen with a unique proton tag assignment.
!
!
  if (pi_IOplotProtons .and. pi_IOaddDetectorScreens) then
      call pi_IOdetectorScreens ()
      call Timers_start    ("IO Beam Protons -> Disk")
      call IO_writeProtons (pi_IOprotonCount,     &
                            pi_IOprotonTags,      &
                            pi_IOprotonPoints,    &
                            pi_IOprotonPointCount )
      call Timers_stop     ("IO Beam Protons -> Disk")
  end if
!
!
!     ...Do some final chores.
!
!
  call pi_closeDetectorFiles  ()
  if (pi_timeResolvedProtonImaging) then
      call pi_closeDiskProtonFile ('old')
      call pi_closeDiskProtonFile ('new')
  end if

  call Timers_stop   ("Transport Beam Protons")
  call Logfile_stamp ('Finished Transporting Beam Protons','[Diagnostics]')
!
!
!     ...Done for now!
!
!
  return
end subroutine pi_transportBeamProtons
