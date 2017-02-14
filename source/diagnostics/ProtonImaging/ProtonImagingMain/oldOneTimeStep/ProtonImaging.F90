!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonImaging
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
!!  Launches batches of protons onto the domain and records them on detector
!!  screen(s). This is the driver routine for proton imaging at a particular
!!  time during the simulation. The domain structure and properties are assumed
!!  to remain frozen during the traversion of the protons.
!!
!!  Currently the protons do not interact with the domain, hence no domain update
!!  is necessary.
!!
!!  The code consists of the following basic steps:
!!
!!         1) Create the protons on the domain surface
!!         2) Follow the protons motion through the domain
!!         3) Record each proton on the detector screen(s)
!!         4) Write each screen proton to the detector file(s)
!!
!!  The code can handle a larger number of protons than there is memory available
!!  by sending batches of protons and reusing available memory.
!!
!! ARGUMENTS
!!
!!  blockCount      : Number of blocks on current processor
!!  blockList       : All block ID numbers
!!  timeStep        : The current time step duration
!!  timeSimulation  : current simulation time
!!
!! NOTES
!!          
!!  The current implementation can handle more than one proton beam and more
!!  than one detector screen. The proton paths are calculated using classical
!!  Newton mechanics and are deflected due to average electrical and magnetic
!!  fields in each cell.
!!
!!***

subroutine ProtonImaging (blockCount, blockList, timeStep, timeSimulation)

  use ProtonImaging_data, ONLY : pi_beams,                      &
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
                                 pi_tagMax,                     &
                                 pi_useIOprotonPlot,            &
                                 pi_useProtonImaging

  use Logfile_interface,  ONLY : Logfile_stamp

  use Driver_interface,   ONLY : Driver_abortFlash

  use Timers_interface,   ONLY : Timers_start, &
                                 Timers_stop

  use IO_interface,       ONLY : IO_checkForPlot,     &
                                 IO_endProtonWrite,   &
                                 IO_startProtonWrite, &
                                 IO_writeProtons

  use pi_interface,       ONLY : pi_closeDetectorFiles,         &
                                 pi_createProtons,              &
                                 pi_createProtonTags,           &
                                 pi_flushScreenProtons2Disk,    &
                                 pi_globalMe,                   &
                                 pi_IOdetectorScreens,          &
                                 pi_IOprotonsCapsule2Domain,    &
                                 pi_openDetectorFiles,          &
                                 pi_printProtonsData,           &
                                 pi_traceProtons,               &
                                 pi_updateProtons

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
!     ...Check, if proton imaging is needed at all. If not, return at once.
!
!
  if (.not. pi_useProtonImaging) then
       return
  end if

  call Logfile_stamp ('Starting Proton Imaging ...','[Diagnostics]')
  call Timers_start  ("ProtonImaging")
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
      call pi_openDetectorFiles (timeSimulation)
  end if
!
!
!     ...Do the logistics for plotting IO protons (if requested).
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
!     ...Loop over the active beams.
!
!
  pi_tagMax = 0
  cycleNrBeams = 0

  do while (activeBeams)

     cycleNrBeams = cycleNrBeams + 1

     moveProtons  = .false.  ! because proton block ID's are valid at creation

     call Timers_start        ("Create Protons")
     call pi_createProtons    (blockCount, blockList, timeSimulation)
     call pi_updateProtons    (moveProtons) 
     call pi_createProtonTags ()
     call Timers_stop         ("Create Protons")

     if (pi_printProtons) then
         write (charCycleNr,'(I4.4)') cycleNrBeams
         fileLabel = "printProtonsBeamCyleNr" // charCycleNr // "ProcNr"
         call pi_printProtonsData (fileLabel, pi_globalMe)
     end if

     if (pi_IOplotProtons .and. pi_IOaddProtonsCapsule2Domain) then
         call pi_IOprotonsCapsule2Domain ()
         call Timers_start    ("IO Protons -> Disk")
         call IO_writeProtons (pi_IOprotonCount,     &
                               pi_IOprotonTags,      &
                               pi_IOprotonPoints,    &
                               pi_IOprotonPointCount )
         call Timers_stop     ("IO Protons -> Disk")
     end if
!
!
!     ...Transport the current batch of protons created through the domain and
!        record them on the detector screens as screen protons. Plot also the
!        IO protons (if requested).
!
!
     moveProtons   = .true.
     activeProtons = .true.

     cycleNrProtons = 0

     do while (activeProtons)

        cycleNrProtons = cycleNrProtons + 1

        call Timers_start     ("Transport Protons")
        call pi_traceProtons  ()
        call pi_updateProtons (moveProtons)
        call Timers_stop      ("Transport Protons")

        if (pi_IOplotProtons) then
            call Timers_start    ("IO Protons -> Disk")
            call IO_writeProtons (pi_IOprotonCount,     &
                                  pi_IOprotonTags,      &
                                  pi_IOprotonPoints,    &
                                  pi_IOprotonPointCount )
            call Timers_stop     ("IO Protons -> Disk")
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
!        to detector-specific files.
!
!
     call Timers_start ("Screen Protons -> Disk")
     call pi_flushScreenProtons2Disk ()
     call Timers_stop  ("Screen Protons -> Disk")
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
!     ...If the user requested a plot of all detector screens, do it at this stage.
!        The screens are plotted by having a proton travel around the edges each
!        screen with a unique proton tag assignment.
!
!
  if (pi_IOplotProtons .and. pi_IOaddDetectorScreens) then
      call pi_IOdetectorScreens ()
      call Timers_start    ("IO Protons -> Disk")
      call IO_writeProtons (pi_IOprotonCount,     &
                            pi_IOprotonTags,      &
                            pi_IOprotonPoints,    &
                            pi_IOprotonPointCount )
      call Timers_stop     ("IO Protons -> Disk")
  end if
!
!
!     ...Do some final chores.
!
!
  if (pi_IOplotProtons) then
      call IO_endProtonWrite ()
  end if

  call pi_closeDetectorFiles ()

  call Timers_stop   ("ProtonImaging")
  call Logfile_stamp ('Finished Proton Imaging','[Diagnostics]')
!
!
!     ...Done for now!
!
!
  return
end subroutine ProtonImaging
