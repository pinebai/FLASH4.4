!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserPulses/ed_setupPulses
!!
!! NAME
!!
!!  ed_setupPulses
!!
!! SYNOPSIS
!!
!!  call ed_setupPulses ()
!!
!! DESCRIPTION
!!
!!  Sets up and characterizes the pulses to be used in the simulation.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!  All needed info is read in as runtime parameters.
!!
!!***

subroutine ed_setupPulses ()

  use Driver_interface,            ONLY : Driver_abortFlash

  use EnergyDeposition_data,       ONLY : ed_Joule2erg,             &
                                          ed_numberOfPulses,        &
                                          ed_pulseNumberOfSections, &
                                          ed_pulses,                &
                                          ed_pulsesAreSetup

  use Logfile_interface,           ONLY : Logfile_stampMessage

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "EnergyDeposition.h"
#include "Flash.h"
#include "constants.h"

  integer :: numberOfSections
  integer :: pulse
  integer :: section

  real    :: oldTime
  real    :: pulsePower
  real    :: pulseTime

  character (len = MAX_STRING_LENGTH) :: parameterString

  if (ed_numberOfPulses < 1) then
      call Driver_abortFlash ("ed_setupPulses: No pulse defined!")
  end if

  if (ed_numberOfPulses > ED_MAXPULSES) then
      call Logfile_stampMessage ("ed_setupPulses: ERROR")
      call Logfile_stampMessage ("# of pulses > maximum # of pulse runtime parameters!")
      call Logfile_stampMessage ("Not enough pulse runtime parameters created.")
      call Logfile_stampMessage ("Increase the maximum number of pulse runtime parameters")
      call Logfile_stampMessage ("using the ed_maxPulses=<number of pulses> setup option.")
      call Driver_abortFlash    ("ed_setupPulses: Not enough pulse runtime parameters. See Log File!")
  end if
!
!
!     ...Allocate the pulses array(s) and read in the data.
!        Check for nonsense (negative time or power, nonincreasing times).
!
!  
  allocate (ed_pulseNumberOfSections                        (1:ed_numberOfPulses))
  allocate (ed_pulses                (1:ED_MAXPULSESECTIONS, 1:ed_numberOfPulses))
  
  do pulse = 1, ed_numberOfPulses

     write (parameterString,'(a,i0)') "ed_numberOfSections_", pulse
     call RuntimeParameters_get (parameterString, numberOfSections)

     if (numberOfSections < 2) then
         call Driver_abortFlash ("ed_setupPulses: No timeframe for current pulse!")
     end if

     if (numberOfSections > ED_MAXPULSESECTIONS) then
         call Logfile_stampMessage ("ed_setupPulses: ERROR")
         call Logfile_stampMessage ("There is not enough storage for the laser pulses sections.")
         call Logfile_stampMessage ("Increase the maximum number of pulses sections using the")
         call Logfile_stampMessage ("ed_maxPulseSections=<number of sections> setup option")
         call Driver_abortFlash    ("ed_setupPulses: Not enough storage for pulse sections! See Log File!")
     end if

     ed_pulseNumberOfSections (pulse) = numberOfSections

     write (parameterString,'(a,i0,a,i0)') "ed_time_", pulse, "_", 1
     call RuntimeParameters_get (parameterString, pulseTime)

     write (parameterString,'(a,i0,a,i0)') "ed_power_", pulse, "_", 1
     call RuntimeParameters_get (parameterString, pulsePower)

     if (pulsePower < 0.0) then
         call Driver_abortFlash ("ed_setupPulses: negative pulse power found!")
     end if

     ed_pulses (1,pulse) % pulseTime  = pulseTime
     ed_pulses (1,pulse) % pulsePower = pulsePower * ed_Joule2erg             ! Watt (J/s) -> erg/s

     oldTime = pulseTime

     do section = 2, numberOfSections

        write (parameterString,'(a,i0,a,i0)') "ed_time_", pulse, "_", section
        call RuntimeParameters_get (parameterString, pulseTime)
        
        write (parameterString,'(a,i0,a,i0)') "ed_power_", pulse, "_", section
        call RuntimeParameters_get (parameterString, pulsePower)

        if (pulseTime <= oldTime) then
            call Driver_abortFlash ("ed_setupPulses: Wrong time order for pulse sections!")
        end if

        if (pulsePower < 0.0) then
            call Driver_abortFlash ("ed_setupPulses: Negative pulse power specified!")
        end if

        ed_pulses (section, pulse) % pulseTime  = pulseTime
        ed_pulses (section, pulse) % pulsePower = pulsePower * ed_Joule2erg

        oldTime = pulseTime

     end do
        
  enddo
!
!
!     ...Set pulses status indicator.
!
!
  ed_pulsesAreSetup = .true.
!
!
!     ...Ready!
!
!
  return
end subroutine ed_setupPulses
