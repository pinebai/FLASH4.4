!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_printPulsesData
!!
!! NAME
!!
!!  ed_printPulsesData
!!
!! SYNOPSIS
!!
!!  call ed_printPulsesData ()
!!
!! DESCRIPTION
!!
!!  Utility routine, which prints detailed info regarding the generated data
!!  for all pulses to a text file. The information is written out to a file named
!!  <basenm>LaserPulsesPrint.txt, where <basenm> is the runtime parameter for
!!  output file names.
!!
!!***

subroutine ed_printPulsesData ()

  use EnergyDeposition_data,  ONLY : ed_baseName,              &
                                     ed_globalMe,              &
                                     ed_numberOfPulses,        &
                                     ed_pulseNumberOfSections, &
                                     ed_pulses

  implicit none
   
#include "Flash.h"
#include "constants.h"

  character (len = MAX_STRING_LENGTH) :: fileName

  integer :: fileUnit
  integer :: numberOfSections
  integer :: pulse
  integer :: section
  integer :: ut_getFreeFileUnit

  real    :: pulsePower
  real    :: pulseTime
!
!
!   ...Do the printout only on the master processor.
!
!
  if (ed_globalMe /= MASTER_PE) then
      return
  end if
!
!
!   ...Open the printout file.
!
!
  fileUnit = ut_getFreeFileUnit ()
  fileName = trim (ed_baseName) // "LaserPulsesPrint.txt"

  open (fileUnit, file = fileName)
!
!
!     ...Loop over all pulses.
!
!
  do pulse = 1, ed_numberOfPulses

     write (fileUnit,'(/)'   )
     write (fileUnit,'(a,i2)') "              LASER PULSE NR ",pulse
     write (fileUnit,'(/)'   )
     write (fileUnit,'(a)'   ) "         Time            Power (erg/s)"
     write (fileUnit,'(a)'   ) " ----------------------------------------"

     numberOfSections = ed_pulseNumberOfSections (pulse)

     do section = 1, numberOfSections

        pulseTime  = ed_pulses (section, pulse) % pulseTime
        pulsePower = ed_pulses (section, pulse) % pulsePower

        write (fileUnit,'(2es20.12)') pulseTime , pulsePower

     end do
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
end subroutine ed_printPulsesData
