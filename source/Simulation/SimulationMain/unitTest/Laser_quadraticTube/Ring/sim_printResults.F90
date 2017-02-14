!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/Ring/sim_printResults
!!
!!  NAME 
!!
!!   sim_printResults
!!
!!  SYNOPSIS
!!
!!   sim_printResults ()
!!
!!  DESCRIPTION
!!
!!   This routine prints some general info about the simulation. It prints the obtained maximum
!!   square radial deviation and the maximum absolute power deviation to a .dat file.
!!   Only the master processor does the printout.
!!
!!***

subroutine sim_printResults ()

  use Simulation_data,    ONLY : sim_A,                  &
                                 sim_baseName,           &
                                 sim_beamTargetRadius,   &
                                 sim_focalPointRadius,   &
                                 sim_globalMe,           &
                                 sim_maxDeltaP,          &
                                 sim_maxDeltaRsqr,       &
                                 sim_nc,                 &
                                 sim_numRaysAnalyzed,    &
                                 sim_numRaysLaunched,    &
                                 sim_nuw,                &
                                 sim_nw,                 &
                                 sim_powerDecayFactor,   &
                                 sim_refinementLevel,    &
                                 sim_Tw,                 &
                                 sim_xc,                 &
                                 sim_xw,                 &
                                 sim_yfocal,             &
                                 sim_Z,                  &
                                 sim_zc,                 &
                                 sim_zw,                 &
                                 sim_totalNumberOfBoxes, &
                                 sim_boxPower,           &
                                 sim_boxRadius,          &
                                 sim_powerBox,           &
                                 sim_radialBox

  implicit none

#include "Flash.h"
#include "constants.h"

  character (len = MAX_STRING_LENGTH) :: fileName

  integer  :: box
  integer  :: fileUnit
  integer  :: nEntries, nTotalEntries
  integer  :: ray
  integer  :: ut_getFreeFileUnit

  real     :: power
  real     :: radius
!
!
!   ...Do the printout only on the master processor.
!
!
  if (sim_globalMe /= MASTER_PE) then
      return
  end if
!
!
!   ...Open the printout file.
!
!
  fileUnit = ut_getFreeFileUnit ()
  fileName = trim (sim_baseName) // "Results.dat"

  open (fileUnit, file = fileName)
!
!
!   ...Print out the main title. 
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "      LASER QUADRATIC TUBE RING SIMULATION PRINTOUT"
  write (fileUnit,'(/)')
!
!
!     ...Print out the general info.
!
!
  write (fileUnit,'(1x,a,es20.10)') "                       Tube center x-coordinate (cm) = ",sim_xw
  write (fileUnit,'(1x,a,es20.10)') "                       Tube center z-coordinate (cm) = ",sim_zw
  write (fileUnit,'(1x,a,es20.10)') "                     Tube critical x-coordinate (cm) = ",sim_xc
  write (fileUnit,'(1x,a,es20.10)') "                     Tube critical z-coordinate (cm) = ",sim_zc
  write (fileUnit,'(1x,a,es20.10)') "                             Focal y-coordinate (cm) = ",sim_yfocal
  write (fileUnit,'(1x,a,es20.10)') "            Electron density at tube center (#/cm^3) = ",sim_nw
  write (fileUnit,'(1x,a,es20.10)') "                  Critical electron density (#/cm^3) = ",sim_nc
  write (fileUnit,'(1x,a,es20.10)') "        Quadratic electron density A factor (#/cm^5) = ",sim_A
  write (fileUnit,'(1x,a,es20.10)') "             Electron temperature at tube center (K) = ",sim_Tw
  write (fileUnit,'(1x,a,es20.10)') "   Inverse-Bremsstrahlung rate at tube center (/sec) = ",sim_nuw
  write (fileUnit,'(1x,a,es20.10)') "             Analytic power at focal point (erg/sec) = ",sim_powerDecayFactor
  write (fileUnit,'(1x,a,i2)'     ) "                         Simulation refinement level = ",sim_refinementLevel
  write (fileUnit,'(1x,a,i10)'    ) "                             Number of rays launched = ",sim_numRaysLaunched
  write (fileUnit,'(1x,a,i10)'    ) "                             Number of rays analyzed = ",sim_numRaysAnalyzed
!
!
!   ...Print out the results title. 
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "                     RESULTS"
  write (fileUnit,'(/)')
!
!
!     ...Print out the results.
!
!
  write (fileUnit,'(1x,a,es20.10)') "                             Beam target radius (cm) = ",sim_beamTargetRadius
  write (fileUnit,'(1x,a,es20.10)') "                             Focal point radius (cm) = ",sim_focalPointRadius
  write (fileUnit,'(1x,a,es20.10)') "              Maximum square radial deviation (cm^2) = ",sim_maxDeltaRsqr
  write (fileUnit,'(1x,a,es20.10)') "          Maximum absolute power deviation (erg/sec) = ",sim_maxDeltaP

  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "    RADIAL BOX CONTENT"
  write (fileUnit,'(/)')
!
!
!     ...Print out the power boxes.
!
!
  nTotalEntries = 0

  do box = 1, sim_totalNumberOfBoxes
     nEntries = sim_powerBox (box)
     nTotalEntries = nTotalEntries + nEntries
     write (fileUnit,'(1x,a,i8,a,i8,a)') " Box Nr ",box," has ",nEntries," entries"
  end do
  write (fileUnit,'(/)')
  write (fileUnit,'(1x,a,i12)') " Total number of entries in all boxes = ",nTotalEntries

  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "    POWER BOX CONTENT PRINTOUT FOR DATA FILE (Power / Count pairs)"
  write (fileUnit,'(/)')

  do box = 1, sim_totalNumberOfBoxes
     nEntries = sim_powerBox (box)
     power    = sim_boxPower (box)
     write (fileUnit,'(1x,a,es20.10,i8)') " ",power,nEntries
  end do
!
!
!     ...Print out the radial boxes.
!
!
  if (1 == 0) then

  nTotalEntries = 0

  do box = 1, sim_totalNumberOfBoxes
     nEntries = sim_radialBox (box)
     nTotalEntries = nTotalEntries + nEntries
     write (fileUnit,'(1x,a,i8,a,i8,a)') " Box Nr ",box," has ",nEntries," entries"
  end do
  write (fileUnit,'(/)')
  write (fileUnit,'(1x,a,i12)') " Total number of entries in all boxes = ",nTotalEntries

  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "    RADIAL BOX CONTENT PRINTOUT FOR DATA FILE (Radius / Count pairs)"
  write (fileUnit,'(/)')

  do box = 1, sim_totalNumberOfBoxes
     nEntries = sim_radialBox (box)
     radius   = sim_boxRadius (box)
     write (fileUnit,'(1x,a,es20.10,i8)') " ",radius,nEntries
  end do

  end if
!
!
!   ...Close the printout file.
!
!
  close (fileUnit)
!
!
!     ...Ready!
!
!
  return
end subroutine sim_printResults
