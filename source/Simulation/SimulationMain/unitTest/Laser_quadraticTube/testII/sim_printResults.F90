!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testII/sim_printResults
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
!!   square radial deviation and the maximum absolute partition power deviation to a .dat file.
!!   Only the master processor does the printout.
!!
!!***

subroutine sim_printResults ()

  use Simulation_data,    ONLY : sim_A,                &
                                 sim_baseName,         &
                                 sim_beamTargetRadius, &
                                 sim_focalPointRadius, &
                                 sim_globalMe,         &
                                 sim_maxDeltaPZ,       &
                                 sim_maxDeltaRsqr,     &
                                 sim_nc,               &
                                 sim_numRaysAnalyzed,  &
                                 sim_numRaysLaunched,  &
                                 sim_nuw,              &
                                 sim_nw,               &
                                 sim_percentFocus,     &
                                 sim_powerDecayFactor, &
                                 sim_powerPartition,   &
                                 sim_refinementLevel,  &
                                 sim_tcross,           &
                                 sim_Tw,               &
                                 sim_xc,               &
                                 sim_xw,               &
                                 sim_yfocal,           &
                                 sim_Z,                &
                                 sim_zc,               &
                                 sim_zw

  implicit none

#include "Flash.h"
#include "constants.h"

  character (len = MAX_STRING_LENGTH) :: fileName

  integer  :: fileUnit
  integer  :: ray
  integer  :: ut_getFreeFileUnit
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
  write (fileUnit,'(a)') "      LASER QUADRATIC TUBE (TEST II) SIMULATION PRINTOUT"
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
  write (fileUnit,'(1x,a,es20.10)') "                     Rays domain crossing time (sec) = ",sim_tcross
  write (fileUnit,'(1x,a,es20.10)') "     Beam power partition function Z (dimensionless) = ",sim_powerPartition
  write (fileUnit,'(1x,a,es20.10)') "    Analytic Z-scaled power at focal point (erg/sec) = ",sim_powerDecayFactor
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
  write (fileUnit,'(1x,a,f8.2)'   ) "                                   Percent focus (%) = ",sim_percentFocus
  write (fileUnit,'(1x,a,es20.10)') "              Maximum square radial deviation (cm^2) = ",sim_maxDeltaRsqr
  write (fileUnit,'(1x,a,es20.10)') " Maximum absolute Z-scaled power deviation (erg/sec) = ",sim_maxDeltaPZ
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
