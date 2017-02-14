!!****if* source/Simulation/SimulationMain/unitTest/ProtonImaging/CircleDeflection/sim_printResults
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
!!   This routine prints some general info about the simulation. Only the master processor
!!   does the printout.
!!
!!***

subroutine sim_printResults ()

  use Simulation_data,     ONLY : sim_baseName,           &
                                  sim_globalMe,           &
                                  sim_refinementLevel,    &
                                  sim_cellSizeX,          &
                                  sim_cellSizeY,          &
                                  sim_cellSizeZ
                                     

  use ProtonImaging_data,  ONLY : pi_protonCharge, &
                                  pi_protonMass,   &
                                  pi_speedOfLight
  implicit none

#include "Flash.h"
#include "constants.h"

  character (len = MAX_STRING_LENGTH) :: fileName

  integer  :: fileUnit
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
  write (fileUnit,'(a)') "      PROTON IMAGING: SIMULATION PRINTOUT"
  write (fileUnit,'(/)')
!
!
!     ...Print out the general info.
!
!
  write (fileUnit,'(1x,a,es20.12)') "                         Cell size x-coordinate (cm) = ",sim_cellSizeX
  write (fileUnit,'(1x,a,es20.12)') "                         Cell size y-coordinate (cm) = ",sim_cellSizeY
  write (fileUnit,'(1x,a,es20.12)') "                         Cell size z-coordinate (cm) = ",sim_cellSizeZ
  write (fileUnit,'(1x,a,es20.12)') "                                     Proton mass (g) = ",pi_protonMass
  write (fileUnit,'(1x,a,es20.12)') "                                 Proton charge (esu) = ",pi_protonCharge
  write (fileUnit,'(1x,a,es20.12)') "                               Speed of Light (cm/s) = ",pi_speedOfLight
  write (fileUnit,'(1x,a,i2)'     ) "                         Simulation refinement level = ",sim_refinementLevel
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
