!!****if* source/Simulation/SimulationMain/ReinickeMeyer/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!
!!***

subroutine Simulation_init()

  use Simulation_data 
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Logfile_interface,           ONLY : Logfile_stampMessage
  use Driver_interface,            ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"

  character(len=MAX_STRING_LENGTH) :: str

  call PhysicalConstants_get('avogadro',sim_avo)
  call PhysicalConstants_get('Boltzmann',sim_boltz)

  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('tinitial',sim_tInitial)    
  call RuntimeParameters_get('eos_singleSpeciesA', sim_singleSpeciesA)
  call RuntimeParameters_get('sim_rfInit', sim_rfInit)
  call RuntimeParameters_get('cond_DensityExponent', sim_DensityExponent)
  call RuntimeParameters_get('cond_TemperatureExponent', sim_TemperatureExponent)
  call RuntimeParameters_get('cond_K0', sim_K0)
  call RuntimeParameters_get('smallt', sim_smallt)
  call RuntimeParameters_get('smlrho', sim_smlrho)

  call RuntimeParameters_get ("geometry", str)
  call RuntimeParameters_mapStrToInt(str, sim_geometry)

  ! This simulation is only set up to run in one of three geometries
  ! at this time:
  !
  ! 1) 1D Spherical
  ! 2) 2D R-z
  ! 3) 3D Cartesian
  !
  ! It could be extended to run in other geometries, but this will
  ! require tweaking Simulation_computeAnalytical
  if(  (NDIM == 1 .and. sim_geometry /= SPHERICAL)   .or. &
       (NDIM == 2 .and. sim_geometry /= CYLINDRICAL) .or. &
       (NDIM == 3 .and. sim_geometry /= CARTESIAN) ) then

     call Logfile_stampMessage("[Simulation_init] ERROR")
     call Logfile_stampMessage("  This simulation only supports the following geometries:")
     call Logfile_stampMessage("    - 1D Spherical")
     call Logfile_stampMessage("    - 2D Cylindrical")
     call Logfile_stampMessage("    - 3D Cartesian")
     call Driver_abortFlash("[Simulation_init] Unsupported Geomtry SEE LOG")
  end if

end subroutine Simulation_init
