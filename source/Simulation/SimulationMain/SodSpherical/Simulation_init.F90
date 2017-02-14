!!****if* source/Simulation/SimulationMain/SodSpherical/Simulation_init
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
!!
!! DESCRIPTION
!!
!!  Setup a sod-like problem in spherical coordinates to test whether a planar
!!  shock in spherical coordinates stays planar.  This effectively tests the 
!!  fictitous forces in force().  If the forces are setup right, then the planar
!!  shock should stay planar.
!!
!!  Right now, this is setup to do the problem in 2-d spherical coordinates.  
!!  sim_idir = 1 is the x-direction.  sim_idir = 2 is the z-direction.  
!!
!! 
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!  sim_rhoLeft    Density on the the left side of the interface
!!  sim_rhoRight   Density on the right side of the interface
!!  sim_pLeft      Pressure on the left side of the interface
!!  sim_pRight     Pressure on the righ side of the interface
!!  sim_shockpos   Point of intersection between the shock surface and the x-axis
!!  sim_idir       the direction along which to propagate the shock.
!!
!!  gamma          ideal gas adiabatic index Cp/Cv
!!  cfl            Courant-Friedrichs-Loewy factor (as used by Hydro)
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use Logfile_interface, ONLY : Logfile_stamp
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get

  implicit none

#include "Flash.h"
#include "constants.h"

  


  character(len=MAX_STRING_LENGTH) :: str_geometry
  integer :: meshGeom

  call RuntimeParameters_get('smallx', sim_smallX) 
  
!!  call RuntimeParameters_get('gamma', sim_gamma)           ! unused
  
  call RuntimeParameters_get('sim_rhoLeft', sim_rhoLeft)
  call RuntimeParameters_get('sim_rhoRight', sim_rhoRight)
  
  call RuntimeParameters_get('sim_pLeft', sim_pLeft)
  call RuntimeParameters_get('sim_pRight', sim_pRight)
  
  call RuntimeParameters_get('sim_shockpos', sim_shockpos)
  call RuntimeParameters_get('sim_idir', sim_idir)

#ifdef CFL_VAR
  call RuntimeParameters_get('cfl', sim_hydroCfl)
#endif

  call RuntimeParameters_get("geometry", str_geometry)
  call RuntimeParameters_mapStrToInt(str_geometry, meshGeom)
  if (meshGeom /= SPHERICAL .AND. NDIM /= 2) then 
     call Driver_abortFlash("ERROR: invalid geometry for SodSpherical")
  endif

  call Logfile_stamp( "initializing SodSpherical problem",  &
       "[Simulation_init]")
     
  ! convert the shock angle paramters
!!  sim_xAngle = sim_xAngle * 0.0174532925 ! Convert to radians.
!!  sim_yAngle = sim_yAngle * 0.0174532925

!!  sim_xCos = cos(sim_xAngle)
  
  if (NDIM == 1) then
!!     sim_xCos = 1.
!!     sim_yCos = 0.
!!     sim_zCos = 0.
     
  elseif (NDIM == 2) then
!!     sim_yCos = sqrt(1. - sim_xCos*sim_xCos)
!!     sim_zCos = 0.
     
  elseif (NDIM == 3) then
!!     sim_yCos = cos(sim_yAngle)
!!     sim_zCos = sqrt( max(0., 1. - sim_xCos*sim_xCos - sim_yCos*sim_yCos) )
  endif
end subroutine Simulation_init







