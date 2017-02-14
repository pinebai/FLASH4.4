!!****if* source/Simulation/SimulationMain/unitTest/Multipole/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  Simulation_init ()
!!
!! DESCRIPTION
!!
!!  Prepares initialization for the multipole solver unit test, based on the MacLaurin
!!  spheroid problem.
!!
!!  References:  Maclaurin, C. 1742,
!!               Chandrasekhar, S. 1987, Ellipsoidal Figures of Equilibrium
!!
!! ARGUMENTS
!!
!!***

subroutine Simulation_init ()

  use Simulation_data
  use Driver_data,                 ONLY : dr_globalMe
  use Grid_interface,              ONLY : Grid_getGeometry
  use Driver_interface,            ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                          RuntimeParameters_set
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Logfile_interface,           ONLY : Logfile_stamp

  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: geometry

  call PhysicalConstants_get ("pi",                        sim_pi)
  call PhysicalConstants_get ("Newton",                    sim_Newton)

  call RuntimeParameters_get ("eccentricity",              sim_e)
  call RuntimeParameters_get ("equatorialSemimajorAxis",   sim_a1)
  call RuntimeParameters_get ("density",                   sim_density)
  call RuntimeParameters_get ("smlrho",                    sim_smallRho)
  call RuntimeParameters_get ("xctr",                      sim_xctr)
  call RuntimeParameters_get ("yctr",                      sim_yctr)
  call RuntimeParameters_get ("zctr",                      sim_zctr)
  call RuntimeParameters_get ("nsubzones",                 sim_nsubzones)
  call RuntimeParameters_get ("passTolerance",             sim_passTolerance)
! 
!
!   ...Compute derived quantities for the MacLaurin spheroid.
!
!
  sim_nsubinv = 1./sim_nsubzones
  sim_a3      = sim_a1 * sqrt(1.0-sim_e**2)
  sim_a1inv   = 1.0/sim_a1
  sim_a3inv   = 1.0/sim_a3
! 
!
!   ...Set the geometry.
!
!
  call Grid_getGeometry (geometry)

  if ((NDIM == 3) .and. (geometry == CARTESIAN)) then
      sim_initGeometry = GRID_3DCARTESIAN
  else if ((NDIM == 3) .and. (geometry == CYLINDRICAL)) then
      sim_initGeometry = GRID_3DCYLINDRICAL
  else if ((NDIM == 2) .and. (geometry == CYLINDRICAL)) then
      sim_initGeometry = GRID_2DCYLINDRICAL
  else if ((NDIM == 2) .and. (geometry == SPHERICAL)) then
      sim_initGeometry = GRID_2DSPHERICAL
  else if ((NDIM == 1) .and. (geometry == SPHERICAL)) then
      sim_initGeometry = GRID_1DSPHERICAL
  else
      call Driver_abortFlash ('MacLaurin unit test: unsupported geometry!')
  endif
! 
!
!   ...Write a message to stdout describing the MacLaurin spheroid.
!
!
  if (dr_globalMe .eq. MASTER_PE) then

     call Logfile_stamp ("Initializing for multipole unit test (MacLaurin problem)","[Simulation_init]")

     write (*,*)
     write (*,*) "MacLaurin spheroid data."
     write (*,*) "------------------------"
     write (*,*)
     write (*,*) "Density      = ", sim_density
     write (*,*) "eccentricity = ", sim_e
     write (*,*) "a1           = ", sim_a1
     write (*,*) "a3           = ", sim_a3
     write (*,*)

  endif
!
!
!   ...Ready!
!
!
  return
end subroutine Simulation_init
