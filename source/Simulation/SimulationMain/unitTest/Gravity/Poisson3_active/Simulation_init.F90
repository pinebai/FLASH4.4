!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson3_active/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes data for the Maclaurin spheroid problem.
!!  Analytical solution given in Ricker, P. 2006 "A Direct Multigrid Poisson Solver for Oct-tree adaptive meshes"
!!
!!  References:  Maclaurin, C. 1742, 
!!               Chandrasekhar, S. 1987, Ellipsoidal Figures of Equilibrium
!!
!! ARGUMENTS
!!
!!  
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use Grid_interface, ONLY : Grid_getGeometry
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
       RuntimeParameters_set
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

#include "constants.h"
#include "Flash.h"

  
  integer :: geometry
  real    :: amean, p0, mass, AA3


  call Driver_getMype(MESH_COMM, sim_meshMe)

  call PhysicalConstants_get("pi", sim_pi)
  call PhysicalConstants_get("Newton", sim_Newton)

  call RuntimeParameters_get("eccentricity", sim_e)
  call RuntimeParameters_get("equatorial_semimajor_axis", sim_a1)
  call RuntimeParameters_get("angular_velocity", sim_Omega1)
  call RuntimeParameters_get("density", sim_density)
  call RuntimeParameters_get("gamma", sim_gamma)
  call RuntimeParameters_get("smallx", sim_smallX)
  call RuntimeParameters_get("smlrho", sim_smallRho)
  call RuntimeParameters_get("smallp", sim_smallP)
  call RuntimeParameters_get("smalle", sim_smallE)
  call RuntimeParameters_get("xctr", sim_xctr)
  call RuntimeParameters_get("yctr", sim_yctr)
  call RuntimeParameters_get("zctr", sim_zctr)
  call RuntimeParameters_get("nsubzones", sim_nsubzones)

 call RuntimeParameters_get("pass_tolerance", sim_passTolerance)
 call RuntimeParameters_get("basenm",sim_baseName)
  

  ! Compute derived quantities.

  sim_nsubinv = 1./sim_nsubzones

  sim_a3  = sim_a1 * sqrt(1.0-sim_e**2)   ! second axis of spheroid
  sim_a1inv  = 1.0/sim_a1
  sim_a3inv  = 1.0/sim_a3

  if (sim_e > 1.E-10) then                ! ellipsoid
     AA3 = (2.0*sqrt(1.0-sim_e**2)/sim_e**2) * &
          (1.0/sqrt(1.0-sim_e**2) - asin(sim_e)/sim_e)
  else
     AA3 = 1.0                            ! sphere
  endif

  amean  = sim_a1 * (1.0-sim_e**2)**(1.0/6.0)

  ! if density is -1, calculate a density that will produce a mass of 1.0 in the spheriod
  if (abs(sim_density + 1.0) < tiny(0.0)) then
     sim_density = 0.75 / ( sim_pi * sim_a1**2 * sim_a3)
     call RuntimeParameters_set("density",sim_density)
     if (sim_meshMe == MASTER_PE) print *, 'sim_density set to ',sim_density,' to generate a mass of 1.0'
     call Logfile_stamp(sim_density,'[Simulation_init] Reset density to ')
  end if

  ! some sort of additional angular velocity
  sim_Omega2 = sqrt(sim_pi*sim_Newton*sim_density) * sim_Omega1

  ! pressure in center
  p0     = (2.0/3.0)*sim_pi*sim_Newton*sim_density**2*amean**2
  ! constant pressure around surface of spheroid
  sim_Pconst = p0 * AA3 * (1.0-sim_e**2)**(2.0/3.0)
  ! mass of total spheroid
  mass   = (4.0/3.0)*sim_pi*sim_density*sim_a1**3*sqrt(1.0-sim_e**2)


  ! Determine location of spheroid center depending on grid geometry.
  ! Allowed geometries:  2D:  axisymmetric; 3D:  Cartesian.
  ! Eventually should also support 3D cylindrical and spherical.

  call Grid_getGeometry(geometry)
  if ((NDIM == 2) .and. (geometry == POLAR)) then                      
     ! 2D axisymmetric
     sim_initGeometry = geometry
     sim_xctr = 0.0
     sim_yctr = 0.0
     sim_zctr = 0.0
     call RuntimeParameters_set("xctr", sim_xctr)
     call RuntimeParameters_set("yctr", sim_yctr)
     call RuntimeParameters_set("zctr", sim_zctr)

  else if ((NDIM == 3) .and. (geometry == CARTESIAN)) then
     ! 3D Cartesian
     sim_initGeometry = geometry

  else            ! unsupported geometry
     call Driver_abortFlash('Simulation_init:  unsupported geometry.  Needs CYLINDRICAL or 3D SPHERICAL')
  endif


  ! Write a message to stdout describing the problem setup.

  if (sim_meshMe == MASTER_PE) then

     write (*,*)
     call Logfile_stamp( "initializing for unitTest/Gravity/Poisson3 problem", & 
          "[Simulation_init]")
     write (*,*) 'Simulation_init:  initializing for unitTest/Gravity/Poisson3 problem.'
     write (*,*)
     write (*,*) "density      = ", sim_density
     write (*,*) "a1           = ", sim_a1
     write (*,*) "eccentricity = ", sim_e
     write (*,*) "Omega1       = ", sim_Omega1
     write (*,*)
     write (*,*) "a3           = ", sim_a3
     write (*,*) "amean        = ", amean
     write (*,*) "p0           = ", p0
     write (*,*) "mass         = ", mass
     write (*,*) "omega2       = ", sim_omega2
     write (*,*)
     write (*,*) "xctr         = ", sim_xctr
     write (*,*) "yctr         = ", sim_yctr
     write (*,*) "zctr         = ", sim_zctr
     write (*,*)
     write (*,*) "nsubzones    = ", sim_nsubzones
     write (*,*)

  endif

  return
end subroutine Simulation_init
