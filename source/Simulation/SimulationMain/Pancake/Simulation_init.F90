!!****if* source/Simulation/SimulationMain/Pancake/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get routine for initialization.
!!  Initializes initial conditions for Particle Orbit problem.
!!
!! ARGUMENTS
!!
!!   
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, & 
                                          RuntimeParameters_set
  use PhysicalConstants_interface, ONLY :  PhysicalConstants_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Cosmology_data, ONLY : csm_scaleFactor

  implicit none

#include "constants.h"
#include "Flash.h"

  
  real :: Curv, Hubinit, vol, scale
  character(len=128) :: str_buffer

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('xmin',sim_xMin)
  call RuntimeParameters_get('ymin',sim_yMin)
  call RuntimeParameters_get('zmin',sim_zMin)
  call RuntimeParameters_get('xmax',sim_xMax)
  call RuntimeParameters_get('ymax',sim_yMax)
  call RuntimeParameters_get('zmax',sim_zMax)

  call RuntimeParameters_get('HubbleConstant',sim_hubble)
  call RuntimeParameters_get('CosmologicalConstant',sim_lambda)
  call RuntimeParameters_get('OmegaMatter',sim_omegaMatter)
  call RuntimeParameters_get('OmegaBaryon',sim_omegaBaryon)
  call RuntimeParameters_get('useParticles',sim_useParticles)
  call RuntimeParameters_get('zfiducial',sim_zfiducial)
  call RuntimeParameters_get('Tfiducial',sim_Tfiducial)
  call RuntimeParameters_get('zcaustic',sim_zcaustic)
  call RuntimeParameters_get('lambda',sim_wavelength)
  call RuntimeParameters_get('xangle',sim_xangle)
  call RuntimeParameters_get('yangle',sim_yangle)

  call RuntimeParameters_get('smlrho',sim_smallRho)
#ifdef PRES_VAR
  call RuntimeParameters_get('smalle',sim_smallE)
  call RuntimeParameters_get('smallp',sim_smallP)
  call RuntimeParameters_get('gamma',sim_gamma)
#endif

  call RuntimeParameters_get('pt_numX', sim_nxp)
  call RuntimeParameters_get('pt_numY', sim_nyp)
  call RuntimeParameters_get('pt_numZ', sim_nzp)

  call PhysicalConstants_get('Newton', sim_Newton)
  call PhysicalConstants_get('pi', sim_pi)
  call PhysicalConstants_get('ideal gas constant', sim_gascon)


! Do some checking.

! If dark matter is included, we must have particles compiled in and turned on.
! If we don't find particles, treat everything as gas.

  if ((sim_omegaBaryon < sim_omegaMatter) .and. (.not. sim_useParticles)) then
    if (sim_meshMe == MASTER_PE) then
      call Logfile_stamp &
        (sim_meshMe, 'WARNING:  sim_omegaBaryon < sim_omegaMatter, but no particles!', & 
               '[SIMULATION Simulation_init]')
      write (str_buffer, *)  &
        'setting sim_omegaMatter = sim_omegaBaryon = ', sim_omegaBaryon
      call Logfile_stamp( str_buffer, '[SIMULATION Simulation_init]')
    endif
    sim_omegaMatter = sim_omegaBaryon
    call RuntimeParameters_set("OmegaMatter", sim_omegaMatter)
  endif

! If gas (baryonic matter) is included, we must have hydrodynamics compiled in.
! There is no 'hydro switch' at the moment, so we check for the existence of
! a mesh variable required for hydro:  the pressure.  If we don't find gas,
! treat everything as dark matter.

#ifndef PRES_VAR
  if (sim_omegaBaryon > 0.) then
     if (sim_meshMe == MASTER_PE) then
        call Logfile_stamp( 'WARNING: sim_omegaBaryon > 0, but no gas!', &
             '[SIMULATION Simulation_init]')
        write (str_buffer, *)  'setting sim_omegaBaryon = 0'
        call Logfile_stamp( str_buffer, '[SIMULATION Simulation_init]')
     endif
     sim_omegaBaryon = 0.
     call RuntimeParameters_set("OmegaBaryon", sim_omegaBaryon)
  endif
#endif

! There can't be more density in baryons than there is in matter.  Limit
! sim_omegaBaryon to be no larger than sim_omegaMatter.
        
  if (sim_omegaBaryon > sim_omegaMatter) then
     if (sim_meshMe == MASTER_PE) then
        write (*,*) 'Simulation_init:  WARNING:  sim_omegaBaryon > sim_omegaMatter!'
        write (*,*) 'Simulation_init:  setting sim_omegaBaryon = sim_omegaMatter = ', &
             sim_omegaMatter
        call Logfile_stamp( 'WARNING:  sim_omegaBaryon > sim_omegaMatter!', & 
             '[SIMULATION Simulation_init]')
        write (str_buffer,*) &
             'setting sim_omegaBaryon = sim_omegaMatter = ', sim_omegaMatter
        call Logfile_stamp( str_buffer, '[SIMULATION Simulation_init]')
     endif
     sim_omegaBaryon = sim_omegaMatter
     call RuntimeParameters_set("OmegaBaryon", sim_omegaBaryon)
  endif
  
! Compute some derived quantities that we will need.

  scale = csm_scaleFactor
  sim_zinitial = (1.0/scale) - 1.0

  Curv = sim_omegaMatter + sim_lambda - 1.
  Hubinit = sim_hubble * sqrt(sim_omegaMatter/scale**3 - Curv/scale**2 + &
       sim_lambda)
  sim_critden  = (3.*sim_hubble**2 / (8.*sim_pi*sim_Newton))
  sim_meanden  = sim_omegaMatter * sim_critden
  sim_bmeanden = sim_omegaBaryon * sim_critden
  sim_dmeanden = (sim_omegaMatter-sim_omegaBaryon)*sim_critden

  sim_kmag = 2.*sim_pi / sim_wavelength

  sim_xangle = sim_xangle * 0.0174532925 ! Convert to radians.
  sim_yangle = sim_yangle * 0.0174532925
      
  sim_xcos = cos(sim_xangle)
      
  if (NDIM .eq. 1) then
    sim_xcos = 1.
    sim_ycos = 0.
    sim_zcos = 0.
    sim_yMax = 1.
    sim_yMin = 0.
    sim_zMax = 1.
    sim_zMin = 0.
    call RuntimeParameters_set("ymin", sim_yMin)
    call RuntimeParameters_set("ymax", sim_yMax)
    call RuntimeParameters_set("zmin", sim_zMin)
    call RuntimeParameters_set("zmax", sim_zMax)
  elseif (NDIM .eq. 2) then
    sim_ycos = sqrt(1. - sim_xcos*sim_xcos)  
    sim_zcos = 0.
    sim_zMax = 1.
    sim_zMin = 0.
    call RuntimeParameters_set("zmin", sim_zMin)
    call RuntimeParameters_set("zmax", sim_zMax)
  elseif (NDIM .eq. 3) then
    sim_ycos = cos(sim_yangle)
    sim_zcos = sqrt( max(0., 1. - sim_xcos*sim_xcos - sim_ycos*sim_ycos) )
  endif

! Particles stuff:

  if (NDIM == 3) then 
     vol=(sim_xMax-sim_xMin)*(sim_yMax-sim_yMin)*(sim_zMax-sim_zMin)
  else if (NDIM == 2) then 
     vol=(sim_xMax-sim_xMin)*(sim_yMax-sim_yMin)
     sim_nzp = 1
  else if (NDIM == 1) then 
     vol=(sim_xMax-sim_xMin)
     sim_nzp = 1
     sim_nyp = 1
  else 
     call Driver_abortFlash("Simulation_init: NDIM is not valid")
  endif
  sim_mr = sim_dmeanden*vol / (float(sim_nxp)*float(sim_nyp)*float(sim_nzp))

! Factors for Lagrangian coordinate computation

  sim_temp1 = (1.0+sim_zcaustic) / (1.0+sim_zinitial)
  sim_temp2 = sim_temp1/sim_kmag
  sim_temp3 = (1.0+sim_zcaustic) / (1.0+sim_zfiducial)

! Let the user know what's going on.

  if (sim_meshMe == MASTER_PE) then
    write (*,*)
    write (*,*) 'Simulaiton_init:  initializing for zeldovich pancake problem.'
    write (*,*)
    write (*,*) 'k    = ', sim_kmag
    write (*,*)
    write (*,*) 'initial scale factor                 = ', scale
    write (*,*) 'initial Hubble parameter             = ', Hubinit
    write (*,*) 'present comoving critical density    = ', sim_critden
    write (*,*) 'initial comoving mean density        = ', sim_meanden
    write (*,*) 'initial comoving mean baryon density = ', sim_bmeanden
    write (*,*)
    write (*,*) 'number of particles                  = ', sim_nxp*sim_nyp*sim_nzp
    write (*,*) 'particle mass resolution             = ', sim_mr
    
    call Logfile_stamp( 'initializing for zeldovich pancake problem.', &
                             '[SIMULATION Simulation_init]')
    write (str_buffer,*) 'k    = ', sim_kmag
    call Logfile_stamp( str_buffer, '[SIMULATION Simulation_init]')
    write (str_buffer,*) 'initial scale factor                 = ', scale
    call Logfile_stamp( str_buffer, '[SIMULATION Simulation_init]')
    write (str_buffer,*) 'initial Hubble parameter             = ', Hubinit
    call Logfile_stamp( str_buffer, '[SIMULATION Simulation_init]')
    write (str_buffer,*) 'present comoving critical density    = ', sim_critden
    call Logfile_stamp( str_buffer, '[SIMULATION Simulation_init]')
    write (str_buffer,*) 'initial comoving mean density        = ', sim_meanden
    call Logfile_stamp( str_buffer, '[SIMULATION Simulation_init]')
    write (str_buffer,*) 'initial comoving mean baryon density = ', sim_bmeanden
    call Logfile_stamp( str_buffer, '[SIMULATION Simulation_init]')
    write (str_buffer,*) 'number of particles                  = ', sim_nxp*sim_nyp*sim_nzp
    call Logfile_stamp( str_buffer, '[SIMULATION Simulation_init]')
    write (str_buffer,*) 'particle mass resolution             = ', sim_mr
    call Logfile_stamp( str_buffer, '[SIMULATION Simulation_init]')
  endif

end subroutine Simulation_init
