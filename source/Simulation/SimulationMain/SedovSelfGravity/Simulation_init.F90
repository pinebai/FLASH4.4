!!****if* source/Simulation/SimulationMain/SedovSelfGravity/Simulation_init
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
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for Sedov Spherical Explosion 
!!  problem.
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!  p_ambient       Initial ambient pressure
!!  rho_Ambient     Initial ambient density
!!  exp_energy      Explosion energy (distributed over 2^dimen central zones)
!!  t_init          Initial time since explosion
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!***
subroutine Simulation_init()
  use Simulation_data, ONLY : sim_nsubzones, sim_inSubzm1,sim_inSubzones,&
                              sim_rProf, sim_vProf, sim_pProf, sim_rhoProf,&
                              sim_drProf,sim_gamma,SIM_NPROFILE,&
                              sim_smlrho, sim_smallp, sim_xn, sim_meshMe
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Grid_interface, ONLY : Grid_getGeometry

  implicit none

#include "constants.h"
#include "Flash.h"

  real, save  :: pi, vctr, smallx
  real, save    :: p_ambient, rho_ambient, exp_energy
  real, save    :: t_init,  p_exp, r_init
  real, save    :: imax
  integer :: i
  integer :: geom
  character(len=80), save :: prof_file
 
  call Grid_getGeometry(geom)
  call Driver_getMype(MESH_COMM,sim_meshMe)

  if (geom /= SPHERICAL)  then
     print *,'Error -- works only for spherical geometry'
     call Driver_abortFlash('This setup requires spherical geometry')
  endif
  
  
  pi = PI
  
  call RuntimeParameters_get( 'p_ambient',   p_ambient)
  call RuntimeParameters_get( 'rho_ambient', rho_ambient)
  call RuntimeParameters_get( 'exp_energy',  exp_energy)
  
  call RuntimeParameters_get( 'r_init',      r_init)
  call RuntimeParameters_get( 't_init',      t_init)
  call RuntimeParameters_get( 'gamma',       sim_gamma)
  
  call RuntimeParameters_get( 'smallx',      smallx)
  call RuntimeParameters_get( 'smlrho',      sim_smlrho)
  call RuntimeParameters_get( 'smallp',      sim_smallp)
  
  call RuntimeParameters_get( 'prof_file',   prof_file)

#if NSPECIES > 0  
  sim_xn(2:NSPECIES)     = smallx
  sim_xn(1) = 1.e0 - (NSPECIES-1)*smallx
#endif
  
  call RuntimeParameters_get('xmax', imax)
  
  call RuntimeParameters_get( 'sim_nsubzones', sim_nsubzones)
  
  if (sim_nsubzones .le. 1) sim_nsubzones = 2
  
  sim_inSubzones = 1./real(sim_nsubzones)
  sim_inSubzm1   = 1./real(sim_nsubzones-1)
  
  !
  !  Calculate the initial volume and interior pressure.
  !
  if ( r_init > 0.e0 ) then
     vctr  = 4./3.*pi*r_init**3
     p_exp = (sim_gamma-1.) * exp_energy / vctr
  else
     p_exp = p_ambient
  end if
  
  !
  !  Write a message to stdout describing the problem setup.
  !
  if (sim_meshMe == MASTER_PE) then
     
     write (*,*)
     call Logfile_stampMessage( "initializing for sedov problem")
     write (*,*) 'flash:  initializing for sedov problem.'
     write (*,*)
     write (*,*) 'p_ambient  = ', p_ambient
     write (*,*) 'rho_ambient= ', rho_ambient
     write (*,*) 'gamma      = ', sim_gamma
     write (*,*) 'exp_energy = ', exp_energy 
     write (*,*) 'r_init     = ', r_init
     write (*,*) 'p_exp      = ', p_exp
     write (*,*) 'nsubzones  = ', sim_nsubzones
     write (*,*)
     
     if (t_init .gt. 0.) then
        call Logfile_stampMessage & 
             ( "t_init > 0:  sedov solution currently broken")
     endif
     
  endif
  sim_drProf = imax / (SIM_NPROFILE-1)
  
  do i = 1, SIM_NPROFILE
     sim_rProf(i)   = (i-1) * sim_drProf
  enddo
  
  !
  !  If t>0, use the analytic Sedov solution to initialize the
  !  code.  Otherwise, just use a top-hat.
  !
  
  if (t_init .gt. 0.) then
     
     call sim_setAnalyticSedov (SIM_NPROFILE, sim_rProf, sim_rhoProf, &
          sim_pProf, sim_vProf, t_init, sim_gamma, exp_energy,&
          p_ambient, rho_ambient)

  else

     do i = 1, SIM_NPROFILE
        sim_rhoProf(i) = rho_ambient
        sim_pProf(i)   = p_ambient
        sim_vProf(i)   = 0.e0
        if (sim_rProf(i) .le. r_init) sim_pProf(i)   = p_exp
        if (sim_rProf(i) .gt. 1.0e0 ) sim_rhoProf(i) = 1.e-15*sim_rhoProf(i)
     enddo
          
  endif

  open(unit=1,file="sedov_ profile",status='unknown')
  write(1,*) sim_drProf*1000. - 1.41562919156465981
  write(1,*) p_exp - 681.16756574500903
  
  do i = 1, SIM_NPROFILE
     write(1,*) i, sim_drProf*i, sim_rhoProf(i),sim_vProf(i),sim_pProf(i)
  enddo
  close(1)
end subroutine Simulation_init
