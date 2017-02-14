!!****if* source/Simulation/SimulationMain/DustCollapse/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the Spherical dust cloud collapse problem
!!  in cartesian coordinates
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!
!!   sim_initRad                  Initial radius of cloud
!!   sim_initDens                 Initial density of cloud
!!   sim_tAmbient                 Initial ambient temperature (everywhere)
!!   sim_iCtr,sim_jCtr, sim_kCtr  Coordinates of the center of the cloud
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stampMessage
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"
  
  
  real, save :: gasConst
  real :: soundSpeed, dr_prof
  integer      :: i,   N_ext

  ! buffers to write init params to logfile
  ! two dimensions refer to name and value of parameter,
  ! which will be expanded to 'name=value' in Logfile unit

  character(len=MAX_STRING_LENGTH), dimension(5,2) :: param_buff
  character(len=MAX_STRING_LENGTH)                 :: num_to_str
     
  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get("xmin", sim_imin)
  call RuntimeParameters_get("xmax", sim_imax)
  call RuntimeParameters_get("ymin", sim_jmin)
  call RuntimeParameters_get("ymax", sim_jmax)
  call RuntimeParameters_get("zmax", sim_kmax)
  call RuntimeParameters_get("zmin", sim_kmin)
  call RuntimeParameters_get("smallp", sim_smallp)
  call RuntimeParameters_get("smalle", sim_smalle)
  call RuntimeParameters_get("sim_initRad", sim_initRad)
  call RuntimeParameters_get("gamma", sim_gamma)
  call RuntimeParameters_get("sim_tAmbient", sim_tAmbient)
  call RuntimeParameters_get("sim_initDens", sim_initDens)
  call RuntimeParameters_get("smlrho", sim_smlrho)
  call RuntimeParameters_get("sim_ictr", sim_ictr)
  call RuntimeParameters_get("sim_jctr", sim_jctr)
  call RuntimeParameters_get("sim_kctr", sim_kctr)
  call PhysicalConstants_get("ideal gas constant", gasConst)
     


  if (sim_meshMe == MASTER_PE) then
     
     write (*,*)
     call Logfile_stampMessage( "initializing for dust collapse problem")
     write (*,*) 'flash:  initializing for dust collapse problem.'
     write (*,*)
     
  endif
  
  !               Construct the radial samples needed for the initialization.
  
  sim_presFrac = gasConst * Sim_tAmbient
  soundSpeed = sqrt( sim_gamma*sim_presFrac )
     
  dr_prof = sqrt( (sim_imax-sim_imin)**2 + K2D*(sim_jmax-sim_jmin)**2 + & 
       &                    K3D*(sim_kmax-sim_kmin)**2 ) / (N_prof-1)
  N_ext = int( Sim_initRad / dr_prof ) + 1
  
  do i = 1, N_prof
     sim_rProf(i) = (i-1) * dr_prof
  enddo
  
  do i = 1, N_ext
     sim_rhoProf(i) = sim_initDens
     sim_pProf(i)   = sim_presFrac * sim_rhoProf(i)
     sim_vProf(i)   = 0.
     !******************** for initial linear radial infall
     !           sim_vProf(i)   = -3.4E9 * sim_rProf(i)/Sim_initRad
  enddo
  
  do i = N_ext+1, N_prof
     sim_rhoProf(i) = sim_rhoProf(N_ext) / 10.
     sim_pProf(i)   = sim_presFrac * sim_rhoProf(i)
     sim_vProf(i)   = 0.
  enddo
  
  if (sim_meshMe == MASTER_PE) then
     
     write (param_buff(1,1), "(A)") 'sound speed'
     write (num_to_str, "(ES20.13)") soundSpeed ! format once corresponded to REAL_FORMAT as defined in flash_defines.fh
     write (param_buff(1,2), "(A)") trim(adjustl(num_to_str))
     
     write (param_buff(2,1), "(A)") 'central dens'
     write (num_to_str, "(ES20.13)") sim_initDens
     write (param_buff(2,2), "(A)") trim(adjustl(num_to_str))
     
     write (param_buff(3,1), "(A)") 'central pres'
     write (num_to_str, "(ES20.13)") sim_presFrac*sim_initDens
     write (param_buff(3,2), "(A)") trim(adjustl(num_to_str))
     
     write (param_buff(4,1), "(A)") 'external dens'
     write (num_to_str, "(ES20.13)") sim_rhoProf(N_ext+1)
     write (param_buff(4,2), "(A)") trim(adjustl(num_to_str))
     
     write (param_buff(5,1), "(A)") 'external pres'
     write (num_to_str, "(ES20.13)") sim_pProf(N_ext+1)
     write (param_buff(5,2), "(A)") trim(adjustl(num_to_str))
     
!!     call Logfile_stamp(param_buff, 5, 2, 'init_param', "type='list' count='0'")
        
     if (sim_rhoProf(N_ext) < 100.*sim_smlrho) then
!!        call Logfile_stamp( 'init_block:  ext dens close to cutoff', 'warning')
     endif
     !           open (99, file="profiles.dat", status="new")
     !           write (99,*) '#         r         rho         p'
     !           do i = 1, N_prof
     !             write (99,'(3(ES15.8,:,2X))')
     !     &         sim_rProf(i), sim_rhoProf(i), sim_pProf(i)
     !           enddo
     !           close (99)
  endif
  
end subroutine Simulation_init



