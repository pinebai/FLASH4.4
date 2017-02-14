!!****if* source/Simulation/SimulationMain/radflaHD/EnergyXchange/Simulation_init
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
!!  Initializes all the parameters needed for a particular simulation
!!
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp

  use Simulation_interface, ONLY: Simulation_mapIntToStr
  
  use Driver_interface, ONLY: Driver_getMype

  implicit none
#include "Flash.h"
#include "constants.h"

  real :: tot_massfrac

  character (len=MAX_STRING_LENGTH) :: rtpar
  character (len=MAX_STRING_LENGTH) :: spec_str

  integer :: i
  integer :: ut_getFreeFileUnit

  call RuntimeParameters_get('sim_rho' , sim_rho )
  call RuntimeParameters_get('sim_tele', sim_tele)
  call RuntimeParameters_get('sim_tion', sim_tion)
  call RuntimeParameters_get('sim_trad', sim_trad)
  call RuntimeParameters_get('smallX', sim_smallX)
  call RuntimeParameters_get("gamma",sim_gamma)

  ! Open file for writing temperature information:
  sim_fileUnitT = ut_getFreeFileUnit()
  open(unit=sim_fileUnitT, file="temperatures.txt", form="formatted", position='append')

  ! Write the file header:
  if(sim_globalME == MASTER_PE) then
     write(sim_fileUnitT,'(a10,7a15)') &
          '#    nstep', 'time (s)', 'dt (s)', 'tion (K)', 'tele (K)', 'trad (K)', 'CV_ele', 'CV_ion'
  end if
  
  ! Open file for writing energy information:
  sim_fileUnitE = ut_getFreeFileUnit()
  open(unit=sim_fileUnitE, file="energies.txt", form="formatted", position='append')

  ! Write the file header:
  if(sim_globalME == MASTER_PE) then
     write(sim_fileUnitE,'(a10,8a15)') &
          '#    nstep', 'time (s)', 'dt (s)', 'eion (erg/cm^3)', 'eele (erg/cm^3)', &
          'eint (erg/cm^3)', 'erad (erg/cm^3)', 'CV_ele', 'CV_ion'
  end if

  ! Get rank:
  call Driver_getMype(GLOBAL_COMM, sim_globalMe)

end subroutine Simulation_init
