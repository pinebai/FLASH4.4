!!****if* source/Simulation/SimulationMain/MGDInfinite/Simulation_init
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

  ! Get the runtime parameter for the mass fraction of each species:
  tot_massfrac = 0.0
  do i = 1, NSPECIES
     ! Get the species name:
     call Simulation_mapIntToStr(i+SPECIES_BEGIN-1,spec_str,MAPBLOCK_UNK)

     ! Form the runtime parameter name:
     write(rtpar,'(3a)') "sim_", trim(spec_str), "MassFrac"
     call RuntimeParameters_get(rtpar, sim_massfracs(i))

     if(sim_massfracs(i) < sim_smallX) then
        call Driver_abortFlash("[Simulation_init] Improper mass fraction. Set " // &
             trim(rtpar) // " to be >= smallX and <= 1.0")
     end if

     tot_massfrac = tot_massfrac + sim_massfracs(i)
  end do

  ! Normalize mass fractions to 1:
  sim_massfracs(:) = sim_massfracs(:) / tot_massfrac


  ! Open file for writing temperature information:
  sim_fileUnit = ut_getFreeFileUnit()
  open(unit=sim_fileUnit, file="temperatures.txt", form="formatted", position='append')

  ! Write the file header:
  if(sim_globalME == MASTER_PE) then
     write(sim_fileUnit,'(a10,5a15)') &
          '#    nstep', 'time (s)', 'dt (s)', 'tion (K)', 'tele (K)', 'trad (K)'
  end if
  
  ! Get rank:
  call Driver_getMype(GLOBAL_COMM, sim_globalMe)

end subroutine Simulation_init
