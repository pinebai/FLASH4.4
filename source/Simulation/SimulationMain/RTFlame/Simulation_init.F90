!!****if* source/Simulation/SimulationMain/RTFlame/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! AUTOGENROBODOC
!!
!!
!!***

! Dean Townsley 2009

subroutine Simulation_init()
  
  use Simulation_data
  use Flame_interface, ONLY : Flame_rhJump
  use Grid_interface, ONLY : Grid_getGlobalIndexLimits, &
    Grid_getMinCellSize
  use hse_interface, ONLY: sim_hse_step, flame_hse
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use fl_effData, ONLY: fl_effDeltae, fl_eff_ye_u, fl_eff_sumy_u, fl_eff_ye_b, fl_eff_sumy_b
  use IO_interface, ONLY: IO_setScalar, IO_getScalar
  use Driver_data, ONLY: dr_restart
  use Driver_interface, ONLY : Driver_getMype
  implicit none
#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  integer, dimension(MDIM) :: globalIndexLimits
  integer  :: nx, myPE
  real     :: deltax

!===========================================================================

  call Driver_getMype(MESH_COMM, myPE)
  call RuntimeParameters_get("gconst",sim_grav)
  call RuntimeParameters_get("vel_pert_amp", sim_vel_pert_amp)
  call RuntimeParameters_get("vel_pert_wavelength1", sim_vel_pert_wavelength1)

  call RuntimeParameters_get("spert_ampl1", sim_spert_ampl1)
  call RuntimeParameters_get("spert_wl1", sim_spert_wl1)
  call RuntimeParameters_get("spert_phase1", sim_spert_phase1)
  call RuntimeParameters_get("spert_ampl2", sim_spert_ampl2)
  call RuntimeParameters_get("spert_wl2", sim_spert_wl2)
  call RuntimeParameters_get("spert_phase2", sim_spert_phase2)

  if (myPE == MASTER_PE) then
     if (dr_restart) then
        call IO_getScalar("last_burned_mass", sim_last_burned_mass)
     else
        sim_last_burned_mass = -1.0
     endif
  endif

  ! calculate some info obout unburned and burned states
  call RuntimeParameters_get("temp_unburned", sim_temp_u)
  call RuntimeParameters_get("dens_unburned", sim_dens_u)
  sim_ye_u = fl_eff_ye_u
  sim_sumy_u = fl_eff_sumy_u
  sim_ye_b = fl_eff_ye_b
  sim_sumy_b = fl_eff_sumy_b

  sim_eosData_u(EOS_DENS) = sim_dens_u
  sim_eosData_u(EOS_TEMP) = sim_temp_u
  sim_eosData_u(EOS_ABAR) = 1.e0 / fl_eff_sumy_u
  sim_eosData_u(EOS_ZBAR) = fl_eff_ye_u * sim_eosData_u(EOS_ABAR)

  sim_eosData_b(EOS_DENS) = sim_dens_u
  sim_eosData_b(EOS_TEMP) = sim_temp_u
  sim_eosData_b(EOS_ABAR) = 1.e0 / fl_eff_sumy_b
  sim_eosData_b(EOS_ZBAR) = fl_eff_ye_b * sim_eosData_b(EOS_ABAR)

  call Flame_rhJump(sim_eosData_u, sim_eosData_b, &
                fl_effDeltae, 0.0, MODE_DENS_TEMP)

  sim_dens_b = sim_eosData_b(EOS_DENS)
  sim_temp_b = sim_eosData_b(EOS_TEMP)

  if (myPE == MASTER_PE) then
     write (6,*) "atwood number = ", (sim_dens_u-sim_dens_b)/(sim_dens_u+sim_dens_b)
  endif

  ! measure domain and dimension array to hold initial profile
  call Grid_getGlobalIndexLimits(globalIndexLimits)
  ! assume the cells are square for now
  call Grid_getMinCellSize(deltax)

  nx = globalIndexLimits(IAXIS)
  allocate(sim_dens_i(nx))
  allocate(sim_temp_i(nx))
  allocate(sim_ye_i(nx))
  allocate(sim_sumy_i(nx))
  allocate(sim_flam_i(nx))

  ! set flame profile
  call RuntimeParameters_get("flame_initial_position", sim_x0)

  call flame_hse(sim_flam_i, sim_dens_i, sim_temp_i, sim_ye_i, sim_sumy_i, &
                 sim_dens_u, sim_temp_u, fl_eff_ye_u, fl_eff_sumy_u, &
                 sim_dens_b, sim_temp_b, fl_eff_ye_b, fl_eff_sumy_b, &
                 sim_x0, 0.5*deltax, deltax, sim_grav, nx)

  call RuntimeParameters_get("refine_uniform_region", sim_refine_uniform_region)
  call RuntimeParameters_get("refine_region_size", sim_refine_region_size)
  call RuntimeParameters_get("refine_region_stepdown_size", sim_refine_region_stepdown_size)
  call RuntimeParameters_get("refine_lead", sim_refine_lead)
  call RuntimeParameters_get("refine_buf", sim_refine_buf)

  call RuntimeParameters_get("sim_ParticleRefineRegion", sim_ParticleRefineRegion)
  call RuntimeParameters_get("sim_ParticleRefineRegionLevel", sim_ParticleRefineRegionLevel)
  call RuntimeParameters_get("sim_ParticleRefineRegionBottom", sim_ParticleRefineRegionBottom)
  call RuntimeParameters_get("sim_ParticleRefineRegionTop", sim_ParticleRefineRegionTop)

end subroutine Simulation_init

