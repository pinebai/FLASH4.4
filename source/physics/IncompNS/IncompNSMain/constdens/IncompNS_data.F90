!!****if* source/physics/IncompNS/IncompNSMain/constdens/IncompNS_data
!!
!! NAME
!!
!!  IncompNS_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE IncompNS_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data and limiter functions that are specific to the Incmp_Navier_Stokes module.
!!
!!***
 
 
module IncompNS_data

#include "Flash.h"
#include "constants.h"
#include "IncompNS.h"


  logical, save :: ins_useIncompNS

  integer, save :: ins_isgs
  integer, save :: ins_cflflg
  integer, save :: ins_intschm
  integer, save :: ins_nstep
  integer, save :: ins_meshMe
  integer, save :: ins_meshNumProcs
  integer, save :: ins_meshComm

  real, save :: ins_cfl
  real, save :: ins_sigma
  real, save :: ins_invRe
  real, save :: ins_dtspec

  integer, parameter :: rkstep = 3

  real, save :: ins_alf(rkstep)
  real, save :: ins_gam(rkstep)
  real, save :: ins_rho(rkstep)
  real, save :: ins_alfa, ins_gama, ins_rhoa

  integer, save :: ins_intschm_type

  logical, save :: ins_prescorr
  real, save :: ins_prescoeff

  logical, save :: ins_restart

  real, save, dimension(LOW:HIGH,MDIM) :: ins_globalDomain
  integer,save,dimension(2,MDIM) :: ins_domainBC

  real, save :: ins_tlevel

  real, save :: ins_Qin
  real, save :: ins_Qout

  logical, save :: ins_predcorrflg

  logical, save :: ins_outflowgridChanged

  real, save :: ins_convvel(LOW:HIGH,MDIM)

  integer, save :: ins_prol_method

  real, save :: ins_meanDivUstar, ins_deltamass

  integer, parameter :: maxblocks_tree = MAXBLOCKS*10 

  ! Arrays for convective outflow BC
  real, save :: uvel_x(NGUARD+1,NYB+2*NGUARD*K2D,NZB+2*NGUARD*K3D,LOW:HIGH,maxblocks_tree)
  real, save :: vvel_x(NGUARD+1,NYB+2*NGUARD*K2D+K2D,NZB+2*NGUARD*K3D,LOW:HIGH,maxblocks_tree)
  real, save :: wvel_x(NGUARD+1,NYB+2*NGUARD*K2D,NZB+2*NGUARD*K3D+K3D,LOW:HIGH,maxblocks_tree)

  real, save :: uvel_y(NGUARD+1,NXB+2*NGUARD+1,NZB+2*NGUARD*K3D,LOW:HIGH,maxblocks_tree)
  real, save :: vvel_y(NGUARD+1,NXB+2*NGUARD,NZB+2*NGUARD*K3D,LOW:HIGH,maxblocks_tree)
  real, save :: wvel_y(NGUARD+1,NXB+2*NGUARD,NZB+2*NGUARD*K3D+K3D,LOW:HIGH,maxblocks_tree)

#if NDIM == 3
  real, save :: uvel_z(NGUARD+1,NXB+2*NGUARD+1,NYB+2*NGUARD*K2D,LOW:HIGH,maxblocks_tree)
  real, save :: vvel_z(NGUARD+1,NXB+2*NGUARD,NYB+2*NGUARD*K2D+1,LOW:HIGH,maxblocks_tree)
  real, save :: wvel_z(NGUARD+1,NXB+2*NGUARD,NYB+2*NGUARD*K2D,LOW:HIGH,maxblocks_tree)
#endif

  ! Container for old timesteps
  real, dimension(rkstep), save :: ins_vardt(-rkstep:0)

  ! Gravitational acceleration:
  real, save    :: ins_gravX, ins_gravY, ins_gravZ 

  ! Fixed pressure gradient, constant mass vars for channel simulations: 
  ! Constant mass Flow in the Z direction.
  real, save    :: ins_dpdx,ins_dpdy,ins_dpdz
  logical, save :: ins_constmass
  real, save    :: ins_WBREF  
  real, save    :: ins_WB, ins_WBold 
  real, save    :: ins_area_solids

end module IncompNS_data
