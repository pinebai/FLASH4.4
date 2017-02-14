!!****if* source/Simulation/SimulationMain/PreHeatPlug/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  Use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data
!!
!!  2014/11/10  Add the material for filled gas
!!  2015/4/29   Add the material for LEH window and Washer (named as Wash)
!!              WashThickness > TargetThickness and WashRadius > TargetRadius are needed
!! 
!!***
module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!  
  real, save :: sim_targetRadius
  real, save :: sim_targetHeight
  real, save :: sim_vacuumHeight
  real, save :: sim_targetThickness         !2014/11/10, by Po-Yu
  real, save :: sim_gasRadius               !2014/11/10, by Po-Yu
  real, save :: sim_windowsThickness        !2014/11/10, by Po-Yu
  real, save :: sim_windowsRadius           !2014/11/10, by Po-Yu  
  real, save :: sim_washThickness         !2015/4/30, by Po-Yu
  real, save :: sim_washRadius            !2015/4/30, by Po-Yu    
  real, save :: sim_plugThickness         !2015/9/12/, by Po-Yu

  real,    save :: sim_rhoTarg  
  real,    save :: sim_teleTarg 
  real,    save :: sim_tionTarg 
  real,    save :: sim_tradTarg 
  real,    save :: sim_zminTarg
  integer, save :: sim_eosTarg

  real,    save :: sim_rhoCham  
  real,    save :: sim_teleCham 
  real,    save :: sim_tionCham 
  real,    save :: sim_tradCham 
  integer, save :: sim_eosCham  
  
  real,    save :: sim_rhoGas               !2014/11/10, by Po-Yu  
  real,    save :: sim_teleGas              !2014/11/10, by Po-Yu
  real,    save :: sim_tionGas              !2014/11/10, by Po-Yu
  real,    save :: sim_tradGas              !2014/11/10, by Po-Yu
  integer, save :: sim_eosGas               !2014/11/10, by Po-Yu
  
  real,    save :: sim_rhoLEH               !2015/4/29, by Po-Yu  
  real,    save :: sim_teleLEH              !2015/4/29, by Po-Yu  
  real,    save :: sim_tionLEH              !2015/4/29, by Po-Yu  
  real,    save :: sim_tradLEH              !2015/4/29, by Po-Yu  
  integer, save :: sim_eosLEH               !2015/4/29, by Po-Yu    
  
  real,    save :: sim_rhoWash               !2015/4/29, by Po-Yu  
  real,    save :: sim_teleWash              !2015/4/29, by Po-Yu  
  real,    save :: sim_tionWash              !2015/4/29, by Po-Yu  
  real,    save :: sim_tradWash              !2015/4/29, by Po-Yu  
  integer, save :: sim_eosWash               !2015/4/29, by Po-Yu    

  real, save :: sim_smallX
  character(len=MAX_STRING_LENGTH), save :: sim_initGeom


end module Simulation_data


