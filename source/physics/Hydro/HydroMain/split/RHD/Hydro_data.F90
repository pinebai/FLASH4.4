!!****if* source/physics/Hydro/HydroMain/split/RHD/Hydro_data
!!
!! NAME
!!   
!!  Hydro_data 
!!
!!
!! SYNOPSIS
!!
!!  use Hydro_data
!!
!!
!! DESCRIPTION
!!
!!  Unit scope data for the relativistic Hydro implementation
!!
!!***

Module Hydro_data

#include "Flash.h"


  ! Energy switch
  real, PARAMETER :: eswitch = 1.e-18

  real,   save :: hy_cfl
  integer,save :: hy_renorm
  integer,save :: hy_reconType

!!$  ! System of units used -- DEV doesn't appear to be ever used
!!$  character(4),save :: hy_units

  integer, save,  DIMENSION(5) :: hy_dtminloc
  real,save :: hy_dtmin

  ! Everybody should know these, right?
  integer, save :: hy_meshMe, hy_meshNumProcs 
 
  ! Constants for non-dimensionalization
  real,save :: hy_xref
  real,save :: hy_tref
  real,save :: hy_dref
  real,save :: hy_vref
  real,save :: hy_pref
  real,save :: hy_eref
  real,save :: hy_gref
  real,save :: hy_nref
  real,save :: hy_gamma  

  integer, save :: hy_meshGeom
  integer, save :: hy_eosMode

  logical, save :: hy_fluxCorrect, hy_useGravity

end module Hydro_data
