!!****if* source/physics/Hydro/HydroMain/split/MHD_8Wave/Hydro_data
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
!!  Placeholder module containing private 8Wave MHD solver-specific data
!!
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

module Hydro_data

#include "Flash.h"

  implicit none

  integer, save :: hy_irenorm
  real,    save :: hy_Rconst
  real,    save :: hy_cfl
  real,    save :: hy_smallpres, hy_smalldens
  logical, save :: hy_killdivb,      &
                   hy_useGravity,    &
                   hy_fluxCorrect,   &
                   hy_useResistivity,&
                   hy_useDiffuse,    &
                   hy_useViscosity,  &
                   hy_useConductivity,&
                   hy_useMagneticResistivity

  real,    save :: hy_maxMagDiff

  ! Energy switch
  real, PARAMETER :: hy_eswitch = 1.e-6

  ! System of units used
  character(4), save :: hy_units

  ! Storage for timestep calculation
  integer,save, DIMENSION(5) :: hy_dtminloc
  real, save :: hy_dtmin

  ! Everybody should know these!
  integer, save :: hy_meshNumProcs, hy_meshMe
  integer, save :: hy_gcMaskSize=NUNK_VARS+NDIM*NFACE_VARS

  ! Constants for non-dimensionalization
  real, save :: hy_xref
  real, save :: hy_tref
  real, save :: hy_dref
  real, save :: hy_vref
  real, save :: hy_pref
  real, save :: hy_eref
  real, save :: hy_qref
  real, save :: hy_bref
  real, save :: hy_gref
  real, save :: hy_mref
  real, save :: hy_nref
  real, save :: hy_kref

  ! Hall parameters
  ! Hall MHD is not supported yet!
  real, save :: hy_hall_parameter, hy_hyperResistivity

end module Hydro_data
