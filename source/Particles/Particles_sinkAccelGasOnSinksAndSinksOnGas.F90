!!****f* source/Particles/Particles_sinkAccelGasOnSinksAndSinksOnGas
!!
!! NAME
!!
!!  Particles_sinkAccelGasOnSinksAndSinksOnGas
!!
!! SYNOPSIS
!!
!!  call Particles_sinkAccelGasOnSinksAndSinksOnGas(integer, OPTIONAL, intent(IN) :: accelProps(MDIM),
!!                                                  integer, OPTIONAL ,intent(IN) :: accelVars(MDIM))
!!
!! DESCRIPTION
!!
!!  Computes gas -> sinks and sinks -> gas gravitational accelerations
!!  by direct summation over all sink particles and grid cells.
!!  For cosmology, will also want to get contribution from PDE
!!  (mapped DM delegate particle density).
!!
!! ARGUMENTS
!!
!!  accelProps : optionally give the indices of the sink particle properties
!!               into which the gas-on-sink accelerations should be stored.
!!               Default is ACCX_PART_PROP, ACCY_PART_PROP, ACCZ_PART_PROP.
!!
!!   accelVars : optionally give the indices of the UNK variables
!!               into which the sink-on-gas accelerations should be stored.
!!               Default is SGAX_VAR,SGAY_VAR,SGAZ_VAR.
!!
!! NOTES
!!
!!   written by Christoph Federrath, 2008-2015
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   added optional argument for accel particle properties - Klaus Weide, 2014
!!   merged with Particles_sinkAccelSinksOnGas to speed up computation (Christoph Federrath, 2015)
!!
!! If accelProps is given but contains an invalid index (e.g., 0), the routine
!! returns without updating any sink particle accelerations, but gas solution variables
!! will still be updated for the sink-on-gas accelerations (and pt_sinkGatherGlobal will
!! still have been called).
!!
!!***

subroutine Particles_sinkAccelGasOnSinksAndSinksOnGas(accelProps,accelVars)
  implicit none
#include "constants.h"
  integer, intent(in), OPTIONAL :: accelProps(MDIM)
  integer, intent(in), OPTIONAL :: accelVars(MDIM)

end subroutine Particles_sinkAccelGasOnSinksAndSinksOnGas
