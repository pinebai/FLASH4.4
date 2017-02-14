!!****if* source/physics/Eos/EosMain/Helmholtz/ExternalAbarZbar/eos_externalComputeAbarZbar
!!
!! NAME
!!
!!  eos_externalComputeAbarZbar
!!
!! SYNOPSIS
!!
!!  call eos_externalComputeAbarZbar(real(in), dimension(:,:)  :: solnscalars,
!!                                   real(out), dimension(:)  :: abardata,
!!                                   real(out), dimension(:)  :: zbardata)
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!! Klaus Weide   2013
!!
!!  This routine private to the Eos unit serves to implement callbacks
!!  for the Eos/Helmholtz/ExternelAbarZbar EOS implementation.
!!  Code units that implement ways for computing Abar and Zbar from
!!  the solutions state (or otherwise) are polled here by calling their
!!  public UNITNAME_computeAbarZbar interfaces.
!!
!!  It is assumed that not more than one polled units have implementations
!!  included in the simulation that actually provide Abar and Zbar.
!!  If several units do provide the data, the last one polled here will win.
!!
!! ARGUMENTS
!!
!!   solnscalars : scalars of the solution 
!!
!!   abardata : abar info
!!
!!   zbardata : zbar info
!!
!!
!!
!!***

#include "Flash.h"
subroutine eos_externalComputeAbarZbar(solnScalars, abarData, zbarData)

  use Flame_interface, ONLY: Flame_computeAbarZbar
  use Burn_interface,  ONLY: Burn_computeAbarZbar

  implicit none

  real, intent(in),  dimension(:,:)  :: solnScalars
  real, intent(out), dimension(:)    :: abarData, zbarData

  call Flame_computeAbarZbar(solnScalars, abarData, zbarData)

#ifdef FLASH_SOURCEBURN
  call Burn_computeAbarZbar(solnScalars, abarData, zbarData)
#endif

end subroutine eos_externalComputeAbarZbar
