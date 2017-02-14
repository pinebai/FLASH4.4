!!****if* source/Particles/ParticlesMain/active/charged/HybridPIC/pt_picData
!!
!! NAME
!!    pt_picData
!!
!! SYNOPSIS
!!    pt_picData
!!
!! DESCRIPTION
!!    Module to hold local variables and data types for the PIC implementation
!!
!!***

module pt_picData

  implicit none

#include "constants.h"
  
  character(len=MAX_STRING_LENGTH) :: pt_picPname_1
  real, save    :: pt_picPmass_1, pt_picPcharge_1
  real, save    :: pt_picPdensity_1
  real, save    :: pt_picPtemp_1
  real, save    :: pt_picPvelx_1, pt_picPvely_1, pt_picPvelz_1
  integer, save :: pt_picPpc_1

  character(len=MAX_STRING_LENGTH) :: pt_picPname_2
  real, save    :: pt_picPmass_2, pt_picPcharge_2
  real, save    :: pt_picPdensity_2
  real, save    :: pt_picPtemp_2
  real, save    :: pt_picPvelx_2, pt_picPvely_2, pt_picPvelz_2
  integer, save :: pt_picPpc_2

  real, save    :: pt_picTe
  real, save    :: pt_picResistivity, pt_picResistivityHyper
  real, save    :: pt_picGam

  real, save    :: pt_picCdensMin
  integer, save :: pt_picNsub

  integer, save :: pt_picRng_seed
  
  integer, dimension(LOW:HIGH, MDIM), save :: pt_picDomainBC
  real, dimension(LOW:HIGH, MDIM), save :: pt_picDomainBoundBox

end module pt_picData
