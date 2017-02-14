!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_prim2con
!!
!! NAME
!!
!!  hy_uhd_prim2con
!!
!! SYNOPSIS
!!
!!  hy_uhd_prim2con( real(IN)  :: V(HY_VARINUM2),
!!                   real(OUT) :: CU(HY_VARINUM))
!!
!! ARGUMENTS
!!
!! V  - primitive variables + GAMC,GAME
!! CU - conservative variables 
!!
!! DESCRIPTION
!!
!!  This routine calculates conversions from primitive variables to conservative variables.
!!
!!***

Subroutine hy_uhd_prim2con(V,CU)

  implicit none

#include "Flash.h"
#include "UHD.h"

  !! Arguments type declaration -----------
  real ,dimension(HY_VARINUM2), intent(IN)  :: V
  real ,dimension(HY_VARINUM),  intent(OUT) :: CU
  !! --------------------------------------

  real  :: u2,B2

  u2 = dot_product(V(HY_VELX:HY_VELZ),V(HY_VELX:HY_VELZ))
  B2 = 0.
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  B2 = dot_product(V(HY_MAGX:HY_MAGZ),V(HY_MAGX:HY_MAGZ))
#endif

  CU(HY_DENS) = V(HY_DENS)
  CU(HY_XMOM:HY_ZMOM) = V(HY_DENS)*V(HY_VELX:HY_VELZ)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  CU(HY_MAGX:HY_MAGZ) = V(HY_MAGX:HY_MAGZ)
#endif

#ifdef FLASH_UGLM_MHD
  CU(HY_GLMP) = V(HY_GLMP)      ! GLM is currently unsupported
#endif

  CU(HY_ENER) = 0.5*V(HY_DENS)*u2 + V(HY_PRES)/(V(HY_GAME)-1.) + 0.5*B2

End subroutine hy_uhd_prim2con
