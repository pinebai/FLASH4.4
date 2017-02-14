!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_prim2flx
!!
!! NAME
!!
!!  hy_uhd_prim2flx
!!
!! SYNOPSIS
!!
!!  hy_uhd_prim2flx( integer(IN)           :: dir,
!!                   real(IN),dimension(*) :: V(HY_VARINUM2),
!!                   real(OUT)             :: F(HY_VARINUM1))
!!
!! ARGUMENTS
!!
!! dir- directional index
!! V  - primitive variables  + GAMC,GAME
!! F  - flux
!!
!! DESCRIPTION
!!
!!  This routine calculates conversion from primitive variables to fluxes.
!!
!!***

Subroutine hy_uhd_prim2flx(dir,V,F)

#include "Flash.h"
#include "UHD.h"

  use Hydro_data,       ONLY : fP => hy_fPresInMomFlux

#ifdef FLASH_UGLM_MHD
  use Hydro_data, only : hy_C_hyp
#endif

  implicit none

  !! Arguments type declaration -----------
  integer, intent(IN)  :: dir
  real, dimension(*),           intent(IN)  :: V
  real, dimension(HY_VARINUM1), intent(OUT) :: F
  !! --------------------------------------

  real  :: u2,E
  real  :: B2,UB,ptot


  ! initialize with zero
  B2 = 0.
  UB = 0.
  ptot = V(HY_PRES)

  u2 = dot_product(V(HY_VELX:HY_VELZ),V(HY_VELX:HY_VELZ))
  E   = 0.5*V(HY_DENS)*u2 + V(HY_PRES)/(V(HY_GAME)-1.)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  B2 = dot_product(V(HY_MAGX:HY_MAGZ),V(HY_MAGX:HY_MAGZ))
  UB = dot_product(V(HY_VELX:HY_VELZ),V(HY_MAGX:HY_MAGZ))
  ptot= V(HY_PRES) + 0.5 *B2
  E   = E + 0.5*B2
#endif


  F(HY_P_FLUX) = ptot

  if (dir==DIR_X) then
     F(F01DENS_FLUX) = V(HY_DENS)*V(HY_VELX)
#ifdef FLASH_UHD_HYDRO
     F(F02XMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELX) +fP*V(HY_PRES)
     F(F03YMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELY)
     F(F04ZMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELZ)
     F(F05ENER_FLUX) =  (E+V(HY_PRES))*V(HY_VELX)
#else
     F(F02XMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELX)-V(HY_MAGX)*V(HY_MAGX) +fP*ptot
     F(F03YMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELY)-V(HY_MAGX)*V(HY_MAGY)
     F(F04ZMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELZ)-V(HY_MAGX)*V(HY_MAGZ)
     F(F05ENER_FLUX) = (E+ptot)*V(HY_VELX)-V(HY_MAGX)*UB
#ifdef FLASH_USM_MHD
     F(F06MAGX_FLUX) = 0.0
#endif
     F(F07MAGY_FLUX) = V(HY_VELX)*V(HY_MAGY)-V(HY_VELY)*V(HY_MAGX)
     F(F08MAGZ_FLUX) = V(HY_VELX)*V(HY_MAGZ)-V(HY_VELZ)*V(HY_MAGX)
#endif

  elseif (dir==DIR_Y) then
     F(F01DENS_FLUX) = V(HY_DENS)*V(HY_VELY)
#ifdef FLASH_UHD_HYDRO
     F(F02XMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELX)
     F(F03YMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELY) +fP*V(HY_PRES)
     F(F04ZMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELZ)
     F(F05ENER_FLUX) =  (E+V(HY_PRES))*V(HY_VELY)
#else
     F(F02XMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELX)-V(HY_MAGY)*V(HY_MAGX)
     F(F03YMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELY)-V(HY_MAGY)*V(HY_MAGY) +fP*ptot
     F(F04ZMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELZ)-V(HY_MAGY)*V(HY_MAGZ)
     F(F05ENER_FLUX) = (E+ptot)*V(HY_VELY)-V(HY_MAGY)*UB
     F(F06MAGX_FLUX) = V(HY_VELY)*V(HY_MAGX)-V(HY_VELX)*V(HY_MAGY)
#ifdef FLASH_USM_MHD
     F(F07MAGY_FLUX) = 0.0
#endif
     F(F08MAGZ_FLUX) = V(HY_VELY)*V(HY_MAGZ)-V(HY_VELZ)*V(HY_MAGY)
#endif

  elseif (dir==DIR_Z) then
     F(F01DENS_FLUX) = V(HY_DENS)*V(HY_VELZ)
#ifdef FLASH_UHD_HYDRO
     F(F02XMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELX)
     F(F03YMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELY)
     F(F04ZMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELZ) +fP*V(HY_PRES)
     F(F05ENER_FLUX) =  (E+V(HY_PRES))*V(HY_VELZ)
#else
     F(F02XMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELX)-V(HY_MAGZ)*V(HY_MAGX)
     F(F03YMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELY)-V(HY_MAGZ)*V(HY_MAGY)
     F(F04ZMOM_FLUX) = F(F01DENS_FLUX)*V(HY_VELZ)-V(HY_MAGZ)*V(HY_MAGZ) +fP*ptot
     F(F05ENER_FLUX) = (E+ptot)*V(HY_VELZ)-V(HY_MAGZ)*UB
     F(F06MAGX_FLUX) = V(HY_VELZ)*V(HY_MAGX)-V(HY_VELX)*V(HY_MAGZ)
     F(F07MAGY_FLUX) = V(HY_VELZ)*V(HY_MAGY)-V(HY_VELY)*V(HY_MAGZ)
#ifdef FLASH_USM_MHD
     F(F08MAGZ_FLUX) = 0.0
#endif
#endif

  endif



#ifdef FLASH_UGLM_MHD
  if (dir==DIR_X) then
     F(F06MAGX_FLUX) = V(HY_GLMP)
     F(F09GLMP_FLUX) = V(HY_MAGX)*hy_C_hyp*hy_C_hyp
  elseif (dir==DIR_Y) then
     F(F07MAGY_FLUX) = V(HY_GLMP)
     F(F09GLMP_FLUX) = V(HY_MAGY)*hy_C_hyp*hy_C_hyp
  elseif (dir==DIR_Z) then
     F(F08MAGZ_FLUX) = V(HY_GLMP)
     F(F09GLMP_FLUX) = V(HY_MAGZ)*hy_C_hyp*hy_C_hyp
  endif
#endif

End Subroutine hy_uhd_prim2flx
