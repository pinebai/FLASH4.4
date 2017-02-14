!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_unitConvert
!!
!! NAME
!!
!!  hy_uhd_unitConvert
!!
!! SYNOPSIS
!!
!!  hy_uhd_unitConvert( integer (IN) :: blockID, 
!!                      integer (IN) :: convertDir)                    
!!
!! DESCRIPTION
!!
!!  This routine converts physical unit of measure (CGS, SI, or NONE) 
!!  for the unsplit MHD/Hydro units.
!!
!! ARGUMENTS
!!
!!  blockID    - local block ID
!!  convertDir - a direction of conversion
!!
!!***

!!REORDER(4): U,B[xyz]

Subroutine hy_uhd_unitConvert(blockID,convertDir)

  use Hydro_data,     ONLY : hy_dref, hy_eref, hy_pref, &
                             hy_vref, hy_bref
  use Grid_interface, ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Argument list -------------------------
  integer, intent(IN) :: blockID, convertDir
  !! ---------------------------------------

  real, pointer, dimension(:,:,:,:) :: U

#ifdef FLASH_USM_MHD
#if NFACE_VARS > 0
#if NDIM > 1
  real, pointer, dimension(:,:,:,:) :: Bx,By,Bz
#endif
#endif
#endif

  call Grid_getBlkPtr(blockID,U,CENTER)

#ifdef FLASH_USM_MHD
#if NFACE_VARS > 0
#if NDIM > 1
  call Grid_getBlkPtr(blockID,Bx,FACEX)
  call Grid_getBlkPtr(blockID,By,FACEY)
  if (NDIM == 3) call Grid_getBlkPtr(blockID,Bz,FACEZ)
#endif
#endif
#endif

  select case(convertDir)
  case(FWDCONVERT)
     U(DENS_VAR,:,:,:) = U(DENS_VAR,:,:,:)/hy_dref
     U(ENER_VAR,:,:,:) = U(ENER_VAR,:,:,:)/hy_eref
     U(PRES_VAR,:,:,:) = U(PRES_VAR,:,:,:)/hy_pref
     U(VELX_VAR:VELZ_VAR,:,:,:) = U(VELX_VAR:VELZ_VAR,:,:,:)/hy_vref

#ifdef FLASH_USM_MHD
     U(MAGX_VAR:MAGZ_VAR,:,:,:) = U(MAGX_VAR:MAGZ_VAR,:,:,:)/hy_bref
#if NFACE_VARS > 0
#if NDIM > 1
     if (NDIM >= 2) then
        Bx(MAG_FACE_VAR,:,:,:) = Bx(MAG_FACE_VAR,:,:,:)/hy_bref
        By(MAG_FACE_VAR,:,:,:) = By(MAG_FACE_VAR,:,:,:)/hy_bref
        if (NDIM==3) Bz(MAG_FACE_VAR,:,:,:) = Bz(MAG_FACE_VAR,:,:,:)/hy_bref
     endif
#endif
#endif
#endif

  case(BWDCONVERT)
     U(DENS_VAR,:,:,:) = U(DENS_VAR,:,:,:)*hy_dref
     U(ENER_VAR,:,:,:) = U(ENER_VAR,:,:,:)*hy_eref
     U(VELX_VAR:VELZ_VAR,:,:,:) = U(VELX_VAR:VELZ_VAR,:,:,:)*hy_vref

#ifdef FLASH_USM_MHD
     U(MAGX_VAR:MAGZ_VAR,:,:,:) = U(MAGX_VAR:MAGZ_VAR,:,:,:)*hy_bref
#if NFACE_VARS > 0
#if NDIM > 1
     if (NDIM >= 2) then
        Bx(MAG_FACE_VAR,:,:,:) = Bx(MAG_FACE_VAR,:,:,:)*hy_bref
        By(MAG_FACE_VAR,:,:,:) = By(MAG_FACE_VAR,:,:,:)*hy_bref
        if (NDIM==3) Bz(MAG_FACE_VAR,:,:,:) = Bz(MAG_FACE_VAR,:,:,:)*hy_bref
     endif
#endif
#endif
#endif

  end select

  call Grid_releaseBlkPtr(blockID,U,CENTER)
#ifdef FLASH_USM_MHD
#if NFACE_VARS > 0
#if NDIM > 1
  call Grid_releaseBlkPtr(blockID,Bx,FACEX)
  call Grid_releaseBlkPtr(blockID,By,FACEY)
  if (NDIM==3) call Grid_releaseBlkPtr(blockID,Bz,FACEZ)
#endif
#endif
#endif

End Subroutine hy_uhd_unitConvert

  
