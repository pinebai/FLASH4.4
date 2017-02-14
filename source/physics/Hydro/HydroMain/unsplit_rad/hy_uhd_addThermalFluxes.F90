!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_addThermalFluxes
!!
!! NAME
!!
!!  hy_uhd_addThermalFluxes
!!
!!
!! SYNOPSIS
!!
!!  hy_uhd_addThermalFluxes(integer(IN) :: blockID,
!!                          integer(IN) :: blkLimitsGC(LOW:HIGH,MDIM),
!!                          integer(IN) :: ix,
!!                          integer(IN) :: iy,
!!                          integer(IN) :: iz,
!!                          real(IN)    :: Flux(HY_VARINUM),
!!                          real(IN)    :: kappa(:,:,:),
!!                          integer(IN) :: sweepDir)
!!
!! DESCRIPTION
!!
!!  Adds thermal flux contributions to total (advection) fluxes when using
!!  the explicit, flux-based time integration for heat conduction.
!!
!! ARGUMENTS
!!
!!  blockID     - a local blockID
!!  blkLimitsGC - an array that holds the lower and upper indices of the section 
!!                of block with the guard cells 
!!  ix,iy,iz    - indices of the line along which the sweep is made
!!  Flux        - array containing fluxes
!!  kappa       - thermal conductivity
!!  sweepDir    - direction of sweep
!!
!!***

!!REORDER(4): U

Subroutine hy_uhd_addThermalFluxes(blockID,blkLimitsGC,ix,iy,iz,Flux,kappa,sweepDir)

  use Hydro_data,     ONLY : hy_qref
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getDeltas

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Argument List ----------------------------------------------------------
  integer, INTENT(IN) :: blockID,ix,iy,iz
  integer, dimension(LOW:HIGH,MDIM),intent(IN)  :: blkLimitsGC 
  real,    dimension(HY_VARINUM), intent(INOUT) :: Flux
#ifdef FIXEDBLOCKSIZE 
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, & 
                  GRID_JLO_GC:GRID_JHI_GC, & 
                  GRID_KLO_GC:GRID_KHI_GC),intent(IN) :: kappa
#else 
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), & 
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), & 
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),& 
                  intent(IN) :: kappa
#endif 
  integer, INTENT(IN) :: sweepDir
  !! ----------------------------------------------------------------------

  real    :: idx,idy,idz,kappa_loc
  real, dimension(MDIM) :: del
  real, pointer, dimension(:,:,:,:) :: U


  !! Get deltas
  call Grid_getDeltas(blockID,del)

  idx=1./del(DIR_X)
  if (NDIM >= 2) then
     idy=1./del(DIR_Y)
     if (NDIM == 3) then
        idz=1./del(DIR_Z)
     endif
  endif

  !! Get pointer
  call Grid_getBlkPtr(blockID,U,CENTER)


  select case(sweepDir)
  case(DIR_X)
     !! Take a spatial average of visc at each interface 
     kappa_loc = 0.5*(kappa(ix-1,iy,iz)+kappa(ix,iy,iz))

     Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX)-&
          idx*kappa_loc*(U(TEMP_VAR,ix,iy,iz)-U(TEMP_VAR,ix-1,iy,iz))

#if NDIM >= 2
  case(DIR_Y)
     !! Take a spatial average of visc at each interface 
     kappa_loc = 0.5*(kappa(ix,iy-1,iz)+kappa(ix,iy,iz))

     Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX)-&
          idy*kappa_loc*(U(TEMP_VAR,ix,iy,iz)-U(TEMP_VAR,ix,iy-1,iz))

#if NDIM == 3
  case(DIR_Z)
     !! Take a spatial average of visc at each interface 
     kappa_loc = 0.5*(kappa(ix,iy,iz-1)+kappa(ix,iy,iz))

     Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX)-&
          idz*kappa_loc*(U(TEMP_VAR,ix,iy,iz)-U(TEMP_VAR,ix,iy,iz-1))
#endif
#endif
  end select

  call Grid_releaseBlkPtr(blockID,U,CENTER)

End Subroutine hy_uhd_addThermalFluxes
