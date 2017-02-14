!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_addViscousFluxes
!!
!! NAME
!!
!!  hy_uhd_addViscousFluxes
!!
!!
!! SYNOPSIS
!!
!!  hy_uhd_addViscousFluxes(integer(IN) :: blockID,
!!                          integer(IN) :: blkLimitsGC(LOW:HIGH,MDIM),
!!                          integer(IN) :: ix,
!!                          integer(IN) :: iy,
!!                          integer(IN) :: iz,
!!                          real(IN)    :: Flux(HY_VARINUM),
!!                          real(IN)    :: mu(:,:,:),
!!                          integer(IN) :: sweepDir)
!!
!!
!! DESCRIPTION
!!
!!  Adds viscous flux contributions to total (advection) fluxes when using
!!  the explicit, flux-based time integration for heat conduction.
!!
!!
!! ARGUMENTS
!!
!!  blockID     - a local blockID
!!  blkLimitsGC - an array that holds the lower and upper indices of the section 
!!                of block with the guard cells 
!!  ix,iy,iz    - indices of the line along which the sweep is made
!!  Flux        - array containing fluxes
!!  mu          - dynamic viscosity
!!  sweepDir    - direction of sweep
!!
!!***

!!REORDER(4): U

Subroutine hy_uhd_addViscousFluxes(blockID,blkLimitsGC,ix,iy,iz,Flux,mu,sweepDir)

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getDeltas

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Argument List ----------------------------------------------------------
  integer, INTENT(IN) :: blockID,ix,iy,iz
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC 
  real, dimension(HY_VARINUM), intent(INOUT) :: Flux

#ifdef FIXEDBLOCKSIZE 
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, & 
                  GRID_JLO_GC:GRID_JHI_GC, & 
                  GRID_KLO_GC:GRID_KHI_GC),intent(IN) :: mu
#else 
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), & 
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), & 
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),& 
                  intent(IN) :: mu
#endif 
  integer, INTENT(IN) :: sweepDir
  !! ----------------------------------------------------------------------

  real    :: idx,idy,idz,mu_loc
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

  !! Compute resistive parts and add them to flux components
  select case(sweepDir)
  case(DIR_X)
     !! Take a spatial average of visc at each interface 
     mu_loc = 0.5*(mu(ix-1,iy,iz)+mu(ix,iy,iz))

     !! First consider 1D case: d/dy=d/dz=0
     Flux(F02XMOM_FLUX:F04ZMOM_FLUX) = Flux(F02XMOM_FLUX:F04ZMOM_FLUX)&
          -idx*mu_loc*(U(VELX_VAR:VELZ_VAR,ix,iy,iz)-U(VELX_VAR:VELZ_VAR,ix-1,iy,iz))

     Flux(F02XMOM_FLUX) = Flux(F02XMOM_FLUX)-idx*mu_loc*(U(VELX_VAR,ix,iy,iz)-U(VELX_VAR,ix-1,iy,iz))/3.

     Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX)-0.5*idx*mu_loc* &
          ((4./3.)*U(VELX_VAR, ix ,iy,iz)**2+U(VELY_VAR, ix ,iy,iz)**2+U(VELZ_VAR, ix ,iy,iz)**2- &
           (4./3.)*U(VELX_VAR,ix-1,iy,iz)**2-U(VELY_VAR,ix-1,iy,iz)**2-U(VELZ_VAR,ix-1,iy,iz)**2)

#if NDIM >= 2

     !! d/dz=0
     Flux(F02XMOM_FLUX) = Flux(F02XMOM_FLUX)+idy*mu_loc*&
                 (U(VELY_VAR, ix ,iy+1,iz)-U(VELY_VAR, ix ,iy-1,iz)+ &
                  U(VELY_VAR,ix-1,iy+1,iz)-U(VELY_VAR,ix-1,iy-1,iz))/3.

     Flux(F03YMOM_FLUX) = Flux(F03YMOM_FLUX)-&
          0.5*idy*mu_loc*(U(VELX_VAR, ix ,iy+1,iz)-U(VELX_VAR, ix ,iy-1,iz)+ &
                          U(VELX_VAR,ix-1,iy+1,iz)-U(VELX_VAR,ix-1,iy-1,iz))/2.

     Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX)+0.5*idy*mu_loc* &
          ((2./3.)*(U(VELX_VAR, ix ,iy,iz)*(U(VELY_VAR, ix ,iy+1,iz)-U(VELY_VAR, ix ,iy-1,iz))+ &
                    U(VELX_VAR,ix-1,iy,iz)*(U(VELY_VAR,ix-1,iy+1,iz)-U(VELY_VAR,ix-1,iy-1,iz))) &
                  -(U(VELY_VAR, ix ,iy,iz)*(U(VELX_VAR, ix ,iy+1,iz)-U(VELX_VAR, ix ,iy-1,iz))+ &
                    U(VELY_VAR,ix-1,iy,iz)*(U(VELX_VAR,ix-1,iy+1,iz)-U(VELX_VAR,ix-1,iy-1,iz))))/2.

#if NDIM == 3

     Flux(F02XMOM_FLUX) = Flux(F02XMOM_FLUX)+0.5*idz*mu_loc*&
           (2./3.)*(U(VELZ_VAR, ix ,iy,iz+1)-U(VELY_VAR, ix ,iy,iz-1)+ &
                    U(VELZ_VAR,ix-1,iy,iz+1)-U(VELY_VAR,ix-1,iy,iz-1))/2.

     Flux(F04ZMOM_FLUX) = Flux(F04ZMOM_FLUX)-0.5*idz*mu_loc*&
                   (U(VELX_VAR, ix ,iy,iz+1)-U(VELX_VAR, ix ,iy,iz-1)+ &
                    U(VELX_VAR,ix-1,iy,iz+1)-U(VELX_VAR,ix-1,iy,iz-1))/2.

     Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX)+0.5*idz*mu_loc* &
          ((2./3.)*(U(VELX_VAR, ix ,iy,iz)*(U(VELZ_VAR, ix ,iy,iz+1)-U(VELZ_VAR, ix ,iy,iz-1))+ &
                    U(VELX_VAR,ix-1,iy,iz)*(U(VELZ_VAR,ix-1,iy,iz+1)-U(VELZ_VAR,ix-1,iy,iz-1))) &
                  -(U(VELZ_VAR, ix ,iy,iz)*(U(VELX_VAR, ix ,iy,iz+1)-U(VELX_VAR, ix ,iy,iz-1))+ &
                    U(VELZ_VAR,ix-1,iy,iz)*(U(VELX_VAR,ix-1,iy,iz+1)-U(VELX_VAR,ix-1,iy,iz-1))))/2.
#endif
#endif

#if NDIM >= 2
  case(DIR_Y)
     !! Take a spatial average of visc at each interface 
     mu_loc = 0.5*(mu(ix,iy-1,iz)+mu(ix,iy,iz))

     Flux(F02XMOM_FLUX:F04ZMOM_FLUX) = Flux(F02XMOM_FLUX:F04ZMOM_FLUX)&
           -idy*mu_loc*(U(VELX_VAR:VELZ_VAR,ix,iy,iz)-U(VELX_VAR:VELZ_VAR,ix,iy-1,iz))

     Flux(F02XMOM_FLUX) = Flux(F02XMOM_FLUX)-0.5*idx*mu_loc*&
                   (U(VELY_VAR,ix+1, iy ,iz)-U(VELY_VAR,ix-1, iy ,iz)+ &
                    U(VELY_VAR,ix+1,iy-1,iz)-U(VELY_VAR,ix-1,iy-1,iz))/2.

     Flux(F03YMOM_FLUX) = Flux(F03YMOM_FLUX)-&
            idy*mu_loc*(U(VELY_VAR, ix,  iy, iz)-U(VELY_VAR, ix, iy-1,iz))/3.+ &
        0.5*idx*mu_loc*(U(VELX_VAR,ix+1, iy, iz)-U(VELX_VAR,ix-1, iy ,iz)+ &
                        U(VELX_VAR,ix+1,iy-1,iz)-U(VELX_VAR,ix-1,iy-1,iz))/3.

     Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX)&
        +0.5*idx*mu_loc* &
          ((2./3.)*(U(VELY_VAR,ix, iy ,iz)*(U(VELX_VAR,ix+1, iy ,iz)-U(VELX_VAR,ix-1, iy ,iz))+ &
                    U(VELY_VAR,ix,iy-1,iz)*(U(VELX_VAR,ix+1,iy-1,iz)-U(VELX_VAR,ix-1,iy-1,iz))) &
                  -(U(VELX_VAR,ix, iy ,iz)*(U(VELY_VAR,ix+1, iy ,iz)-U(VELY_VAR,ix-1, iy ,iz))+ &
                    U(VELX_VAR,ix,iy-1,iz)*(U(VELY_VAR,ix+1,iy-1,iz)-U(VELY_VAR,ix-1,iy-1,iz))))/2. &
        -0.5*idy*mu_loc* &
                  ( U(VELX_VAR,ix, iy ,iz)**2+(4./3.)*U(VELY_VAR,ix, iy ,iz)**2+U(VELZ_VAR,ix, iy ,iz)**2 &
                   -U(VELX_VAR,ix,iy-1,iz)**2-(4./3.)*U(VELY_VAR,ix,iy-1,iz)**2-U(VELZ_VAR,ix,iy-1,iz)**2)

#if NDIM == 3
     Flux(F03YMOM_FLUX) = Flux(F03YMOM_FLUX)&
        +0.5*idz*mu_loc*&
          (2./3.)* (U(VELZ_VAR,ix, iy, iz+1)-U(VELZ_VAR,ix, iy, iz-1)+ &
                    U(VELZ_VAR,ix,iy-1,iz+1)-U(VELZ_VAR,ix,iy-1,iz-1))/2.

     Flux(F04ZMOM_FLUX) = Flux(F04ZMOM_FLUX)&
        -0.5*idz*mu_loc*&
                   (U(VELY_VAR,ix, iy, iz+1)-U(VELY_VAR,ix, iy, iz-1)+ &
                    U(VELY_VAR,ix,iy-1,iz+1)-U(VELY_VAR,ix,iy-1,iz-1))/2.

     Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX)&
        +0.5*idz*mu_loc* &
          ((2./3.)*(U(VELY_VAR,ix, iy ,iz)*(U(VELZ_VAR,ix, iy ,iz+1)-U(VELZ_VAR,ix, iy ,iz-1))+ &
                    U(VELY_VAR,ix,iy-1,iz)*(U(VELZ_VAR,ix,iy-1,iz+1)-U(VELZ_VAR,ix,iy-1,iz-1))) &
                  -(U(VELZ_VAR,ix, iy ,iz)*(U(VELY_VAR,ix, iy ,iz+1)-U(VELY_VAR,ix, iy ,iz-1))+ &
                    U(VELZ_VAR,ix,iy-1,iz)*(U(VELY_VAR,ix,iy-1,iz+1)-U(VELY_VAR,ix,iy-1,iz-1))))/2.
#endif

#if NDIM == 3
  case(DIR_Z)
     !! Take a spatial average of visc at each interface 
     mu_loc = 0.5*(mu(ix,iy,iz-1)+mu(ix,iy,iz))

     Flux(F02XMOM_FLUX:F04ZMOM_FLUX) = Flux(F02XMOM_FLUX:F04ZMOM_FLUX)&
           -idz*mu_loc*(U(VELX_VAR:VELZ_VAR,ix,iy,iz)-U(VELX_VAR:VELZ_VAR,ix,iy,iz-1))

     Flux(F02XMOM_FLUX) = Flux(F02XMOM_FLUX)&
        -0.5*idx*mu_loc*&
                   (U(VELZ_VAR,ix+1,iy, iz )-U(VELZ_VAR,ix-1,iy, iz )+ &
                    U(VELZ_VAR,ix+1,iy,iz-1)-U(VELZ_VAR,ix-1,iy,iz-1))/2.

     Flux(F03YMOM_FLUX) = Flux(F03YMOM_FLUX)&
        -0.5*idy*mu_loc*&
                   (U(VELZ_VAR,ix,iy+1, iz )-U(VELZ_VAR,ix,iy-1, iz )+ &
                    U(VELZ_VAR,ix,iy+1,iz-1)-U(VELZ_VAR,ix,iy-1,iz-1))/2.

     Flux(F04ZMOM_FLUX) = Flux(F04ZMOM_FLUX)-&
            idz*mu_loc*(U(VELZ_VAR, ix,  iy,  iz )-U(VELZ_VAR, ix,  iy, iz-1))/3.+ &
        0.5*idx*mu_loc*(U(VELX_VAR,ix+1, iy,  iz )-U(VELX_VAR,ix-1, iy,  iz )+ &
                        U(VELX_VAR,ix+1, iy, iz-1)-U(VELX_VAR,ix-1, iy, iz-1))/3.+&
        0.5*idy*mu_loc*(U(VELY_VAR, ix, iy+1, iz )-U(VELY_VAR, ix, iy-1, iz )+ &
                        U(VELY_VAR, ix, iy+1,iz-1)-U(VELY_VAR,ix-1,iy-1,iz-1))/3.

     Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX)&
      -0.5*idx*mu_loc*( U(VELX_VAR,ix,iy, iz )*(U(VELZ_VAR,ix+1,iy, iz )-U(VELZ_VAR,ix-1,iy, iz )) + &
                        U(VELX_VAR,ix,iy,iz-1)*(U(VELZ_VAR,ix+1,iy,iz-1)-U(VELZ_VAR,ix-1,iy,iz-1))   &
                   -(2./3.)*&
                       (U(VELZ_VAR,ix,iy, iz )*(U(VELX_VAR,ix+1,iy, iz )-U(VELX_VAR,ix-1,iy, iz )) + &
                        U(VELZ_VAR,ix,iy,iz-1)*(U(VELX_VAR,ix+1,iy,iz-1)-U(VELX_VAR,ix-1,iy,iz-1))))/2.&
      -0.5*idy*mu_loc*( U(VELY_VAR,ix,iy, iz )*(U(VELZ_VAR,ix,iy+1, iz )-U(VELZ_VAR,ix,iy-1, iz )) + &
                        U(VELY_VAR,ix,iy,iz-1)*(U(VELZ_VAR,ix,iy+1,iz-1)-U(VELZ_VAR,ix,iy-1,iz-1))   &
                   -(2./3.)*&
                       (U(VELZ_VAR,ix,iy, iz )*(U(VELY_VAR,ix,iy+1, iz )-U(VELY_VAR,ix,iy-1, iz )) + &
                        U(VELZ_VAR,ix,iy,iz-1)*(U(VELY_VAR,ix,iy+1,iz-1)-U(VELY_VAR,ix,iy-1,iz-1))))/2.&
      -0.5*idz*mu_loc*( U(VELX_VAR,ix,iy, iz )**2+U(VELY_VAR,ix,iy, iz )**2+(4./3.)*U(VELZ_VAR,ix,iy, iz )**2&
                       -U(VELX_VAR,ix,iy,iz-1)**2-U(VELY_VAR,ix,iy,iz-1)**2-(4./3.)*U(VELZ_VAR,ix,iy,iz-1)**2)


#endif
#endif
  end select


  !! Release pointer
  call Grid_releaseBlkPtr(blockID,U,CENTER)

End Subroutine hy_uhd_addViscousFluxes
