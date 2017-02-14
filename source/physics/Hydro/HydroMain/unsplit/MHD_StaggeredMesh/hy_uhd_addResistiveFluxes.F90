!!****if* source/physics/Hydro/HydroMain/unsplit/MHD_StaggeredMesh/hy_uhd_addResistiveFluxes
!!
!! NAME
!!
!!  hy_uhd_addResistiveFluxes
!!
!!
!! SYNOPSIS
!!
!!  hy_uhd_addResistiveFluxes(blockID,blkLimitsGC,ix,iy,iz,Flux,magVisc,sweepDir)
!!
!!  hy_uhd_addResistiveFluxes(integer(IN) :: blockID,
!!                            integer(IN) :: blkLimitsGC(LOW:HIGH,MDIM),
!!                            integer(IN) :: ix,
!!                            integer(IN) :: iy,
!!                            integer(IN) :: iz,
!!                            real(INOUT)    :: Flux,
!!                            real(IN)    :: eta,
!!                            integer(IN) :: sweepDir)
!!
!!
!! DESCRIPTION
!!
!!  Adds resistive flux contributions to total MHD fluxes
!!
!!
!! ARGUMENTS
!!
!!  blockID     - a local blockID
!!  blkLimitsGC - an array that holds the lower and upper indices of the section 
!!                of block with the guard cells 
!!  ix,iy,iz    - indices of the line along which the sweep is made
!!  Flux        - array containing MHD fluxes
!!  eta         - magnetic resistivity
!!  sweepDir    - direction of sweep
!!
!!***

!!REORDER(4): U

Subroutine hy_uhd_addResistiveFluxes(blockID,blkLimitsGC,ix,iy,iz,Flux,eta,sweepDir)

  use Grid_interface, ONLY : Grid_getBlkPtr, &
                             Grid_releaseBlkPtr, &
                             Grid_getDeltas, &
                             Grid_getCellCoords
                             
  use hy_uhd_slopeLimiters, ONLY : mc, &
                                   minmod, &
                                   vanLeer, &
                                   signum, &
                                   get_upwind

  use Hydro_data, ONLY: hy_geometry, hy_conserveAngField

  
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
                  GRID_KLO_GC:GRID_KHI_GC),intent(IN) :: eta
#else 
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), & 
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), & 
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),& 
                  intent(IN) :: eta
#endif 
  integer, INTENT(IN) :: sweepDir
  !! ----------------------------------------------------------------------
  real :: dxBy, dxBz, dyBx, dyBz, dzBx, dzBy
  real    :: jx,jy,jz,idx,idy,idz,eta_loc
  real, dimension(MDIM) :: del
  real, pointer, dimension(:,:,:,:) :: U !,bx,by
  real :: inv_dVr !2/(rp*abs(rp) - rm*abs(rm))
#ifdef FIXEDBLOCKSIZE  
  real, dimension(GRID_IHI_GC) :: xCenter  
#else  
  real, dimension(blkLimitsGC(HIGH,IAXIS)) :: xCenter  
#endif   


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
!!$  call Grid_getBlkPtr(blockID,bx,FACEX)
!!$  call Grid_getBlkPtr(blockID,by,FACEY)

  !! Get cell x-coords for this block  
  if (hy_geometry /= CARTESIAN) then
     call Grid_getCellCoords(IAXIS,blockID, CENTER,.true.,xCenter, blkLimitsGC(HIGH,IAXIS))
  endif

  if (hy_geometry == CARTESIAN) then

  !! Compute resistive parts and add them to flux components
     select case(sweepDir)
     case(DIR_X)
     !! Take an upwinding at each interface 
!     eta_loc = get_upwind(FLUX(F01DENS_FLUX),eta(ix-1,iy,iz),eta(ix,iy,iz))
        eta_loc = (eta(ix-1,iy,iz) + eta(ix,iy,iz))*0.5
      
     !! 1D case : d/dy=d/dz=0
     !! jy = -dBz/dx
!     jy = -minmod(U(MAGZ_VAR,ix+1,iy,iz)-U(MAGZ_VAR,ix-1,iy,iz),&
!                  U(MAGZ_VAR,ix,  iy,iz)-U(MAGZ_VAR,ix-2,iy,iz))*0.5*idx
        jy = -(U(MAGZ_VAR,ix,iy,iz)-U(MAGZ_VAR,ix-1,iy,iz))*idx
     !! jz = dBy/dx
!     jz =  minmod(U(MAGY_VAR,ix+1,iy,iz)-U(MAGY_VAR,ix-1,iy,iz),&
!                  U(MAGY_VAR,ix,  iy,iz)-U(MAGY_VAR,ix-2,iy,iz))*0.5*idx
        jz =  (U(MAGY_VAR,ix,iy,iz)-U(MAGY_VAR,ix-1,iy,iz))*idx

#if NDIM >= 2
     !! 2D case : d/dy .ne. 0 but d/dz=0
     !! jz = dBy/dx - dBx/dy
!     jz = jz - minmod(U(MAGX_VAR,ix,  iy+1,iz)-U(MAGX_VAR,ix,  iy,  iz)&
!                     +U(MAGX_VAR,ix-1,iy+1,iz)-U(MAGX_VAR,ix-1,iy,  iz),&
!                      U(MAGX_VAR,ix,  iy,  iz)-U(MAGX_VAR,ix,  iy-1,iz)&
!                     +U(MAGX_VAR,ix-1,iy,  iz)-U(MAGX_VAR,ix-1,iy-1,iz))*0.5*idy
        jz = jz - (U(MAGX_VAR,ix,  iy+1,iz) - U(MAGX_VAR,ix,  iy-1,iz) &
                +  U(MAGX_VAR,ix-1,iy+1,iz) - U(MAGX_VAR,ix-1,iy-1,iz))*0.25*idy 
#if NDIM == 3
     !! jy = dBx/dz - dBz/dx
!     jy = jy + minmod(U(MAGX_VAR, ix,  iy, iz+1)-U(MAGX_VAR, ix,  iy, iz  )&
!                     +U(MAGX_VAR, ix-1,iy, iz+1)-U(MAGX_VAR, ix-1,iy, iz  ),&
!                      U(MAGX_VAR, ix,  iy, iz  )-U(MAGX_VAR, ix,  iy, iz-1)&
!                     +U(MAGX_VAR, ix-1,iy, iz  )-U(MAGX_VAR, ix-1,iy, iz-1))*0.5*idz

        jy = jy + (U(MAGX_VAR, ix,  iy,iz+1) - U(MAGX_VAR, ix,  iy, iz-1)&
                +  U(MAGX_VAR, ix-1,iy,iz+1) - U(MAGX_VAR, ix-1,iy, iz-1))*0.25*idz


#endif
#endif

        Flux(F07MAGY_FLUX) = Flux(F07MAGY_FLUX) - eta_loc*jz
        Flux(F08MAGZ_FLUX) = Flux(F08MAGZ_FLUX) + eta_loc*jy

     !! eta*(By*jz - Bz*jy)
     !! Use upwinding
!     Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX) &
!          -eta_loc*( get_upwind(Flux(F01DENS_FLUX),U(MAGY_VAR,ix-1,iy,iz),U(MAGY_VAR,ix,iy,iz))*jz &
!                    -get_upwind(Flux(F01DENS_FLUX),U(MAGZ_VAR,ix-1,iy,iz),U(MAGZ_VAR,ix,iy,iz))*jy )

        Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX) &
             -eta_loc*( (U(MAGY_VAR,ix-1,iy,iz)+U(MAGY_VAR,ix,iy,iz))*jz*0.5 &
                      - (U(MAGZ_VAR,ix-1,iy,iz)+U(MAGZ_VAR,ix,iy,iz))*jy*0.5 )


#if NDIM >= 2
     case(DIR_Y)
     !! Take an upwinding
!     eta_loc = get_upwind(Flux(F01DENS_FLUX),eta(ix,iy-1,iz),eta(ix,iy,iz))
        eta_loc = (eta(ix,iy,iz) + eta(ix,iy-1,iz))*0.5

     !! jx = dBz/dy
!     jx =  minmod(U(MAGZ_VAR,ix,   iy+1,  iz )-U(MAGZ_VAR,ix,  iy-1, iz),&
!                  U(MAGZ_VAR,ix,   iy,    iz )-U(MAGZ_VAR,ix,  iy-2, iz))*0.5*idy
        jx = (U(MAGZ_VAR,ix,iy,iz) - U(MAGZ_VAR,ix,iy-1,iz))*idy

     !! jz = dBy/dx - dBx/dy
!     jz = -minmod(U(MAGX_VAR,ix,  iy+1,iz )-U(MAGX_VAR,ix,  iy-1,iz ),&
!                  U(MAGX_VAR,ix,  iy,  iz )-U(MAGX_VAR,ix,  iy-2,iz ))*0.5*idy &
!          +minmod(U(MAGY_VAR,ix+1,iy , iz )-U(MAGY_VAR,ix,  iy,  iz )  &
!                 +U(MAGY_VAR,ix+1,iy-1,iz )-U(MAGY_VAR,ix,  iy-1,iz ), &
!                  U(MAGY_VAR,ix,  iy,  iz )-U(MAGY_VAR,ix-1,iy,  iz )  &
!                 +U(MAGY_VAR,ix,  iy-1,iz )-U(MAGY_VAR,ix-1,iy-1,iz ))*0.5*idx
        jz = - (U(MAGX_VAR,ix  ,iy  ,iz ) - U(MAGX_VAR,ix  ,iy-1,iz))*idy &
             + (U(MAGY_VAR,ix+1,iy  ,iz ) - U(MAGY_VAR,ix-1,iy  ,iz)      & 
             +  U(MAGY_VAR,ix+1,iy-1,iz ) - U(MAGY_VAR,ix-1,iy-1,iz))*0.25*idx



#if NDIM == 3
     !! jx = dBz/dy - dBy/dz
!     jx = jx - minmod(U(MAGY_VAR,ix, iy,  iz+1 )-U(MAGY_VAR,ix, iy,  iz  ) &
!                     +U(MAGY_VAR,ix, iy-1,iz+1 )-U(MAGY_VAR,ix, iy-1,iz  ),&
!                      U(MAGY_VAR,ix, iy,  iz   )-U(MAGY_VAR,ix, iy,  iz-1) &
!                     +U(MAGY_VAR,ix, iy-1,iz   )-U(MAGY_VAR,ix, iy-1,iz-1))*0.5*idz

        jx = jx - (U(MAGY_VAR,ix,iy  ,iz+1 ) - U(MAGY_VAR,ix,iy  ,iz-1)        & 
              +    U(MAGY_VAR,ix,iy-1,iz+1 ) - U(MAGY_VAR,ix,iy-1,iz-1) )*0.25*idx


#endif

        Flux(F06MAGX_FLUX) = Flux(F06MAGX_FLUX) + eta_loc*jz
        Flux(F08MAGZ_FLUX) = Flux(F08MAGZ_FLUX) - eta_loc*jx

     !! eta*(Bz*jx - Bx*jz)
     !! Use upwinding
!     Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX) &
!          -0.5*eta_loc*( get_upwind(Flux(F01DENS_FLUX),U(MAGZ_VAR,ix,iy-1,iz),U(MAGZ_VAR,ix,iy,iz))*jx &
!                        -get_upwind(Flux(F01DENS_FLUX),U(MAGX_VAR,ix,iy-1,iz),U(MAGX_VAR,ix,iy,iz))*jz )

        Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX) &
             -eta_loc*( (U(MAGZ_VAR,ix,iy-1,iz) + U(MAGZ_VAR,ix,iy,iz))*jx*0.5 &
                      - (U(MAGX_VAR,ix,iy-1,iz) + U(MAGX_VAR,ix,iy,iz))*jz*0.5 )

#if NDIM == 3
     case(DIR_Z)
     !! Take an upwinding
!     eta_loc = get_upwind(Flux(F01DENS_FLUX),eta(ix,iy,iz-1),eta(ix,iy,iz))
        eta_loc = (eta(ix,iy,iz) + eta(ix,iy,iz-1))*0.5

     !! jx = dBz/dy - dBy/dz
!     jx = minmod(U(MAGZ_VAR,ix, iy+1,iz  )-U(MAGZ_VAR,ix, iy,  iz  ) &
!                +U(MAGZ_VAR,ix, iy+1,iz-1)-U(MAGZ_VAR,ix, iy,  iz-1),&
!                 U(MAGZ_VAR,ix, iy,  iz  )-U(MAGZ_VAR,ix, iy-1,iz  ) &
!                +U(MAGZ_VAR,ix, iy,  iz-1)-U(MAGZ_VAR,ix, iy-1,iz-1))*0.5*idy&
!         -minmod(U(MAGY_VAR,ix, iy,  iz+1)-U(MAGY_VAR,ix,  iy,  iz-1),&
!                 U(MAGY_VAR,ix, iy,  iz  )-U(MAGY_VAR,ix,  iy,  iz-2))*0.5*idz

        jx = (U(MAGZ_VAR,ix,iy+1,iz  ) -U(MAGZ_VAR,ix,iy-1,iz  )           &
            + U(MAGZ_VAR,ix,iy+1,iz-1) -U(MAGZ_VAR,ix,iy-1,iz-1))*0.25*idy &
            -(U(MAGY_VAR,ix,iy  ,iz  ) -U(MAGY_VAR,ix,iy  ,iz-1))*idz


     !! jy = dBx/dz - dBz/dx
!     jy = minmod(U(MAGX_VAR,ix,  iy, iz+1)-U(MAGX_VAR,ix,  iy, iz-1),&
!                 U(MAGX_VAR,ix,  iy, iz  )-U(MAGX_VAR,ix,  iy, iz-2))*0.5*idz&
!         -minmod(U(MAGZ_VAR,ix+1,iy, iz  )-U(MAGZ_VAR,ix,  iy, iz ) &
!                +U(MAGZ_VAR,ix+1,iy, iz-1)-U(MAGZ_VAR,ix,  iy, iz-1),&
!                 U(MAGZ_VAR,ix,  iy, iz  )-U(MAGZ_VAR,ix-1,iy, iz ) &
!                +U(MAGZ_VAR,ix,  iy, iz-1)-U(MAGZ_VAR,ix-1,iy, iz-1))*0.5*idx

        jy = (U(MAGX_VAR,ix  ,iy,iz  )-U(MAGX_VAR,ix  ,iy,iz-1))*idz &
            -(U(MAGZ_VAR,ix+1,iy,iz  )-U(MAGZ_VAR,ix-1,iy,iz  )      &
            + U(MAGZ_VAR,ix+1,iy,iz-1)-U(MAGZ_VAR,ix-1,iy,iz-1))*0.25*idx

        Flux(F06MAGX_FLUX) = Flux(F06MAGX_FLUX) - eta_loc*jy
        Flux(F07MAGY_FLUX) = Flux(F07MAGY_FLUX) + eta_loc*jx

     !! eta*(Bx*jy - By*jx)
!     Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX) &
!          -0.5*eta_loc*( get_upwind(Flux(F01DENS_FLUX),U(MAGX_VAR,ix,iy,iz-1),U(MAGX_VAR,ix,iy,iz))*jy &
!                        -get_upwind(Flux(F01DENS_FLUX),U(MAGY_VAR,ix,iy,iz-1),U(MAGY_VAR,ix,iy,iz))*jx )


        Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX) &
             -eta_loc*( (U(MAGX_VAR,ix,iy,iz-1) + U(MAGX_VAR,ix,iy,iz))*jy*0.5 &
                      - (U(MAGY_VAR,ix,iy,iz-1) + U(MAGY_VAR,ix,iy,iz))*jx*0.5 )

#endif
#endif
     end select
  endif

  if (hy_geometry == CYLINDRICAL) then

  !! Compute resistive parts and add them to flux components
  !! Notice that X == R, Y == Z, Z == PHI. Be aware of signs
  !! when calculating curls
  
     select case(sweepDir)
     case(DIR_X) ! R-direction
        eta_loc = (eta(ix-1,iy,iz) + eta(ix,iy,iz))*0.5
      
     !! 1D case : d/dy=d/dz=0
     !! jy = Jz = 1/r d/dr (rBphi)

  !      dxBz = (xCenter(ix)*U(MAGZ_VAR,ix,iy,iz)&
  !      -xCenter(ix-1)*U(MAGZ_VAR,ix-1,iy,iz))*idx/(0.5*(xCenter(ix)+xCenter(ix-1)))
        
     !! inv_dVr = 2/(rp*abs(rp) - rm*abs(rm))
       
!        if (xCenter(ix-1)<0.0 .and. xCenter(ix) > 0.0) then
        inv_dVr = xCenter(ix)*xCenter(ix) - xCenter(ix-1)*abs(xCenter(ix-1))
        inv_dVr = 2.0/inv_dVr
        
        dxBz = (U(MAGZ_VAR,ix  ,iy,iz)*xCenter(ix) &
             -  U(MAGZ_VAR,ix-1,iy,iz)*abs(xCenter(ix-1)))*inv_dVr
 !       endif
        
        dxBy = (U(MAGY_VAR,ix,iy,iz)-U(MAGY_VAR,ix-1,iy,iz))*idx
      
#if NDIM >= 2
     !! 2D case : d/dy .ne. 0 but d/dz=0
          dyBx =  (U(MAGX_VAR,ix,  iy+1,iz) - U(MAGX_VAR,ix,  iy-1,iz) &
                +  U(MAGX_VAR,ix-1,iy+1,iz) - U(MAGX_VAR,ix-1,iy-1,iz))*0.25*idy 
          dzBx = 0.0
#endif
        jy = dzBx - dxBz
        jz = dxBy - dyBx 

        Flux(F07MAGY_FLUX) = Flux(F07MAGY_FLUX) - eta_loc*jz
        if (hy_conserveAngField) then 
           Flux(F08MAGZ_FLUX) = Flux(F08MAGZ_FLUX) + eta_loc*jy !!this is added as a source if not conserveAngField 
        endif 
     !! eta*(By*jz - Bz*jy)

        Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX) &
            -eta_loc*(  (U(MAGY_VAR,ix-1,iy,iz)+U(MAGY_VAR,ix,iy,iz))*jz*0.5 &
                      - (U(MAGZ_VAR,ix-1,iy,iz)+U(MAGZ_VAR,ix,iy,iz))*jy*0.5 )

#if NDIM >= 2
     case(DIR_Y)
        eta_loc = (eta(ix,iy,iz) + eta(ix,iy-1,iz))*0.5

        dyBx = (U(MAGX_VAR,ix,iy,iz) - U(MAGX_VAR,ix,iy-1,iz))*idy
        dyBz = (U(MAGZ_VAR,ix,iy,iz) - U(MAGZ_VAR,ix,iy-1,iz))*idy
        dxBy = (U(MAGY_VAR,ix+1,iy,iz)   - U(MAGY_VAR,ix-1,iy,iz)      &  
             +  U(MAGY_VAR,ix+1,iy-1,iz) - U(MAGY_VAR,ix-1,iy-1,iz))*0.25*idx
        dzBy = 0.0 
   
        jx = dyBz - dzBy
        jz = dxBy - dyBx 

        Flux(F06MAGX_FLUX) = Flux(F06MAGX_FLUX) + eta_loc*jz
        Flux(F08MAGZ_FLUX) = Flux(F08MAGZ_FLUX) - eta_loc*jx

     !! eta*(Bz*jx - Bx*jz)

        Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX) &
             -eta_loc*( (U(MAGZ_VAR,ix,iy-1,iz) + U(MAGZ_VAR,ix,iy,iz))*jx*0.5 &
                      - (U(MAGX_VAR,ix,iy-1,iz) + U(MAGX_VAR,ix,iy,iz))*jz*0.5 )
#endif
     end select
  endif
  

  !! Release pointer
  call Grid_releaseBlkPtr(blockID,U,CENTER)
!!$  call Grid_releaseBlkPtr(blockID,bx,FACEX)
!!$  call Grid_releaseBlkPtr(blockID,by,FACEY)

End Subroutine hy_uhd_addResistiveFluxes
