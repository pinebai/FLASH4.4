!!****if* source/physics/Hydro/HydroMain/unsplit/MHD_StaggeredMesh/hy_uhd_addOhmicHeating
!!
!! NAME
!!
!!  hy_uhd_addOhmicHeating
!!
!!
!! SYNOPSIS
!!
!!  hy_uhd_addOhmicHeating(blockID,blkLimitsGC,ix,iy,iz,Qohm,eta)
!!
!!  hy_uhd_addOhmicHeating(integer(IN) :: blockID,
!!                         integer(IN) :: blkLimits(LOW:HIGH,MDIM),
!!                         integer(IN) :: ix,
!!                         integer(IN) :: iy,
!!                         integer(IN) :: iz,
!!                         real(INOUT) :: Qohm,
!!                         real(IN)    :: eta)
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
!!  blkLimits   - an array that holds the lower and upper indices of the section 
!!                of block with the guard cells 
!!  ix,iy,iz    - indices of the line along which the sweep is made
!!  Qohm        - heating source for internal energy and electron energy
!!  eta         - magnetic resistivity
!!
!!***

!!REORDER(4): U

Subroutine hy_uhd_addOhmicHeating(blockID,blkLimits,ix,iy,iz,Qohm,eta)

  use Grid_interface, ONLY : Grid_getBlkPtr, &
                             Grid_releaseBlkPtr, &
                             Grid_getDeltas, &
                             Grid_getCellCoords
                             
  use Hydro_data, ONLY: hy_geometry

  
  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Argument List ----------------------------------------------------------
  integer, INTENT(IN) :: blockID,ix,iy,iz
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits
  real, intent(INOUT) :: Qohm
  real, intent(IN)    :: eta

  !! ----------------------------------------------------------------------
  real :: dxBy, dxBz, dyBx, dyBz, dzBx, dzBy
  real    :: Jx,Jy,Jz,idx,idy,idz
  real, dimension(MDIM) :: del
  real, pointer, dimension(:,:,:,:) :: U !,bx,by
  real :: inv_r 
#ifdef FIXEDBLOCKSIZE  
  real, dimension(GRID_ILO:GRID_IHI) :: xCenter  
#else  
  real, dimension(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)) :: xCenter  
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

  !! Get cell x-coords for this block  
  if (hy_geometry /= CARTESIAN) then
     call Grid_getCellCoords(IAXIS,blockID, CENTER,.false.,xCenter, blkLimits(HIGH,IAXIS)- blkLimits(LOW,IAXIS)+1)
  endif

  !! initialize derivatives and currents to zero
  dzBx = 0.0
  dzBy = 0.0
  dyBx = 0.0
  dyBz = 0.0
  dxBy = 0.0
  dxBz = 0.0
  Jx   = 0.0
  Jy   = 0.0
  Jz   = 0.0
    
  if (hy_geometry == CARTESIAN) then

  !! Compute partial derivatives to construct J
  
     !! 1D case : d/dy=d/dz=0
     dxBz = (U(MAGZ_VAR,ix+1,iy,iz)-U(MAGZ_VAR,ix-1,iy,iz))*idx*0.5
     dxBy =  (U(MAGY_VAR,ix+1,iy,iz)-U(MAGY_VAR,ix-1,iy,iz))*idx*0.5

#if NDIM >= 2

     !! 2D case : d/dy .ne. 0 but d/dz=0
     dyBx = (U(MAGX_VAR,ix,  iy+1,iz) - U(MAGX_VAR,ix,  iy-1,iz))*0.5*idy
     dyBz = (U(MAGZ_VAR,ix,  iy+1,iz) - U(MAGZ_VAR,ix,  iy-1,iz))*0.5*idy 

#if NDIM == 3

     dzBx = (U(MAGX_VAR, ix,  iy,iz+1) - U(MAGX_VAR, ix,  iy, iz-1))*0.5*idz
     dzBy = (U(MAGY_VAR, ix,  iy,iz+1) - U(MAGY_VAR, ix,  iy, iz-1))*0.5*idz

#endif
#endif
  endif
  
  if (hy_geometry == CYLINDRICAL) then

  !! Notice that X == R, Y == Z, Z == PHI. Be aware of signs
  !! when calculating curls    
     !! 1D case : d/dy=d/dz=0
     inv_r = 1.0/xCenter(ix)
     dxBy = (U(MAGY_VAR,ix+1,iy,iz) - U(MAGY_VAR,ix-1,iy,iz))*idx*0.5
     dxBz = (U(MAGZ_VAR,ix+1,iy,iz) - U(MAGZ_VAR,ix-1,iy,iz))*idx*0.5 &
          +  U(MAGZ_VAR,ix,iy,iz)*inv_r
          
#if NDIM >= 2

     !! 2D case : d/dy .ne. 0 but d/dz=0
     dyBx = (U(MAGX_VAR,ix,  iy+1,iz) - U(MAGX_VAR,ix,  iy-1,iz))*idy*0.5
     dyBz = (U(MAGZ_VAR,ix,  iy+1,iz) - U(MAGZ_VAR,ix,  iy-1,iz))*idy*0.5

#endif
  endif
  
  !! Get Jx, Jy, Jz
  Jx = dyBz - dzBy
  Jy = dzBx - dxBz
  Jz = dxBy - dyBx
      
  !! Get Qohm = eta J^2   
  Qohm = eta*(Jx*Jx + Jy*Jy +Jz*Jz)
  
  !! Release pointer
  call Grid_releaseBlkPtr(blockID,U,CENTER)

End Subroutine hy_uhd_addOhmicHeating
