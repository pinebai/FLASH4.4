!!****if* source/Particles/ParticlesMain/active/charged/HybridPIC/pt_picCurl
!!
!! NAME
!!
!!  pt_picCurl
!!
!! SYNOPSIS
!!
!!  call pt_picCurl(integer(in) :: ix,
!!            integer(in) :: iy,
!!            integer(in) :: iz,
!!            integer(in) :: jx,
!!            integer(in) :: jy,
!!            integer(in) :: jz,
!!            real(in) :: mul,
!!            logical(in) :: zer)
!!
!! DESCRIPTION
!!
!!   Compute curl(u(ixyz,:,:,:))*mul and add to u(jxyz,:,:,:)
!!   If zer is true, the target is first zeroed
!!
!! ARGUMENTS
!!
!!   ix :  variable index 
!!
!!   iy : variable index 
!!
!!   iz : variable index 
!!
!!   jx : variable index 
!!
!!   jy : variable index 
!!
!!   jz : variable index 
!!
!!   mul : multiplication factor
!!
!!   zer : whether to zero all variables jx/y/z before processing
!!
!!
!!
!!***

subroutine pt_picCurl(ix, iy, iz, jx, jy, jz, mul, zer)
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr,&
       Grid_getBlkPhysicalSize, Grid_getBlkIndexLimits, Grid_releaseBlkPtr, &
       Grid_fillGuardCells, Grid_getBlkCenterCoords


#include "Flash.h"
#include "constants.h"  
#include "Particles.h"
  implicit none

  integer, intent(in) :: ix, iy, iz, jx, jy, jz
  real,    intent(in) :: mul
  logical, intent(in) :: zer

  integer, save ::  localNumBlocks
  integer :: blockList(MAXBLOCKS), blockCount, localSize(MDIM)
  integer :: blkLimits(LOW:HIGH,MDIM), blkLimitsGC(LOW:HIGH,MDIM), guard(MDIM)
  real :: blockSize(MDIM), blockCenter(MDIM), blockLo(MDIM), h(MDIM)
  integer :: i, j, k, block_no
  real, dimension(:,:,:,:), pointer :: u
  real :: dzdy, dydz, dxdz, dzdx, dydx, dxdy

  call Grid_fillGuardCells( CENTER, ALLDIR)

  call Grid_getListOfBlocks(LEAF, blockList, blockCount)
  do block_no = 1, blockCount
     call Grid_getBlkPtr(blockList(block_no), u)
     ! Get block coordinates
     call Grid_getBlkPhysicalSize(blockList(block_no), blockSize)
     call Grid_getBlkCenterCoords(blockList(block_no), blockCenter)
     blockLo = blockCenter - 0.5*blockSize ! lower block corner
     ! Get cell indices
     call Grid_getBlkIndexLimits(blockList(block_no), &
          blkLimits, blkLimitsGC, CENTER)
     localSize=blkLimits(HIGH,:)-blkLimits(LOW,:)+1   ! NXB, NYB and NZB
     guard = blkLimits(LOW,:)-blkLimitsGC(LOW,:)      ! NGUARD

     h = blockSize/localSize  ! cell size

     if (zer) then 
        u(jx,:,:,:) = 0.0
        u(jy,:,:,:) = 0.0
        u(jz,:,:,:) = 0.0
     end if

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              !
              ! field v for cell (k, j, i) in u(ii,i,j,k)
              ! Compute curl(v)*mul and store on grid in (jx, jy, jz)
              !
              if (NDIM == 3) then
                 dzdy = u(iz,i,j+1,k) - u(iz,i,j-1,k) 
                 dydz = u(iy,i,j,k+1) - u(iy,i,j,k-1) 
                 dxdz = u(ix,i,j,k+1) - u(ix,i,j,k-1) 
                 dzdx = u(iz,i+1,j,k) - u(iz,i-1,j,k) 
                 dydx = u(iy,i+1,j,k) - u(iy,i-1,j,k) 
                 dxdy = u(ix,i,j+1,k) - u(ix,i,j-1,k) 
                 u(jx,i,j,k) = u(jx,i,j,k) + mul*0.5*( &
                      dzdy/h(2) - dydz/h(3) )
                 u(jy,i,j,k) = u(jy,i,j,k) + mul*0.5*( &
                      dxdz/h(3) - dzdx/h(1) )
                 u(jz,i,j,k) = u(jz,i,j,k) + mul*0.5*( &
                      dydx/h(1) - dxdy/h(2) )
              else if (NDIM == 2) then
                 dzdy = u(iz,i,j+1,k) - u(iz,i,j-1,k) 
                 dydz = 0.0
                 dxdz = 0.0
                 dzdx = u(iz,i+1,j,k) - u(iz,i-1,j,k) 
                 dydx = u(iy,i+1,j,k) - u(iy,i-1,j,k) 
                 dxdy = u(ix,i,j+1,k) - u(ix,i,j-1,k) 
                 u(jx,i,j,k) = u(jx,i,j,k) + mul*0.5*dzdy/h(2)
                 u(jy,i,j,k) = u(jy,i,j,k) - mul*0.5*dzdx/h(1) 
                 u(jz,i,j,k) = u(jz,i,j,k) + mul*0.5*( &
                      dydx/h(1) - dxdy/h(2) )
              else    ! 1D
                 dzdy = 0.0
                 dydz = 0.0
                 dxdz = 0.0
                 dzdx = u(iz,i+1,j,k) - u(iz,i-1,j,k) 
                 dydx = u(iy,i+1,j,k) - u(iy,i-1,j,k)
                 dxdy = 0.0
                 u(jy,i,j,k) = u(jy,i,j,k) - mul*0.5*dzdx/h(1)
                 u(jz,i,j,k) = u(jz,i,j,k) + mul*0.5*dydx/h(1)
              end if
           end do
        end do
     end do
     call Grid_releaseBlkPtr(blockList(block_no), u)
  end do

  call Grid_fillGuardCells( CENTER, ALLDIR)

end subroutine pt_picCurl
