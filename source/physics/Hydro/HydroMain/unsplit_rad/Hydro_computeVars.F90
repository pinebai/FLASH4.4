!!****if* source/physics/Hydro/HydroMain/unsplit_rad/Hydro_computeVars
!!
!! NAME
!!
!!  Hydro_computeVars
!!
!!
!! SYNOPSIS
!!
!!  Hydro_computeVars()
!!
!! DESCRIPTION
!!
!!  This routine computes any variable user wish to compute anywhere in the code,
!!  not only within the hydro/MHD scopes, but also from Driver, or elsewhere.
!!  An example shown here is to compute the divergence of mangnetic fields for MHD.
!!
!! ARGUMENTS
!!
!!***

!!REORDER(4): U
Subroutine Hydro_computeVars()  

#include "Flash.h"
#include "constants.h"

  use Hydro_data
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,&
                               Grid_getCellCoords,Grid_getBlkIndexLimits,&
                               Grid_getDeltas,Grid_getListOfBlocks

  implicit none

  real, pointer :: U(:,:,:,:),Bx(:,:,:,:),By(:,:,:,:),Bz(:,:,:,:)
  integer :: i, j, k
  real    :: delxinv, delyinv, delzinv
  real    :: rc, rp, rm
  real, dimension(MDIM) :: del
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  real, allocatable, dimension(:) :: xCenter,yCenter,zCenter
  integer, dimension(MAXBLOCKS) :: blockList
  integer :: in,blockID, numLeafBlocks


  call Grid_getListOfBlocks(LEAF,blockList, numLeafBlocks)

#if NDIM > 1 /* compute divB only for 2D or 3D */
#if NFACE_VARS > 0 /* compute divB when facevars are defined */
  do in = 1, numLeafBlocks

     blockID = blockList(in)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getDeltas(blockID,del)
     call Grid_getBlkPtr(blockID,U,CENTER)
     call Grid_getBlkPtr(blockID,Bx,FACEX)
     call Grid_getBlkPtr(blockID,By,FACEY)
     if (NDIM == 3) &
     call Grid_getBlkPtr(blockID,Bz,FACEZ)

     allocate(xCenter(blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1))
     allocate(yCenter(blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1))
     if (NDIM == 3) &
     allocate(zCenter(blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1))


     if (hy_geometry /= CARTESIAN) then
        !get coord info will use this for Areas and Volumes
        call Grid_getCellCoords(IAXIS,blockID, CENTER,.true.,xCenter, blkLimitsGC(HIGH,IAXIS))
        if (NDIM > 1) call Grid_getCellCoords(JAXIS,blockID, CENTER,.true.,yCenter, blkLimitsGC(HIGH,JAXIS))
        if (NDIM > 2) call Grid_getCellCoords(KAXIS,blockID, CENTER,.true.,zCenter, blkLimitsGC(HIGH,KAXIS))
     elseif (hy_geometry == CARTESIAN) then
        rc = 1.
        rm = 1.
        rp = 1.
     endif

     delyinv = 0.
     delzinv = 0.
     delxinv = 1.0/del(IAXIS)
     if (NDIM > 1) delyinv = 1.0/del(JAXIS)
     if (NDIM > 2) delzinv = 1.0/del(KAXIS)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

              if (hy_geometry == CYLINDRICAL) then
                 rc = abs(xCenter(i))
                 rc = 1./rc
                 rm = abs(xCenter(i) - 0.5*del(IAXIS))
                 rp = abs(xCenter(i) + 0.5*del(IAXIS))
              endif

              if (NDIM > 1) then
                 U(DIVB_VAR,i,j,k) = &
                      (rp*Bx(MAG_FACE_VAR,i+1,j,  k  ) - rm*Bx(MAG_FACE_VAR,i,j,k))*delxinv*rc &
                     +(   By(MAG_FACE_VAR,i,  j+1,k  ) -    By(MAG_FACE_VAR,i,j,k))*delyinv

                 if (NDIM == 3) then
                    U(DIVB_VAR,i,j,k) = U(DIVB_VAR,i,j,k) + &
                      (   Bz(MAG_FACE_VAR,i,  j,  k+1) -    Bz(MAG_FACE_VAR,i,j,k))*delzinv
                 endif
              endif !NDIM > 1

           enddo ! i-loop
        enddo !j-loop
     enddo !k-loop



     deallocate(xCenter)
     deallocate(yCenter)
     if (NDIM == 3) &
     deallocate(zCenter)

     call Grid_releaseBlkPtr(blockID,U,CENTER)
     call Grid_releaseBlkPtr(blockID,Bx,FACEX)
     call Grid_releaseBlkPtr(blockID,By,FACEY)
     if (NDIM == 3) &
     call Grid_releaseBlkPtr(blockID,Bz,FACEZ)

  enddo ! do in = 1, numLeafBlocks
#endif /* NFACE_VARS > 0 */
#endif /* NDIM > 1 */

  return

End Subroutine Hydro_computeVars


