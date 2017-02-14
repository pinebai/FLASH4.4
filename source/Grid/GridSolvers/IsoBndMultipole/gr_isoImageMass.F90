!!****if* source/Grid/GridSolvers/IsoBndMultipole/gr_isoImageMass
!!
!! NAME
!!  gr_isoImageMass
!!
!! SYNOPSIS
!!
!!  gr_isoImageMass(integer(IN) :: isoln,
!!                           integer(IN) :: iiden)
!!
!! DESCRIPTION
!!
!!  Given a zero-boundary potential in the solution variable with
!!  index isoln, compute the surface density of the image mass
!!  required to produce this potential.  The result is written into
!!  the solution variable with index iiden.
!!
!! ARGUMENTS
!!
!!  isoln  - variable index of the solution variable 
!!  iiden  - variable index of the density variable 
!!
!!***

!!REORDER(4): solnData

subroutine gr_isoImageMass (isoln, iiden)

  !=======================================================================

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,&
                             Grid_getListOfBlocks,Grid_getBlkBC,&
                             Grid_getDeltas,Grid_getBlkIndexLimits

  implicit none

#include "Flash.h"
#include "constants.h"

  integer,intent(IN)       :: isoln, iiden

  integer       :: i, j, k, lb, blkCount,blockID
  integer,dimension(LOW:HIGH,MDIM)::faces
  integer, dimension(MAXBLOCKS) :: blkList
  real,dimension(MDIM) :: delta
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real, pointer :: solnData(:,:,:,:)
  real          ::  dx, dy, dz
  integer       :: lowGuard,hiGuard

  !=====================================================================

  call Grid_getListOfBlocks(LEAF,blkList,blkCount)
  do lb = 1, blkCount
     blockID=blkList(lb)
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getDeltas(blockID,delta)
     call Grid_getBlkBC(blockID,faces)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     lowGuard=blkLimits(LOW,IAXIS)-1
     hiGuard=blkLimits(HIGH,IAXIS)+1

     solnData(iiden,:,:,:) = 0

     dx = delta(IAXIS)
     
     if (faces(LOW,IAXIS)/=NOT_BOUNDARY) then       ! -x boundary
        solnData(iiden,lowGuard,:,:) = &
             3.5/dx * solnData(isoln,lowGuard+1,:,:)  &
             -0.5/dx * solnData(isoln,lowGuard+2,:,:)
        solnData(iiden,lowGuard+1,&
                 blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                 blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
                 0.5/dx * solnData(iiden,lowGuard,&
                 blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                 blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
     endif
     
     if (faces(HIGH,IAXIS) /=NOT_BOUNDARY) then       ! +x boundary
        solnData(iiden,hiGuard,:,:) = &
             3.5/dx * solnData(isoln,hiGuard-1,:,:)  &
             -0.5/dx * solnData(isoln,hiGuard-2,:,:)
        solnData(iiden,hiGuard-1,&
                 blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                 blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
                 0.5/dx * solnData(iiden,hiGuard,&
                 blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                 blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
     endif
     
     if (NDIM >= 2) then
        dy = delta(JAXIS)
        lowGuard=blkLimits(LOW,JAXIS)-1
        hiGuard=blkLimits(HIGH,JAXIS)+1
        if (faces(LOW,JAXIS) /=NOT_BOUNDARY) then       ! -y boundary
           solnData(iiden,:,lowGuard,:) = &
                3.5/dy * solnData(isoln,:,lowGuard+1,:)  &
                -0.5/dy * solnData(isoln,:,lowGuard+2,:)
           solnData(iiden,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                lowGuard+1,blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
                solnData(iiden,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                lowGuard+1,blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) + &
                0.5/dy * solnData(iiden,&
                blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                lowGuard,blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
        endif
        
        if (faces(HIGH,JAXIS) /=NOT_BOUNDARY) then       ! +y boundary
           solnData(iiden,:,hiGuard,:) = &
                3.5/dy * solnData(isoln,:,hiGuard-1,:)  &
                -0.5/dy * solnData(isoln,:,hiGuard-2,:)
           solnData(iiden,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                hiGuard-1,blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
                solnData(iiden,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                hiGuard-1,blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) + &
                0.5/dy * solnData(iiden,&
                blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),hiGuard,&
                blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
        endif
     endif
     
     if (NDIM == 3) then
        dz = delta(KAXIS)
        lowGuard=blkLimits(LOW,KAXIS)-1
        hiGuard=blkLimits(HIGH,KAXIS)+1

        if (faces(LOW,KAXIS) /=NOT_BOUNDARY) then       ! -z boundary
           solnData(iiden,:,:,lowGuard) = &
                3.5/dz * solnData(isoln,:,:,lowGuard+1)  &
                -0.5/dz * solnData(isoln,:,:,lowGuard+2)
           solnData(iiden,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),lowGuard+1) = &
                solnData(iiden,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),lowGuard+1) + &
                0.5/dz * solnData(iiden,&
                blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),lowGuard)
        endif
        
        if (faces(HIGH,KAXIS) /=NOT_BOUNDARY) then       ! +z boundary
           solnData(iiden,:,:,hiGuard) = &
                3.5/dz * solnData(isoln,:,:,hiGuard-1)  &
                -0.5/dz * solnData(isoln,:,:,hiGuard-2)
           solnData(iiden,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),hiGuard-1) = &
                solnData(iiden,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),hiGuard-1) + &
                0.5/dz * solnData(iiden,&
                blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),hiGuard)
        endif
     endif
     
     call Grid_releaseBlkPtr(blockID, solnData)
     
  enddo
  
  !=======================================================================
  
  return
end subroutine gr_isoImageMass
