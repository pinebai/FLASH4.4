!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/Grid_pfftMapToInput
!!
!! NAME 
!!
!!   Grid_pfftMapToInput
!!
!! SYNOPSIS
!!
!!   Grid_pfftMapToInput(integer(IN) :: gridVar,
!!                 real, dimension(:), (OUT) :: pfft_inArray)
!!
!! DESCRIPTION 
!!
!!  Takes the data correnspoding to the variable gridVar in the mesh data 
!!  structure and redistributes to make it campatible with Pfft requirements
!!  based upon the map determined by the routine Grid_pfftInit.
!!
!! ARGUMENTS
!!
!!  gridVar          - variable on the mesh on which pfft is to be applies
!!  pfft_inArray     - array that is input to pfft
!!
!! NOTE 
!!
!!  Users must make sure that Grid_pfftInput has been called correctly
!!  before calling this routine.
!!
!!***

!!REORDER(4): solnData


subroutine Grid_pfftMapToInput(gridVar,pfft_inArray)

  implicit none
  integer,intent(IN) :: gridVar
  real, dimension(:),intent(OUT):: pfft_inArray

  call gr_pfftGridMapToInputInternal(gridVar,pfft_inArray,size(pfft_inArray))

end subroutine Grid_pfftMapToInput


subroutine Grid_pfftMapToInput3DArr(gridVar,pfft_inArray)

  implicit none
  integer,intent(IN) :: gridVar
  real, dimension(:,:,:),intent(OUT):: pfft_inArray

  call gr_pfftGridMapToInputInternal(gridVar,pfft_inArray,size(pfft_inArray))

end subroutine Grid_pfftMapToInput3DArr

#include "constants.h"
#include "Flash.h"
#ifdef Grid_releaseBlkPtr
! doing drift macro expansion by hand because this file has multiple 'use Grid_interface'
! lines. see: drift
#undef Grid_releaseBlkPtr
#endif

subroutine Grid_pfftMapToInputMultivars(gridVar,pfft_inArray,nVars,varList)

  use gr_pfftData, ONLY : pfft_ndim, pfft_needMap, pfft_usableProc
  use Grid_data, ONLY : gr_oneRefLev
#ifndef FLASH_GRID_UG     
  use tree, ONLY : lrefine_max
#endif
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getListOfBlocks, Grid_getBlkCornerID
  use gr_pfftInterface, ONLY : gr_pfftMapToInput
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer,intent(IN) :: gridVar
  real, dimension(:),intent(OUT):: pfft_inArray
  integer,intent(IN) :: nVars
  integer,dimension(*),intent(IN) :: varList

  real, dimension(:,:,:,:), pointer :: solnData
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: blkList(MAXBLOCKS), numBlks
#ifndef FLASH_GRID_UG     
  integer :: cornerID(MDIM), stride(MDIM)
#endif
  integer :: blockID, i, j, k, n, iblk, is, ivar, fanout

  if (.not.pfft_usableProc) return

  if (.not.pfft_needMap) then
     !This was only applicable to UG, now also works for 1D PM - KW
     call Grid_getListOfBlocks(ALL_BLKS,blkList,numBlks)
     do iblk = 1,numBlks
        blockID = blkList(iblk)
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
#ifdef FLASH_GRID_UG     
        n = 0
        fanout = 1
#else
        call Grid_getBlkCornerID(blockID,cornerID,stride)
        fanout = stride(IAXIS)/ 2 ** (lrefine_max - gr_oneRefLev)
        n = (cornerID(IAXIS) - 1) / 2 ** (lrefine_max - gr_oneRefLev)
#endif
        n = n*nVars
        do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                 do is = 1,fanout
                    do ivar=1,nvars
                       n = n + 1
                       pfft_inArray(n) = solnData(varList(ivar),i,j,k)
                    end do
                 end do
              end do
           end do
        end do

        call Driver_driftSetSrcLoc(__FILE__,__LINE__); call Grid_releaseBlkPtr(blockID,solnData,CENTER)

     end do
  else

     if (pfft_ndim == 1) then
        call Driver_abortFlash("[Grid_pfftMapToInput]: Shouldn't be here for 1D!")
     else
        call gr_pfftMapToInput(gridVar, pfft_inArray)
     end if

  end if

end subroutine Grid_pfftMapToInputMultivars



subroutine gr_pfftGridMapToInputInternal(gridVar,pfft_inArray,arraySize)

  use gr_pfftData, ONLY : pfft_ndim, pfft_needMap, pfft_usableProc
  use Grid_data, ONLY : gr_oneRefLev
#ifndef FLASH_GRID_UG     
  use tree, ONLY : lrefine_max
#endif
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getListOfBlocks, Grid_getBlkCornerID
  use gr_pfftInterface, ONLY : gr_pfftMapToInput
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer,intent(IN) :: gridVar, arraySize
  real, dimension(arraySize),intent(OUT):: pfft_inArray

  real, dimension(:,:,:,:), pointer :: solnData
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: blkList(MAXBLOCKS), numBlks
#ifndef FLASH_GRID_UG     
  integer :: cornerID(MDIM), stride(MDIM)
#endif
  integer :: blockID, i, j, k, n, iblk, is, fanout

  if (.not.pfft_usableProc) return

  if (.not.pfft_needMap) then
     !This was only applicable to UG, now also works for 1D PM - KW
     call Grid_getListOfBlocks(ALL_BLKS,blkList,numBlks)
     do iblk = 1,numBlks
        blockID = blkList(iblk)
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
#ifdef FLASH_GRID_UG     
        n = 0
        fanout = 1
#else
        call Grid_getBlkCornerID(blockID,cornerID,stride)
        fanout = stride(IAXIS)/ 2 ** (lrefine_max - gr_oneRefLev)
        n = (cornerID(IAXIS) - 1) / 2 ** (lrefine_max - gr_oneRefLev)
#endif

        do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                 do is = 1,fanout
                    n = n + 1;
!                    if (is.NE.1) &
!                         print*,'n,i,j,k,is:',n,i,j,k,is
                    pfft_inArray(n) = solnData(gridVar,i,j,k)
                 end do
              end do
           end do
        end do

        call Driver_driftSetSrcLoc(__FILE__,__LINE__); call Grid_releaseBlkPtr(blockID,solnData,CENTER)

     end do
  else

     if (pfft_ndim == 1) then
        call Driver_abortFlash("[Grid_pfftMapToInput]: Shouldn't be here for 1D!")
     else
        call gr_pfftMapToInput(gridVar, pfft_inArray)
     end if

  end if

end subroutine Gr_pfftGridMapToInputInternal
