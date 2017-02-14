!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/Grid_pfftMapFromOutput
!!
!! NAME 
!!
!!   Grid_pfftMapFromOutput
!!
!! SYNOPSIS
!!
!!   Grid_pfftMapFromOutput(integer(IN) :: gridVar,
!!                 real, dimension(:)(IN)   :: pfft_outArray)
!!
!! DESCRIPTION 
!!
!!  The reverse of Grid_pfftMapToInput:
!!  Places the data from the PPFT array pfft_outArray into 
!!  the mesh at variable gridVar.  Once again it uses the 
!!  map determined by the routine Grid_pfftInit.
!!
!! ARGUMENTS
!!
!!  gridVar          - variable on the mesh on which pfft is to be applies
!!  pfft_outArray     - array that is output from pfft
!!
!! NOTE 
!!
!!  Users must make sure that Grid_pfftInput has been called correctly
!!  before calling this routine.
!!
!!***

!!REORDER(4): solnData

subroutine Grid_pfftMapFromOutput(gridVar,pfft_outArray)
  implicit none

  integer,intent(IN) :: gridVar
  real, dimension(:),intent(IN) :: pfft_outArray

  call gr_pfftGridMapFromOutputIntern(gridVar,pfft_outArray,size(pfft_outArray))

end subroutine Grid_pfftMapFromOutput


subroutine Grid_pfftMapFromOutput3DArr(gridVar,pfft_outArray)
  implicit none

  integer,intent(IN) :: gridVar
  real, dimension(:,:,:),intent(IN) :: pfft_outArray

  call gr_pfftGridMapFromOutputIntern(gridVar,pfft_outArray,size(pfft_outArray))

end subroutine Grid_pfftMapFromOutput3DArr


subroutine gr_pfftGridMapFromOutputIntern(gridVar,pfft_outArray,arraySize)
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getListOfBlocks, Grid_getBlkCornerID
#include "constants.h"
#include "Flash.h"
#include "Pfft.h"
  use gr_pfftData, ONLY : pfft_ndim, pfft_needMap, pfft_usableProc
  use Grid_data, ONLY : gr_oneRefLev
#ifndef FLASH_GRID_UG     
  use tree, ONLY : lrefine_max
#endif
  use gr_pfftInterface, ONLY : gr_pfftMapFromOutput
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer,intent(IN) :: gridVar, arraySize
  real, dimension(arraySize),intent(IN) :: pfft_outArray  

  real, dimension(:,:,:,:), pointer :: solnData
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: blkList(MAXBLOCKS), numBlks
#ifndef FLASH_GRID_UG     
  integer :: cornerID(MDIM), stride(MDIM)
#endif
  integer :: blockID, i, j, k, n, iblk, is, fanin

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
        fanin = 1
#else
        call Grid_getBlkCornerID(blockID,cornerID,stride)
        fanin = stride(IAXIS)/ 2 ** (lrefine_max - gr_oneRefLev)
        n = (cornerID(IAXIS) - 1) / 2 ** (lrefine_max - gr_oneRefLev)
#endif

        do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                 solnData(gridVar,i,j,k) = 0.0
                 do is = 1,max(1,fanin)
                    n = n + 1
                    solnData(gridVar,i,j,k) = solnData(gridVar,i,j,k) + pfft_outArray(n)
                 end do
                 if (fanin>1) &
                      solnData(gridVar,i,j,k) = solnData(gridVar,i,j,k) / fanin
              end do
           end do
        end do

        call Grid_releaseBlkPtr(blockID,solnData,CENTER)

     end do
  else

     if (pfft_ndim == 1) then
        call Driver_abortFlash("[Grid_pfftMapFromOutput]: Shouldn't be here for 1D!")
     else
        call gr_pfftMapFromOutput(gridVar, pfft_outArray)
     end if

  end if

end subroutine Gr_pfftGridMapFromOutputIntern
