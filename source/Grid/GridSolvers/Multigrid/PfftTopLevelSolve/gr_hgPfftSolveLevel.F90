!!****if* source/Grid/GridSolvers/Multigrid/PfftTopLevelSolve/gr_hgPfftSolveLevel
!!
!! NAME
!!
!!  gr_hgPfftSolveLevel
!!
!! SYNOPSIS
!!
!!  gr_hgPfftSolveLevel(integer(IN) :: iSrc, &
!!                      integer(IN) :: iSoln, &
!!                      integer(IN) :: level)
!!
!! DESCRIPTION
!!
!! This subroutine must enclose the work of the PFFT (directSolver) version of 
!! Grid_solvePoisson.  We must do this because by using Multigrid 
!! we require Multigrid's implementation of Grid_solvePoisson.
!!
!! ARGUMENTS
!!
!! iSrc  - Specifies the source grid element.
!! iSoln - Specifies the solution grid element.
!! level - The refinement level at which we are solving.
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "Pfft.h"  

subroutine gr_hgPfftSolveLevel (iSrc, iSoln, level)

  use Grid_interface, ONLY : Grid_pfftMapToInput,&
       Grid_pfftInit, Grid_pfft, Grid_pfftMapFromOutput
  use gr_pfftInterface, ONLY : gr_pfftDerivs, gr_pfftPoissonDirect
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_hgPfftData, ONLY : gr_hgPfftInArray, gr_hgPfftOutArray, &
       gr_hgPfftTranArray, gr_hgPfftSolveFlag, gr_hgPfftPoisfact
  use gr_pfftData, ONLY : pfft_transformType, pfft_inLen, pfft_globalLen, &
       pfft_usableProc

  implicit none
  integer, intent(in)    :: iSrc, iSoln, level

  integer, dimension(MDIM) :: localSize,globalSize,transformType
  real :: poisfact
  integer :: inSize, solveflag

  if (pfft_usableProc) then

     !First store module level information in local variables.
     transformType(:) = pfft_transformType(:)
     localSize(:) = pfft_inLen(:)
     globalSize(:) = pfft_globalLen(:)
     solveFlag = gr_hgPfftSolveFlag
     inSize=pfft_inLen(IAXIS)*pfft_inLen(JAXIS)*pfft_inLen(KAXIS)
     poisfact = gr_hgPfftPoisfact


     !! Here's the real work of the fft
     call Grid_pfftMapToInput(iSrc,gr_hgPfftInArray) 
  
     ! Call gr_pfftPoisson Direct Forward:
     call gr_pfftPoissonDirect (PFFT_FORWARD, solveflag, inSize, localSize, globalSize, &
          transformType, gr_hgPfftInArray, gr_hgPfftOutArray)

     ! Call gr_pfftPoisson Direct Inverse: 
     call gr_pfftPoissonDirect (PFFT_INVERSE, solveflag, inSize, localSize, globalSize, & 
          transformType, gr_hgPfftInArray, gr_hgPfftOutArray)

     ! Now multiply by the poisson factor
     gr_hgPfftOutArray(1:inSize) = gr_hgPfftOutArray(1:inSize)*poisfact
  
     ! Map back to the non-uniform mesh
     call Grid_pfftMapFromOutput(iSoln,gr_hgPfftOutArray)

  end if

  call CopyFromUnkToWork(iSoln, iSrc, level)

end subroutine gr_hgPfftSolveLevel




!-------------------------------------------------------------------


subroutine CopyFromUnkToWork(iSoln, iSource, level)

  use Grid_data, ONLY: gr_meshMe
  use Grid_interface, ONLY : Grid_getListOfBlocks
  use gr_hgInterface, ONLY : gr_hgGuardCell
  use workspace, ONLY : work
  use physicaldata, ONLY : unk

  implicit none

  integer, intent(IN) :: iSoln, iSource, level

  integer :: blockList(MAXBLOCKS)
  integer :: i, j, k, blk, blockID, blockCount

  call Grid_getListOfBlocks(REFINEMENT,blockList,blockCount,level)
  
  do blk = 1, blockCount
     blockID = blockList(blk)
     do k = 1, NZB
        do j = 1, NYB
           do i = 1, NXB
              work(i+NGUARD,j+K2D*NGUARD,k+K3D*NGUARD,blockID,1) = &
                   unk(iSoln,i+NGUARD,j+K2D*NGUARD,k+K3D*NGUARD,blockID) 
           enddo
        enddo
     enddo
!!$     work(NGUARD+1:NXB+NGUARD,NGUARD+1:NYB+NGUARD,NGUARD+1:NZB+NGUARD,blockID,1) = &
!!$     unk(iSoln,NGUARD+1:NXB+NGUARD,NGUARD+1:NYB+NGUARD,NGUARD+1:NZB+NGUARD,blockID)
  end do

  call gr_hgGuardCell(gr_meshMe, 1, .false.)

  return

end subroutine CopyFromUnkToWork
