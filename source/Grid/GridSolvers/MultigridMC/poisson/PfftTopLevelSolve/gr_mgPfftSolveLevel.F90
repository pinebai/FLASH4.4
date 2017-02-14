!!****if* source/Grid/GridSolvers/MultigridMC/poisson/PfftTopLevelSolve/gr_mgPfftSolveLevel
!!
!! NAME
!!
!!  gr_mgPfftSolveLevel
!!
!! SYNOPSIS
!!
!!  gr_mgPfftSolveLevel(integer(IN) :: iSrc, &
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

subroutine gr_mgPfftSolveLevel (iSrc, iSoln, level)

  use Grid_interface, ONLY : Grid_pfftMapToInput,&
       Grid_pfftInit, Grid_pfft, Grid_pfftMapFromOutput
  use gr_pfftInterface, ONLY : gr_pfftDerivs, gr_pfftPoissonDirect
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_mgPfftData, ONLY : gr_mgPfftInArray, gr_mgPfftOutArray, &
       gr_mgPfftTranArray, gr_mgPfftSolveFlag, gr_mgPfftPoisfact

  use Grid_data,ONLY : gr_domainBC

  use gr_pfftData, ONLY : pfft_transformType, pfft_inLen, pfft_globalLen, &
       pfft_usableProc
  use physicaldata, ONLY : unk

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
     solveFlag = gr_mgPfftSolveFlag
     inSize=pfft_inLen(IAXIS)*pfft_inLen(JAXIS)*pfft_inLen(KAXIS)
     poisfact = gr_mgPfftPoisfact


     !! Here's the real work of the fft
     call Grid_pfftMapToInput(iSrc,gr_mgPfftInArray) 
  
     ! Call gr_pfftPoisson Direct Forward:
     call gr_pfftPoissonDirect (PFFT_FORWARD, solveflag, inSize, localSize, globalSize, &
          transformType, gr_mgPfftInArray, gr_mgPfftOutArray)

     ! Call gr_pfftPoisson Direct Inverse: 
     call gr_pfftPoissonDirect (PFFT_INVERSE, solveflag, inSize, localSize, globalSize, & 
          transformType, gr_mgPfftInArray, gr_mgPfftOutArray)

     ! Now multiply by the poisson factor
     gr_mgPfftOutArray(1:inSize) = gr_mgPfftOutArray(1:inSize)*poisfact
  
     ! Map back to the non-uniform mesh
     call Grid_pfftMapFromOutput(iSoln,gr_mgPfftOutArray)

  end if

end subroutine gr_mgPfftSolveLevel
