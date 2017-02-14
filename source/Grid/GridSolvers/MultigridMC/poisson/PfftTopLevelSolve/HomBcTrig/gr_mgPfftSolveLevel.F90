!!****if* source/Grid/GridSolvers/MultigridMC/poisson/PfftTopLevelSolve/HomBcTrig/gr_mgPfftSolveLevel
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
  use gr_pfftInterface, ONLY : gr_pfftDerivs
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_mgPfftData, ONLY : gr_mgPfftInArray, gr_mgPfftOutArray, &
       gr_mgPfftTranArray, gr_mgPfftPoisfact
  use Grid_data,ONLY : gr_domainBC
  use gr_pfftData, ONLY : pfft_transformType, pfft_inLen, pfft_outLen, &
       pfft_globalLen, pfft_usableProc
  use physicaldata, ONLY : unk

  implicit none
  integer, intent(in)    :: iSrc, iSoln, level

  real, allocatable, dimension(:) :: tranArray
  integer, dimension(MDIM) :: localSize,globalSize,transformType
  real :: poisfact
  integer :: inSize,tranSize

  if (pfft_usableProc) then
 
     !First store module level information in local variables.
     transformType(:) = pfft_transformType(:)
     localSize(:) = pfft_inLen(:)
     globalSize(:) = pfft_globalLen(:)
     inSize=pfft_inLen(IAXIS)*pfft_inLen(JAXIS)*pfft_inLen(KAXIS)
     tranSize=2*pfft_outLen(IAXIS)*pfft_outLen(JAXIS)*pfft_outLen(KAXIS)
     poisfact = gr_mgPfftPoisfact

     allocate(tranArray(tranSize))

     !! Here's the real work of the fft
     call Grid_pfftMapToInput(iSrc,gr_mgPfftInArray) 
  
     ! Forward transform of density 
     call Grid_pfft(PFFT_FORWARD,gr_mgPfftInArray,tranArray)

     ! Calculates the transform of iSoln = GPOT
     !  which is the transform of the delSquared(u) = rho in Poisson equation
     call gr_pfftDerivs(tranArray)

     ! Inverse transform of GPOT
     call Grid_pfft(PFFT_INVERSE,tranArray,gr_mgPfftOutArray)
     deallocate(tranArray)

     ! Now multiply by the poisson factor
     gr_mgPfftOutArray(1:inSize) = gr_mgPfftOutArray(1:inSize)*poisfact
  
     ! Map back to the non-uniform mesh
     call Grid_pfftMapFromOutput(iSoln,gr_mgPfftOutArray)

  end if

end subroutine gr_mgPfftSolveLevel


