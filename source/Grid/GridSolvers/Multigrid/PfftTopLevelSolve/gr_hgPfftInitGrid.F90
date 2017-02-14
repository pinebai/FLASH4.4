!!****if* source/Grid/GridSolvers/Multigrid/PfftTopLevelSolve/gr_hgPfftInitGrid
!!
!! NAME
!!
!!  gr_hgPfftInitGrid
!!
!! SYNOPSIS
!!
!!  gr_hgPfftInitGrid(integer(IN) :: refinementLevel, &
!!                    integer(IN) :: gridChanged, &
!!                    integer(IN) :: poisfact)
!!
!! DESCRIPTION
!!
!!  This routine initializes the data necessary for the Multigrid PFFT extensions.
!!
!! ARGUMENTS
!!
!! refinementLevel - Specifies a block refinement level which is passed
!!                   to Grid_pfftInit.  Grid_pfftInit will then generate 
!!                   a map according to this refinement level.
!! gridChanged     - Whether the grid has changed: 0 = NO, 1 = YES.
!! poisfact        - The poisson factor.
!!
!!***

subroutine gr_hgPfftInitGrid(refinementLevel, gridChanged, poisfact)

  use Grid_interface, ONLY : Grid_getGlobalIndexLimits, Grid_pfftInit
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_hgInterface, ONLY : gr_hgPfftFinalizeGrid
  use gr_pfftInterface, ONLY : gr_pfftSpecifyTransform
  use gr_hgPfftData, ONLY : gr_hgPfftInArray, gr_hgPfftOutArray, &
       gr_hgPfftTranArray, gr_hgPfftSolveFlag, gr_hgPfftPoisFact, &
       gr_hgPfftLastMappedLevel, gr_hgBcTypes
  use Grid_data, ONLY : gr_meshNumProcs, gr_meshMe, gr_domainBC

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Pfft.h"

  integer, intent(IN) :: refinementLevel, gridChanged
  real, intent(IN) :: poisfact

  !--------------------------------------------------------------------------
  integer, dimension(MDIM) :: localSize, globalSize, transformType
  integer, dimension(0:MDIM):: baseDatType
  integer :: inSize, error
  logical :: needMap

  gr_hgPfftPoisFact = poisfact
  

  !The variable "gridChanged" takes the value 1 when the 
  !grid has changed.  Otherwise, it has a value 0 when the grid is the same.
  !Recall that gr_hgPfftLastMappedLevel is set to NONEXISTENT on first time step.

  if ( (gr_hgPfftLastMappedLevel == refinementLevel) .and. &
       (gridChanged == 0) ) then 

     if (gr_meshMe == 0) then
        print *, "[gr_hgPfftInitGrid]: Able to retain PFFT grid from last time."
     end if

     !We can re-use the old communication pattern and  
     !the previously allocated space.
     return

  else

     !It may be the case that we allocated space for the following 
     !arrays during a previous call to gr_hgSolve.  The 
     !variable gr_hgPfftLastMappedLevel is initialised to NONEXISTENT 
     !in gr_hgPfftInit(), so we don't deallocate on the first call.
     call gr_hgPfftFinalizeGrid() !deallocates gr_hgPfftInArray and gr_hgPfftOutArray also.


     call Grid_getGlobalIndexLimits(globalSize)

     call gr_pfftSpecifyTransform(transformType, baseDatType, gr_hgBcTypes) !gr_hgBcTypes initialized in gr_hgPfftInit

     needMap=.true.

     !We need to refer to the solve flag in order that we 
     !initialise the PFFT grid appropriately.
     if (gr_hgPfftSolveFlag .eq. 1) then   ! Blktri Solver Requires Processor Grid distributed in Y
        call Grid_pfftInit(NDIM,needMap,globalSize,&
             localSize,transformType, baseDatType=baseDatType, refinementLevel=refinementLevel, &
             jProcs=PFFT_ALL_PROCS,kProcs=PFFT_ONE_PROC)        
     else
        call Grid_pfftInit(NDIM,needMap,globalSize,&
             localSize,transformType, baseDatType=baseDatType, refinementLevel=refinementLevel)
     endif


     inSize = localSize(IAXIS) * localSize(JAXIS) * localSize(KAXIS)
     if (inSize > 0) then        
        !inSize will be 0 on any processor that is not part of the PFFT
        !MPI communicator (as happens if there are more processors than blocks
        !on a particular level).
        allocate(gr_hgPfftInArray(inSize+2), gr_hgPfftOutArray(inSize+2), STAT=error)
        if (error /= 0) then
           call Driver_abortFlash("[gr_hgPfftInitGrid]: Memory cannot be allocated!")
        end if
     end if

  end if

  gr_hgPfftLastMappedLevel = refinementLevel

end subroutine gr_hgPfftInitGrid
