!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhGetBlockWorkload
!!
!! NAME
!!
!!  gr_bhGetBlockWorkload
!!
!!
!! SYNOPSIS
!!
!!   call gr_bhGetTBlockWorkload
!!        integer(in) :: block
!!        )
!!
!! DESCRIPTION
!!
!!   Returns average workload done on a given block by the tree solver.
!!   Used by amr_reorder_block for load balancing.
!!
!! ARGUMENTS
!!
!!  block  - block number
!!
!! RESULT
!!
!!   Average workload done on a given block by the tree solver.
!!
!!***


real function gr_bhGetBlockWorkload(block)
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits
  use gr_bhData, ONLY : gr_bhWrklMinGlob, gr_bhWrklMaxGlob, gr_bhLoadBalancing, &
    gr_bhMaxBlkWeight
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
#include "constants.h"
#include "Flash.h"
  integer,intent(in) :: block
  real, POINTER, DIMENSION(:,:,:,:) :: solnData
  integer, dimension(2,MDIM)   :: blkLimits,blkLimitsGC
  real :: wl

  if (gr_bhLoadBalancing) then
#ifdef GRAV_TREE_LB
    call Grid_getBlkIndexLimits(block,blkLimits,blkLimitsGC) 
    call Grid_getBlkPtr(block, solnData, CENTER) 
    wl = (solnData(WRKL_VAR,blkLimitsGC(LOW,IAXIS),blkLimitsGC(LOW,JAXIS) &
    &  , blkLimitsGC(LOW,KAXIS)) - gr_bhWrklMinGlob) &
    &  / (gr_bhWrklMaxGlob - gr_bhWrklMinGlob)
    wl = 2.0 + (gr_bhMaxBlkWeight - 2.0)*wl
    !print *, "WL: ", block, solnData(WRKL_VAR,blkLimitsGC(LOW,IAXIS) &
    !& ,blkLimitsGC(LOW,JAXIS),blkLimitsGC(LOW,KAXIS)), wl, &
    !& gr_bhWrklMaxGlob, gr_bhWrklMinGlob, solnData(WRKL_VAR &
    !&,blkLimitsGC(LOW,IAXIS),blkLimitsGC(LOW,JAXIS),blkLimitsGC(LOW,KAXIS))
    call Grid_releaseBlkPtr(block, solnData, CENTER)
    gr_bhGetBlockWorkload = wl
#else
    print *, "gr_bhLoadBalancing is True, however, WRKL_VAR does not exist."
    print *, "Run setu script with bhtreeLB=1 and recompile"
    call Driver_abortFlash ('[gr_bhGetBlockWorkload]: Load balancing not compiled in.')

#endif
  else
    gr_bhGetBlockWorkload = 2.0
  endif

  return
end function gr_bhGetBlockWorkload


