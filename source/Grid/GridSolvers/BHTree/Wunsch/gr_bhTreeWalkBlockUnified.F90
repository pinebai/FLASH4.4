!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhTreeWalkBlockUnified
!!
!! NAME
!!
!!  gr_bhTreeWalkBlock
!!
!!
!! SYNOPSIS
!!
!!   gr_bhTreeWalkBlock(
!!          integer,intent(in) :: block
!!          )
!!
!! DESCRIPTION
!!
!!   Traverses the tree in a specific block. In this, "unified tree walk",
!!   the MAC for all interactions is evaluated for each cell, as in the original
!!   Barnes&Hut algorithm.
!!
!! ARGUMENTS
!!
!!  block    - ID of a block where the contribution is calculated
!!
!!***



subroutine gr_bhTreeWalkBlockUnified(block)

  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
  use gr_bhLocalInterface, ONLY : gr_bhTreeWalkPoint, &
    gr_bhStartBlock, gr_bhFinalizeBlock, &
    gr_bhCheckAccuracy
  use Timers_interface, ONLY : Timers_start, Timers_stop

  use gr_bhData, ONLY : gr_bhLocCoords, gr_bhTreeDcount, gr_bhTreeBS, &
    gr_bhTreeNewTW, gr_bhTreeOldAccept, gr_bhWrklMin, gr_bhWrklMax, &
    gr_bhLoadBalancing, gr_bhAcceptAccurateOld, gr_bhWrklDecayI, &
    gr_bhNTestCells, gr_bhTestCells
  !use Gravity_interface, ONLY : Gravity_bhOld2NewBlock

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(IN) :: block
  integer       :: i, ii, jj, kk, ii_mbLp1, jj_mbLp1, kk_mbLp1
  real          :: dcount(1:3)
  real, POINTER, DIMENSION(:,:,:,:) :: solnData
  real, dimension(UNK_VARS_BEGIN:UNK_VARS_END) :: tempCell
  integer, dimension(2,MDIM)   :: blkLimits,blkLimitsGC
  integer, dimension(MDIM) :: point
  integer, dimension(MDIM,gr_bhNTestCells) :: testCellsGC
  logical :: accept
  

  ! get information about the block
  ! indeces necessary to write values into solnData
  call Grid_getBlkIndexLimits(block,blkLimits,blkLimitsGC) 
  ! pointer to the density and gp field
  call Grid_getBlkPtr(block,solnData,CENTER)

  ! call physics routines to do block initialization
  call gr_bhStartBlock(block, blkLimits, solnData)

  dcount(:) = 0.0
  
  ! check test cells
  call Timers_start("check accuracy")
  if (gr_bhAcceptAccurateOld) then
    testCellsGC(IAXIS,:) = gr_bhTestCells(IAXIS,:) + blkLimits(LOW,IAXIS) - 1
    testCellsGC(JAXIS,:) = gr_bhTestCells(JAXIS,:) + blkLimits(LOW,JAXIS) - 1
    testCellsGC(KAXIS,:) = gr_bhTestCells(KAXIS,:) + blkLimits(LOW,KAXIS) - 1

    do i = 1, gr_bhNTestCells
      ! save the grid cell before TreeWalk rewrites it
      tempCell(UNK_VARS_BEGIN:UNK_VARS_END) = &
      & solnData(:, testCellsGC(IAXIS,i), testCellsGC(JAXIS,i), testCellsGC(KAXIS,i))
      call gr_bhTreeWalkPoint(gr_bhLocCoords(gr_bhTestCells(IAXIS,i), IAXIS, block), &
      &                       gr_bhLocCoords(gr_bhTestCells(JAXIS,i), JAXIS, block), &
      &                       gr_bhLocCoords(gr_bhTestCells(KAXIS,i), KAXIS, block), &
      & block, testCellsGC(:,i), blkLimits, solnData, dcount)
      accept = gr_bhCheckAccuracy(block, testCellsGC(:,i), blkLimits, solnData)
      solnData(:, testCellsGC(IAXIS,i), testCellsGC(JAXIS,i), testCellsGC(KAXIS,i)) = &
      & tempCell(UNK_VARS_BEGIN:UNK_VARS_END)
      ! restore the grid cell
      if (.not. accept) exit
    enddo

  else
    accept = .false.
  endif
  call Timers_stop("check accuracy")

  if (accept) then
    ! experimental feature; switched off in this version
    !call Timers_start("old accepted")
    !gr_bhTreeOldAccept = gr_bhTreeOldAccept + 1
    !call Gravity_bhOld2NewBlock(solnData)
    !call Timers_stop("old accepted")
  else
    ! prepare the block for new calculation (i.e. erase rubish created at test cells)
    call Timers_start("new calculated")
    gr_bhTreeNewTW = gr_bhTreeNewTW + 1
    call Timers_start("start block")
    call gr_bhStartBlock(block, blkLimits, solnData)
    call Timers_stop("start block")

    ! recalculate
  
    do kk = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS) ! this block
      point(KAXIS) = kk
      kk_mbLp1 = kk - blkLimits(LOW, KAXIS) + 1
      do jj = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        point(JAXIS) = jj
        jj_mbLp1 = jj - blkLimits(LOW, JAXIS) + 1
        do ii = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
          point(IAXIS) = ii
          ii_mbLp1 = ii - blkLimits(LOW, IAXIS) + 1
 
          call gr_bhTreeWalkPoint( gr_bhLocCoords(ii_mbLp1, IAXIS, block), &
          &                        gr_bhLocCoords(jj_mbLp1, JAXIS, block), &
          &                        gr_bhLocCoords(kk_mbLp1, KAXIS, block), &
          & block, point, blkLimits, solnData, dcount)
     
          !print *, "SS = ", sqrt(ss)
        enddo
      enddo
    enddo
    !print *, "TWBU: new sol calculated, ", block
    call Timers_stop("new calculated")
  endif

  ! call physics routines to do block finalization
  call gr_bhFinalizeBlock(block, blkLimits, solnData)

  ! Load balancing
#ifdef WRKL_VAR
  ! evolve workload using Ornstein-Uhlenbeck scheme
  if (gr_bhLoadBalancing) then
    solnData(WRKL_VAR,:,:,:) = solnData(WRKL_VAR,:,:,:)*exp(-gr_bhWrklDecayI) &
    &  + (dcount(1) + dcount(2) + dcount(3))*(1.-exp(-gr_bhWrklDecayI))
    if (solnData(WRKL_VAR,1,1,1) > gr_bhWrklMax) gr_bhWrklMax = solnData(WRKL_VAR,1,1,1)
    if (solnData(WRKL_VAR,1,1,1) < gr_bhWrklMin) gr_bhWrklMin = solnData(WRKL_VAR,1,1,1)
  endif
#endif

  ! release the block pointer
  call Grid_releaseBlkPtr(block,solnData, CENTER)
  gr_bhTreeDcount = gr_bhTreeDcount + dcount


  return
end subroutine gr_bhTreeWalkBlockUnified

