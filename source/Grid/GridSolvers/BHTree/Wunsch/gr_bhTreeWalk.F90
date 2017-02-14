!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhTreeWalk
!!
!! NAME
!!
!!  gr_bhTreeWalk
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhTreeWalk(logical(OUT) :: iterate)
!!
!! DESCRIPTION
!!
!!   Executes the tree walk. Calls gr_treeTreeWalkBlock
!!   for all LEAF blocks.
!!
!! ARGUMENTS
!!
!!  iterate - specifies whether another iteration is needed 
!!            to reached desired accuracy (will be used by 
!!            TreeRay)
!!
!!***



subroutine gr_bhTreeWalk(iterate)

  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
      Grid_getListOfBlocks, Grid_updateRefinement, &
      Grid_fillGuardCells
  use Gravity_interface, ONLY : Gravity_bhTreeWalkEnd
  use TreeRay_interface, ONLY : TreeRay_bhTreeWalkEnd
  use gr_bhData, ONLY : gr_bhTreeDcount, gr_bhTreeZones, &
    gr_bhTreeMyPE, gr_bhComm, gr_bhTreeBS, &
    gr_bhTreeLimangle, gr_bhPhysMACTW_step, gr_bhUseUnifiedTW, &
    gr_bhTWType, GR_TREE_TWSTD, GR_TREE_TWUNI, GR_TREE_TWPQ, &
    gr_bhOAMin, gr_bhOAAvg, gr_bhOAMax, gr_bhOACnt, &
    gr_bhTreeNewTW, gr_bhTreeOldAccept, &
    gr_bhWrklMinGlob, gr_bhWrklMaxGlob, gr_bhWrklMin, gr_bhWrklMax
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp
  use gr_bhLocalInterface, ONLY : gr_bhTreeWalkBlock, &
    gr_bhTreeWalkBlockUnified, gr_bhTreeWalkBlockPQ
      
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  logical, intent(OUT) :: iterate

  integer :: blockID, blockCount
  integer :: blockList(MAXBLOCKS)
  
  integer :: tot_zones, ierr
  integer :: tot_newTW = 0, tot_oldAccept = 0
  real    :: tot_dcount(1:3), tot_OAAvg, tot_OAMin, tot_OAMax, tot_OACnt
  character(len=MAX_STRING_LENGTH) :: strBuff

  call Timers_start("treewalk")

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  ! for calculating numebrs of interactions and average limit angle
  gr_bhTreeDcount = 0.0
  gr_bhTreeZones = 0
  gr_bhOAAvg = 0.0
  gr_bhOAMin = 10.0
  gr_bhOAMax = 0.0
  gr_bhOAcnt = 0
  !gr_bhTotMass = 0.0

  gr_bhTreeNewTW = 0
  gr_bhTreeOldAccept = 0
  tot_newTW = 0
  tot_oldAccept = 0
  gr_bhWrklMax = 0.0
  gr_bhWrklMin = 1.0d99

  do blockID = 1, blockCount
    if (gr_bhTWType == GR_TREE_TWSTD) then
      call gr_bhTreeWalkBlock(blockList(blockID))
    else if (gr_bhTWType == GR_TREE_TWUNI) then
      call gr_bhTreeWalkBlockUnified(blockList(blockID))
    else if (gr_bhTWType == GR_TREE_TWPQ) then
      call gr_bhTreeWalkBlockPQ(blockList(blockID))
    else
    endif
    gr_bhTreeZones = gr_bhTreeZones + gr_bhTreeBS*gr_bhTreeBS*gr_bhTreeBS
  enddo
  !print *, "TW: ", gr_bhTreeMyPE, gr_bhTreeZones, gr_bhTreeDcount
  !print *, "TW totmass: ", gr_bhTreeMyPE, gr_bhTotMass

  ! derivatives of GPOT are calculated for ACEI array
  ! it may be good idea to fill guard cells
  !call Timers_start("mpi reduce")
  call Grid_fillGuardCells(CENTER, ALLDIR)

  call MPI_Reduce(gr_bhTreeDcount,tot_dcount,3,FLASH_REAL,FLASH_SUM,MASTER_PE,gr_bhComm, ierr)
  call MPI_Reduce(gr_bhTreeZones,tot_zones,1,FLASH_INTEGER,FLASH_SUM,MASTER_PE,gr_bhComm, ierr)
  call MPI_Reduce(gr_bhOAAvg,tot_OAAvg,1,FLASH_REAL,FLASH_SUM,MASTER_PE,gr_bhComm, ierr)
  call MPI_Reduce(gr_bhOAMin,tot_OAMin,1,FLASH_REAL,MPI_MIN,MASTER_PE,gr_bhComm, ierr)
  call MPI_Reduce(gr_bhOAMax,tot_OAMax,1,FLASH_REAL,MPI_MAX,MASTER_PE,gr_bhComm, ierr)
  call MPI_Reduce(gr_bhOACnt,tot_OACnt,1,FLASH_REAL,FLASH_SUM,MASTER_PE,gr_bhComm, ierr)

  ! these experimental features are switched off in this version
  !call MPI_Reduce(gr_bhTreeNewTW,tot_newTW,1,FLASH_INTEGER,FLASH_SUM,MASTER_PE,gr_bhComm, ierr)
  !call MPI_Reduce(gr_bhTreeOldAccept,tot_oldAccept,1,FLASH_INTEGER,FLASH_SUM,MASTER_PE,gr_bhComm, ierr)
  !call MPI_AllReduce(gr_bhWrklMax,gr_bhWrklMaxGlob,1,FLASH_REAL,MPI_MAX,gr_bhComm, ierr)
  !call MPI_AllReduce(gr_bhWrklMin,gr_bhWrklMinGlob,1,FLASH_REAL,MPI_MIN,gr_bhComm, ierr)

  !call Timers_stop("mpi reduce")

  if (gr_bhTreeMyPE == MASTER_PE) then
     write (strBuff, '("cell-cell distances: ", e10.3, ", per zone: ", f8.1)') &
     & tot_dcount(1), tot_dcount(1)/tot_zones
     call Logfile_stamp( strBuff, "[BHTree]")
     write (strBuff, '("cell-node distances: ", e10.3, ", per zone: ", f8.1)') &
     & tot_dcount(2), tot_dcount(2)/tot_zones
     call Logfile_stamp( strBuff, "[BHTree]")
     write (strBuff, '("cell-block distances: ", e10.3, ", per zone: ", f8.1)') &
     & tot_dcount(3), tot_dcount(3)/tot_zones
     call Logfile_stamp( strBuff, "[BHTree]")
     tot_OAAvg = tot_OAAvg / tot_OACnt
     write (strBuff, '("opening angle (min,avg,max): ", e10.3, e10.3, e10.3)') &
     & tot_OAMin, tot_OAAvg, tot_OAMax
     call Logfile_stamp( strBuff, "[BHTree]")
     !write (strBuff, '("Blks new calculates, old accepted: ", i6, i6)') &
     !& tot_newTW, tot_oldAccept
     !call Logfile_stamp( strBuff, "[BHTree]")
     !write (strBuff, '("count check: (dcount vs. OACnt): ", e20.10, e20.10)') &
     !& tot_dcount(1)+tot_dcount(2)+tot_dcount(3), tot_OACnt
     !call Logfile_stamp( strBuff, "[BHTree]")
  end if


  call Gravity_bhTreeWalkEnd()
  call TreeRay_bhTreeWalkEnd(iterate)

  call Timers_stop("treewalk")

  return
end subroutine gr_bhTreeWalk


