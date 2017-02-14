!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhBuildTree
!!
!! NAME
!!
!!  gr_bhBuildTree
!!
!!
!! SYNOPSIS
!!
!!   call gr_bhBuildTree()
!!
!! DESCRIPTION
!!
!!   Builds the tree for the Poisson solver. The first part builds trees 
!!   in all LEAF blocks. The second part distributes masses and mass 
!!   centres of all blocks to all CPUs (array gr_bhTreeParentTree). The 
!!   third part calculates masses and mass centres of all parent blocks
!!   (also in the array gr_bhTreeParentTree).
!!
!! ARGUMENTS
!!
!!
!!***



subroutine gr_bhBuildTree()

  use Grid_interface, ONLY : Grid_getListOfBlocks
  use gr_bhData, ONLY : gr_bhTreeParentTree, gr_bhTreeMaxlnblocks, &
    gr_bhTreeNumProcs, gr_bhTreeArray, gr_bhTreeBS, gr_bhTreeLrefine, &
    gr_bhTreeChild, gr_bhTreeLnblocks, gr_bhTreeMyPE, gr_bhTreeNodetype, &
    GR_TREE_NSIZE, gr_bhTreeLevels, gr_bhTreeLoff, &
    GR_TREE_BTP_POS, GR_TREE_BTP_LEV, gr_bhTreeNodeSize, gr_bhBlockTreePos, &
    gr_bhTreeLrefineMax
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use gr_bhLocalInterface, ONLY : gr_bhBuildTreeBlock, gr_bhComParentTree &
  & , gr_bhAccNode, gr_bhNormalizeNode, gr_bhPostprocNode
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  use tree, ONLY : lrefine, child, mchild, nodetype, lnblocks, maxblocks_tr
  !use Logfile_interface, ONLY : Logfile_stamp
      
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"


  integer :: blockID, blockCount, cpu, tr, lev, i, j, ierr
  integer :: tc, cc, pos
  integer :: blockList(gr_bhTreeMaxlnblocks)
  real    :: m, xm, ym, zm
  integer :: ipr, nstep
  logical       :: gcell = .false.
  real,dimension(1:gr_bhTreeBS) :: xCoords, yCoords, zCoords
  real, allocatable :: accnode(:), subnode(:), node(:)

  !call Logfile_stamp(gr_bhTreeMyPE, "Building tree", "[TREE] ")
  call Timers_start("build_tree")

  ! initialization of tree data structures
  call gr_bhBTInit()
  
  call Timers_start("bt1: det_block_trees")

  call Grid_getListOfBlocks(LEAF,blockList,blockCount) ! get list of LEAF blocks

  ! reset gr_bhTreeArray
  do tr = 1, gr_bhTreeMaxlnblocks
    do cpu = 0, gr_bhTreeNumProcs-1
      nullify(gr_bhTreeArray(cpu, tr)%p) ! allocate crashes if gr_bhTreeArray()%p != NULL
    enddo
  enddo

  ! fill gr_bhTreeParentTree array with -1
  gr_bhTreeParentTree(1:GR_TREE_NSIZE,1:gr_bhTreeMaxlnblocks,0:gr_bhTreeNumProcs-1) = -1

  ! call BuildTreeBlock for each LEAF block
  do blockID = 1, blockCount
    call gr_bhBuildTreeBlock(blockList(blockID))
    !print *, "bt1: ", blockID, blockList(blockID) &
    !& , gr_bhTreeLrefine(blockList(blockID),gr_bhTreeMyPE) &
    !& , gr_bhTreeNodeSize(gr_bhTreeLrefine(blockList(blockID), gr_bhTreeMyPE))
  enddo

  call Timers_stop("bt1: det_block_trees")
  call Timers_start("bt2: com_leaf_blocks")
  ! gr_bhTreeParentTree includes masses and MC positions of all LEAF blocks
  ! we need to communicate them to compute masses and MC positions of parent blocks 
  ! (children can be on different CPU)
  call gr_bhComParentTree(-1)
  call Timers_stop("bt2: com_leaf_blocks")

  call Timers_start("bt3: com_parent_blocks")
  allocate(accnode(1:GR_TREE_NSIZE), stat=ierr)
  if (ierr /= 0) call Driver_abortFlash("could not allocate accnode in gr_bhBuildTree.F90")
  allocate(subnode(1:GR_TREE_NSIZE), stat=ierr)
  if (ierr /= 0) call Driver_abortFlash("could not allocate subnode in gr_bhBuildTree.F90")
  do lev = gr_bhTreeLrefineMax-1, 1, -1
    do cpu = 0, gr_bhTreeNumProcs-1
      do tr = 1, gr_bhTreeLnblocks(cpu)
        if ((gr_bhTreeNodetype(tr, cpu) /= 1) .and. (gr_bhTreeLrefine(tr, cpu) .eq. lev)) then
          accnode(1:GR_TREE_NSIZE) = 0.0
        
          do i = 1,mchild
            tc = gr_bhTreeChild(1,i,tr,cpu)
            cc = gr_bhTreeChild(2,i,tr,cpu)
            do j = 1, GR_TREE_NSIZE
              if (gr_bhTreeParentTree(j, tc, cc) == -1) then
                print *, "missing PT, myPE, tr, cpu = ", gr_bhTreeMyPE, tr, cpu, tc, cc
                print *, gr_bhTreeParentTree(:, tc, cc)
                call Driver_abortFlash("missing Parent Tree :(")
              endif
            enddo
            subnode(1:GR_TREE_NSIZE) = gr_bhTreeParentTree(1:GR_TREE_NSIZE, tc, cc)
            call gr_bhAccNode(subnode, accnode)
          enddo
          call gr_bhNormalizeNode(accnode)
          gr_bhTreeParentTree(1:GR_TREE_NSIZE, tr, cpu) = accnode
        endif
      enddo
    enddo
  enddo
  deallocate(subnode)
  deallocate(accnode)
  call Timers_stop("bt3: com_parent_blocks")
  call Timers_start("bt4: postproc_tree")

  call Timers_stop("bt4: postproc_tree")
  allocate(node(1:GR_TREE_NSIZE), stat=ierr)
  if (ierr /= 0) call Driver_abortFlash("could not allocate node in gr_bhBuildTree.F90")
  ! postprocess nodes in leaf blocks
  do blockID = 1, blockCount
    ! go through all nodes with the exception of the bottom level
    ! it is not neccesary to postprocess bottom level nodes now
    ! not decided if there should be a special subroutine gr_bhPostprocBotNode
    ! or if bottom level nodes should be postprocessed by gr_bhPostprocNode
    do i = 1, (8**gr_bhTreeLevels-1)/7 
      !print *, "bt4: ", blockID, blockList(blockID), i
      pos = gr_bhBlockTreePos(GR_TREE_BTP_POS, i)
      !print *, "bt4: pos = ", pos
      lev = gr_bhBlockTreePos(GR_TREE_BTP_LEV, i)
      !print *, "bt4: lev = ", lev
      node = gr_bhTreeArray(gr_bhTreeMyPE, blockList(blockID))%p(pos:pos+GR_TREE_NSIZE-1)
      !print *, "bt4: node = ", node
      call gr_bhPostprocNode(gr_bhTreeNodeSize(lev & 
      & + gr_bhTreeLrefine(blockList(blockID),gr_bhTreeMyPE)), node)
      gr_bhTreeArray(gr_bhTreeMyPE, blockList(blockID))%p(pos:pos+GR_TREE_NSIZE-1) = node
    enddo
  enddo
  ! postprocess nodes in the parent tree
  do cpu = 0, gr_bhTreeNumProcs-1
    do tr = 1, gr_bhTreeLnblocks(cpu)
      !print *, "bt4: cpu, tr = ", cpu, tr
      !print *, "bt4: lev = ", lev
      node = gr_bhTreeParentTree(1:GR_TREE_NSIZE, tr, cpu)
      !print *, "bt4: node1 = ", node
      !print *, "bt4: Lref = ", gr_bhTreeLrefine(tr,cpu)
      call gr_bhPostprocNode(gr_bhTreeNodeSize(gr_bhTreeLrefine(tr,cpu)), node)
      gr_bhTreeParentTree(1:GR_TREE_NSIZE, tr, cpu) = node
      !print *, "bt4: node2 = ", node
    enddo
  enddo
  deallocate(node)
  call Timers_stop("build_tree")

  !call Logfile_stamp(gr_bhTreeMyPE, "Finished building tree", "[TREE] ")  
  return
end subroutine gr_bhBuildTree

