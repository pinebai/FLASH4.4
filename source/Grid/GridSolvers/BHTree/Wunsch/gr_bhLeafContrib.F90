!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhLeafContrib
!!
!! NAME
!!
!!  gr_bhLeafContrib
!!
!!
!! SYNOPSIS
!!
!!   gr_bhLeafContrib(
!!          integer,intent(in)             :: block, tr, cpu, temp
!!          real,intent(INOUT),dimension() :: solnData
!!          integer, dimension(2,MDIM)     :: blkLimits
!!          real,intent(out)               :: cellcnt, nodecnt
!!          )
!!
!! DESCRIPTION
!!
!!   Computes the contribution of a leaf block at a specific point.
!!
!! ARGUMENTS
!!
!!  block     - ID of a block where the contribution is calculated
!!  tr        - ID of a block/tree which contributes
!!  cpu       - cpu of a block/tree which contributes
!!  temp      - template number (2-27; one of 26 3D-neighbours)
!!  solnData  - pointer to field arrays, contribution added to them
!!  blkLimits -
!!  cellcnt   - number of cell-cell interactions
!!  nodecnt   - number of cell-node interactions
!!
!!***



subroutine gr_bhLeafContrib(block, tr, cpu, solnData, blkLimits, cellcnt, nodecnt)

  use gr_bhLocalInterface, ONLY : gr_bhBotNodeContrib, gr_bhNodeContrib, &
    gr_bhMAC, gr_bhSelfContrib, gr_bhPeriodicDr
  use gr_bhData, ONLY : gr_bhTreeLevels, gr_bhTreeArray, gr_bhTreeBS, gr_bhPhysMACTW_step, &
    gr_bhLocCoords, gr_bhTreeMyPE, gr_bhTreeNodeSize, gr_bhTreeNodeSize2, &
    gr_bhTreeBCen, gr_bhTreeLrefine, gr_bhTreeCellSize, gr_bhBlockTreePos, &
    gr_bhLocRecvTreeLevels, &
    GR_TREE_NSIZE, GR_TREE_BNSIZE, GR_TREE_IX, GR_TREE_IY, GR_TREE_IZ, &
    GR_TREE_BTP_POS ,GR_TREE_BTP_LEV, GR_TREE_BTP_I, GR_TREE_BTP_J, GR_TREE_BTP_K, &
    GR_TREE_BTP_C1, GR_TREE_BTP_C8, &
    gr_bhBlockTreeNodeCen
  use Driver_interface, ONLY : Driver_abortFlash
    
  implicit none
#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"


  integer    :: ii, jj, kk
  integer,intent(in) :: block, tr, cpu
  real,intent(out)   :: cellcnt, nodecnt
  real, dimension(:,:,:,:), POINTER_INTENT_IN :: solnData
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  integer       :: i, j, k, l, level, sp, istat, node_btp, pos
  integer :: stack(1:(8**(gr_bhTreeLevels+1)-1)/7)
  integer, dimension(MDIM) :: point
  real ::  dist2
  real, allocatable  :: botnode(:), node(:)
  real, dimension(MDIM) :: nodecoords
  real, dimension(MDIM+2) :: dr ! one additional element for 1/|dr|
  real, dimension(MDIM) :: drDiff, drGC
  logical :: node_accepted

  allocate(botnode(GR_TREE_BNSIZE), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate botnode in gr_bhLeafContrib.F90")
  allocate(node(GR_TREE_NSIZE), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate node in gr_bhLeafContrib.F90")
  cellcnt = 0.0 ! "distance to cell" counter
  nodecnt = 0.0 ! "distance to tree node" counter
  do kk = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS) ! this block
    point(KAXIS) = kk
    do jj = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
      point(JAXIS) = jj
      do ii = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
        point(IAXIS) = ii

        ! put root node on stack
        stack(1) = 1
        sp = 1
   
        ! walk through the tree
        do
          ! take the node on the bottom of the stack
          node_btp = stack(sp)
          sp = sp - 1

          ! get the position in gr_bhTreeArray and the node level
          pos = gr_bhBlockTreePos(GR_TREE_BTP_POS, node_btp)
          level = gr_bhBlockTreePos(GR_TREE_BTP_LEV, node_btp)

          if (level == gr_bhTreeLevels) then
            ! the lowest level: cell-cell interaction

            ! get indeces of the interacting cell in the block
            i = gr_bhBlockTreePos(GR_TREE_BTP_I, node_btp)
            j = gr_bhBlockTreePos(GR_TREE_BTP_J, node_btp)
            k = gr_bhBlockTreePos(GR_TREE_BTP_K, node_btp)
            
            ! get the node array
            botnode = gr_bhTreeArray(cpu, tr)%p(pos:pos+GR_TREE_BNSIZE-1)

            if ((gr_bhTreeMyPE .eq. cpu) .and. (block .eq. tr) &
            .and. (i .eq. (ii-blkLimits(LOW,IAXIS)+1)) &
            .and. (j .eq. (jj-blkLimits(LOW,JAXIS)+1)) &
            .and. (k .eq. (kk-blkLimits(LOW,KAXIS)+1))) then
              ! add self-contribution of the cell
              call gr_bhSelfContrib(botnode, gr_bhTreeCellSize(gr_bhTreeLrefine(tr,cpu),:) &
              & , block, point, blkLimits, solnData)
            else
                
              ! and then calculate the coords from the global gr_bhTreeBCen array
              nodecoords(IAXIS) = gr_bhTreeBCen(IAXIS, tr, cpu) &
              &                 - 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), IAXIS)*(gr_bhTreeBS-2*i+1)
              nodecoords(JAXIS) = gr_bhTreeBCen(JAXIS, tr, cpu) & 
              &                 - 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), JAXIS)*(gr_bhTreeBS-2*j+1)
              nodecoords(KAXIS) = gr_bhTreeBCen(KAXIS, tr, cpu) &
              &                 - 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), KAXIS)*(gr_bhTreeBS-2*k+1)

              dr(IAXIS) = nodecoords(IAXIS) &
              &         - gr_bhLocCoords(point(IAXIS)-blkLimits(LOW,IAXIS)+1, IAXIS, block)
              dr(JAXIS) = nodecoords(JAXIS) &
              &         - gr_bhLocCoords(point(JAXIS)-blkLimits(LOW,JAXIS)+1, JAXIS, block)
              dr(KAXIS) = nodecoords(KAXIS) &
              &         - gr_bhLocCoords(point(KAXIS)-blkLimits(LOW,KAXIS)+1, KAXIS, block)


              ! correct for periodoc conditions
              call gr_bhPeriodicDr(dr)

              dist2  = dr(IAXIS)*dr(IAXIS) + dr(JAXIS)*dr(JAXIS) + dr(KAXIS)*dr(KAXIS)
              dr(MDIM+1) = dist2
              dr(MDIM+2) = sqrt(1.0 / (dist2 + 1D-99))
   
              call gr_bhBotNodeContrib(botnode, level, gr_bhTreeLrefine(tr,cpu) &
              & , dr, block, point, blkLimits, solnData)
  
              cellcnt = cellcnt + 1 ! update the distance counter
            endif
            
          else
            ! higher levels: cell-node interaction
            
            ! read the mass centre coords from the gr_bhTreeArray
            nodecoords(IAXIS) = gr_bhTreeArray(cpu, tr)%p(pos+GR_TREE_IX-1)
            nodecoords(JAXIS) = gr_bhTreeArray(cpu, tr)%p(pos+GR_TREE_IY-1)
            nodecoords(KAXIS) = gr_bhTreeArray(cpu, tr)%p(pos+GR_TREE_IZ-1)
            
            ! compute square of distance
            dr(IAXIS) = nodecoords(IAXIS) &
            &         - gr_bhLocCoords(point(IAXIS)-blkLimits(LOW,IAXIS)+1, IAXIS, block)
            dr(JAXIS) = nodecoords(JAXIS) &
            &         - gr_bhLocCoords(point(JAXIS)-blkLimits(LOW,JAXIS)+1, JAXIS, block)
            dr(KAXIS) = nodecoords(KAXIS) &
            &         - gr_bhLocCoords(point(KAXIS)-blkLimits(LOW,KAXIS)+1, KAXIS, block)
            call gr_bhPeriodicDr(dr)

            ! construct the node array
            node = gr_bhTreeArray(cpu, tr)%p(pos:pos+GR_TREE_NSIZE-1)

            ! position vector to the node mass centre
            dist2  = dr(IAXIS)*dr(IAXIS) + dr(JAXIS)*dr(JAXIS) + dr(KAXIS)*dr(KAXIS)
            dr(MDIM+1) = dist2
            dr(MDIM+2) = sqrt(1.0 / (dist2 + 1D-99))
   
            ! position vector to the node geometric centre
            ! drDiff is a vector from the node mass centre to its geometrical centre
            drDiff(:) = gr_bhTreeBCen(:, tr, cpu) &
            & + gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), :) &
            & * gr_bhBlockTreeNodeCen(:, node_btp) &
            & - node(GR_TREE_IX:GR_TREE_IZ)
            !drDiff(IAXIS) = gr_bhTreeBCen(IAXIS, tr, cpu) &
            !& + gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), IAXIS) &
            !& * gr_bhBlockTreeNodeCen(IAXIS, node_btp) &
            !& - node(GR_TREE_IX)
            !drDiff(JAXIS) = gr_bhTreeBCen(JAXIS, tr, cpu) &
            !& + gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), JAXIS) &
            !& * gr_bhBlockTreeNodeCen(JAXIS, node_btp) &
            !& - node(GR_TREE_IY)
            !drDiff(KAXIS) = gr_bhTreeBCen(KAXIS, tr, cpu) &
            !& + gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), KAXIS) &
            !& * gr_bhBlockTreeNodeCen(KAXIS, node_btp) &
            !& - node(GR_TREE_IZ)
            drGC = dr(1:MDIM) + drDiff
            !print *, "LC: ", dr, "|", drGC, "|", drDiff, "|", gr_bhTreeBCen(:, tr, cpu) &
            !& , "|", gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), :) &
            !& , "|", gr_bhBlockTreeNodeCen(:, node_btp) &
            !& , "|", node(GR_TREE_IX:GR_TREE_IZ)

            ! test the multipole acceptence criterion
            !print *, "Leaf: MAC will be called"
            node_accepted = gr_bhMAC(gr_bhPhysMACTW_step, node &
            & , gr_bhTreeNodeSize2(level+gr_bhTreeLrefine(tr,cpu)) &
            & , drGC, dr, block, point, blkLimits, solnData)

            if (node_accepted) then

              ! add the contribution of the node
              call gr_bhNodeContrib(node, level, gr_bhTreeLrefine(tr,cpu) &
              & , dr, block, point, blkLimits, solnData)
              nodecnt = nodecnt + 1 ! distance counter
            else
              if (cpu /= gr_bhTreeMyPE) then
                if ((level+1) > gr_bhLocRecvTreeLevels(tr, cpu)) then
                  print *, "Asking for higher level than present on a given cpu: ", &
                  level, gr_bhLocRecvTreeLevels(tr, cpu), tr, cpu, gr_bhTreeMyPE
                  print *, "LFF: ", tr, cpu, gr_bhTreeMyPE, sqrt(dr(MDIM+1)), "|", node
                  print *, "LFF2: ", drGC, "|", drDiff, "|", gr_bhTreeBCen(:, tr, cpu) &
                  & , gr_bhBlockTreeNodeCen(:, node_btp)
                  print *, "LFF3: ", dr, "|", &
                       gr_bhLocCoords(point(IAXIS)-blkLimits(LOW,IAXIS)+1, IAXIS, block), &
                       gr_bhLocCoords(point(JAXIS)-blkLimits(LOW,JAXIS)+1, JAXIS, block), &
                       gr_bhLocCoords(point(KAXIS)-blkLimits(LOW,KAXIS)+1, KAXIS, block)
                  call Driver_abortFlash("TREE_ARRAY ACCESS FAILED!!!")
                endif
              endif
              ! put all children on the stack
              do i = GR_TREE_BTP_C1,GR_TREE_BTP_C8
                sp = sp + 1
                stack(sp) = gr_bhBlockTreePos(i, node_btp)
              enddo
            endif
          endif
              
          ! stack is empty - exiting
          if (sp == 0) exit
        enddo

      enddo
    enddo
  enddo
  deallocate(botnode)
  deallocate(node)
  return
end

