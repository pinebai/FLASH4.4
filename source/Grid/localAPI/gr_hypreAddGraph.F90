!!****if* source/Grid/localAPI/gr_hypreAddGraph
!!
!!  NAME
!!
!!    gr_hypreAddGraph
!!
!!  SYNOPSIS
!!
!!    call gr_hypreAddGraph(integer, intent(IN) :: blockHandle
!!                          integer, intent(IN) :: direction
!!                          integer, intent(IN) :: blockID
!!                          integer, intent(IN) :: blkPartNo
!!                          integer, intent(IN) :: datasize(MDIM)
!!                          integer, intent(IN) :: CornerID(MDIM)
!!                          integer, intent(IN) :: blkStride(MDIM))
!!
!!
!!  DESCRIPTION
!!    This routine handles addition of Graphs (fine-coarse boundaries) across a
!!    a particular face.
!!
!!
!! ARGUMENTS
!!   blockHandle    : Local blockID (within a processor, i.e blockHandle = 1 ...blockCount).
!!   direction      : Face on which Graph objects are to be constructed (ILO_FACE, IHI_FACE, JLO_FACE ....).
!!   blockID        : Global blockID i.e blockList(blockHandle).
!!   blkPartNo      : HYPRE part to which this particular block (leaf) belongs to.
!!   datasize       : No of cells in X,Y,Z direction (i.e NXB, NYB, NZB).
!!   CornerID       : block CornerID.
!!   blkStride      : block Stride.
!!
!! SIDE EFFECTS
!!
!!
!! NOTES:
!!
!!   Uses HYPRE library.
!!
!!***

subroutine gr_hypreAddGraph (blockHandle, blockID, blkPartNo, direction, datasize, CornerID, blkStride, &
     firstHypreVar, numVars)
  
  implicit none
  
#include "constants.h"   
  
  integer, intent(IN) :: blockHandle
  integer, intent(IN) :: direction
  integer, intent(IN) :: blockID   
  integer, intent(IN) :: blkPartNo
  integer, intent(IN) :: datasize(MDIM)
  integer, intent(IN) :: CornerID(MDIM)
  integer, intent(IN) :: blkStride(MDIM)
  integer, intent(IN),OPTIONAL :: firstHypreVar, numVars
  
  
  return
  
end subroutine gr_hypreAddGraph
