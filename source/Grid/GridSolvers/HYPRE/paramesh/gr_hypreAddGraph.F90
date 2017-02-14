!!****if* source/Grid/GridSolvers/HYPRE/paramesh/gr_hypreAddGraph
!!
!!  NAME 
!!
!!    gr_hypreAddGraph
!!
!!  SYNOPSIS
!!
!!    call gr_hypreAddGraph(         integer(IN) :: blockHandle,
!!                                   integer(IN) :: direction,
!!                                   integer(IN) :: blockID,
!!                                   integer(IN) :: blkPartNo,
!!                                   integer(IN) :: datasize(MDIM),
!!                                   integer(IN) :: CornerID(MDIM),
!!                                   integer(IN) :: blkStride(MDIM),
!!                          OPTIONAL,integer(IN) :: firstHypreVar,
!!                          OPTIONAL,integer(IN) :: numVars)
!!
!!
!!  DESCRIPTION 
!!    This routine handles addition of Graphs (fine-coarse boundaries) across a
!!    a particular face.
!! 
!!
!! ARGUMENTS
!!   blockHandle    : Local blockID (within a processor, i.e blockHandle = 1 ... blockCount).
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
  
  
  
  use Grid_interface,   ONLY : Grid_getBlkIndexLimits, Grid_getBlkCornerID  
  
  use gr_hypreData,     ONLY : gr_hypreLower, gr_hypreUpper, gr_hypreVecX, &
                               gr_hypreVecB, gr_hypreMatA, gr_hypreGraph, gr_hypreNParts, &
                               gr_hypreGrid, gr_hypreRefineMIN, gr_hypreRefineMAX, gr_hypreNeghLevels, &
                               gr_hypreSolverType, gr_hypreSurrBlkSum
  
  use Grid_data,        ONLY : gr_meshComm, gr_meshMe  

  use gr_interfaceTypeDecl
  use gr_interface,     ONLY : gr_findAllNeghID, gr_getBlkHandle  
  use tree,             ONLY:  surr_blks, lrefine_max 
  
  
  implicit none
  
#include "Flash.h"
#include "constants.h"   
#include "HYPREf.h"   
#include "Flash_mpi.h"  
  
  integer, intent(IN) :: blockHandle
  integer, intent(IN) :: direction
  integer, intent(IN) :: blockID   
  integer, intent(IN) :: blkPartNo
  integer, intent(IN) :: datasize(MDIM)
  integer, intent(IN) :: CornerID(MDIM)
  integer, intent(IN) :: blkStride(MDIM)
  integer, intent(IN),OPTIONAL :: firstHypreVar, numVars


  integer, dimension(MDIM) :: BoxLow, BoxHigh,dir
  integer, dimension(MDIM) :: from_index, to_index
  integer, dimension(MDIM) :: neghStride
  integer, dimension(MDIM) :: iLower, iUpper   
  integer :: i,j,k, ii, jj, GraphOffset, blkNumNegh
  integer :: GraphDirection, SecondDirection, ThirdDirection
  integer :: eachNegh, ierr, neghRefineLevel, neghPartNO
  integer :: fvar, nvars, var

  if (present(firstHypreVar)) then
     fvar = firstHypreVar
  else
     fvar = 0
  end if
  if (present(numVars)) then
     nvars = numVars
  else
     nvars = 1
  end if
  
  BoxLow(IAXIS)  = gr_hypreLower(blockHandle,IAXIS)
  BoxLow(JAXIS)  = gr_hypreLower(blockHandle,JAXIS)
  BoxLow(KAXIS)  = gr_hypreLower(blockHandle,KAXIS)
  
  BoxHigh(IAXIS) = gr_hypreUpper(blockHandle,IAXIS)
  BoxHigh(JAXIS) = gr_hypreUpper(blockHandle,JAXIS)
  BoxHigh(KAXIS) = gr_hypreUpper(blockHandle,KAXIS)
  
  dir(IAXIS) = CENTER
  dir(JAXIS) = 1 + K2D
  dir(KAXIS) = 1 + K3D

  
  select case (direction)     
  case (ILO_FACE)
     BoxHigh(NDIM-IAXIS+1) = BoxLow (NDIM-IAXIS+1)
     dir(IAXIS) = LEFT_EDGE
     GraphDirection  = IAXIS
     SecondDirection = JAXIS
     ThirdDirection  = KAXIS
     from_index (NDIM-IAXIS+1) = BoxLow (NDIM-IAXIS+1)                           
  case (IHI_FACE)
     BoxLow(NDIM-IAXIS+1)  = BoxHigh(NDIM-IAXIS+1) 
     dir(IAXIS) = RIGHT_EDGE
     GraphDirection  = IAXIS
     SecondDirection = JAXIS
     ThirdDirection  = KAXIS
     from_index (NDIM-IAXIS+1) = BoxHigh(NDIM-IAXIS+1)                           
#if NDIM >=2
  case (JLO_FACE)
     BoxHigh(NDIM-JAXIS+1) = BoxLow (NDIM-JAXIS+1) 
     dir(JAXIS) = LEFT_EDGE
     GraphDirection  = JAXIS
     SecondDirection = IAXIS
     ThirdDirection  = KAXIS
     from_index (NDIM-JAXIS+1) = BoxLow (NDIM-JAXIS+1)                           
  case (JHI_FACE)
     BoxLow(NDIM-JAXIS+1)  = BoxHigh(NDIM-JAXIS+1)
     dir(JAXIS) = RIGHT_EDGE
     GraphDirection  = JAXIS
     SecondDirection = IAXIS
     ThirdDirection  = KAXIS
     from_index (NDIM-JAXIS+1) = BoxHigh(NDIM-JAXIS+1)                           
#if NDIM == 3
  case (KLO_FACE)     
     BoxHigh(NDIM-KAXIS+1) = BoxLow (NDIM-KAXIS+1)
     dir(KAXIS) = LEFT_EDGE
     GraphDirection  = KAXIS
     SecondDirection = IAXIS
     ThirdDirection  = JAXIS
     from_index (NDIM-KAXIS+1) = BoxLow(NDIM-KAXIS+1)                           
  case (KHI_FACE)
     BoxLow(NDIM-KAXIS+1)  = BoxHigh(NDIM-KAXIS+1)
     dir(KAXIS) = RIGHT_EDGE
     GraphDirection  = KAXIS     
     SecondDirection = IAXIS
     ThirdDirection  = JAXIS
     from_index (NDIM-KAXIS+1) = BoxHigh(NDIM-KAXIS+1)                           
#endif
#endif
  end select
  
#if NDIM >= 2 
  SecondDirection = NDIM - SecondDirection + 1
#if NDIM == 3
  ThirdDirection  = NDIM - ThirdDirection + 1
#endif
#endif

  blkNumNegh = gr_hypreSurrBlkSum(blockHandle) % regionInfo(dir(IAXIS),dir(JAXIS),dir(KAXIS)) % numNegh  
  if (blkNumNegh == 0) return
  
  neghRefineLevel = gr_hypreNeghLevels(blockHandle,dir(IAXIS),dir(JAXIS),dir(KAXIS))
  neghPartNO = neghRefineLevel - gr_hypreRefineMIN

  if (neghPartNO == blkPartNO) return
  
  neghStride = 2**(lrefine_max - neghRefineLevel)   
  
  if (direction == ILO_FACE .or. direction == JLO_FACE .or. direction == KLO_FACE) then
     GraphOffset = -1     
     ilower(GraphDirection) = CornerID(GraphDirection)-neghStride(GraphDirection)*(datasize(GraphDirection)+1)      
     ilower(GraphDirection) = ceiling(real(ilower(GraphDirection)) / real(neghStride(GraphDirection)))
     iupper(GraphDirection) = ilower(GraphDirection) + datasize (GraphDirection)                                  
     
     to_index(NDIM-GraphDirection+1) = iupper(GraphDirection)   
     
  else
     GraphOffset = +1          
     ilower(GraphDirection) = CornerID(GraphDirection) + blkStride(GraphDirection)*(datasize(GraphDirection)+1)  
     ilower(GraphDirection) = ceiling(real(ilower(GraphDirection)) / real(neghStride(GraphDirection)))  
     iupper(GraphDirection) = ilower(GraphDirection) + datasize (GraphDirection)   
     
     to_index (NDIM-GraphDirection+1) = ilower(GraphDirection)                                                               
  end if
  
  !! 1 in 1D, 2 in 2D and 4 in 3D  
  do j = BoxLow(ThirdDirection), BoxHigh(ThirdDirection)
     do i = BoxLow(SecondDirection), BoxHigh(SecondDirection)          
        
        from_index (SecondDirection) = i
        from_index (ThirdDirection)  = j        
        
        if (neghPartNO > blkPartNO) then !! coarse->fine graph
           do jj = j-1, j*K3D
              do ii = i - 1, i*K2D
                 to_index (SecondDirection) = i+ii
                 to_index (ThirdDirection)  = j+jj                 

                 do var=fvar,fvar+nvars-1
                    call HYPRE_SStructGraphAddEntries (gr_hypreGraph, blkPartNO, from_index(1:NDIM),var, &
                      neghPartNO,to_index(1:NDIM), var, ierr)
                 end do
              end do
           end do
        else !! fine->coarse graph.           
           to_index (SecondDirection) = ceiling(real(i) / real(neghStride(SecondDirection)/blkStride(SecondDirection)))
           to_index (ThirdDirection)  = ceiling(real(j) / real(neghStride(ThirdDirection)/blkStride(ThirdDirection)))           
           
           do var=fvar,fvar+nvars-1
              call HYPRE_SStructGraphAddEntries (gr_hypreGraph, blkPartNO, from_index(1:NDIM),var, &
                   neghPartNO,to_index(1:NDIM), var, ierr)
           end do
        end if
     end do
  end do
  
end subroutine gr_hypreAddGraph
