!!****if* source/Grid/GridSolvers/HYPRE/gr_hypreUpdateSoln
!!
!!  NAME 
!!
!!  gr_hypreUpdateSoln
!!
!!  SYNOPSIS
!!
!!  call gr_hypreUpdateSoln (integer,intent(IN) :: iVar,
!!                           integer,intent(IN) :: blockCount,
!!                           integer,intent(IN) :: blockList (blockCount))
!!
!!  DESCRIPTION 
!!      This routine updates solution after solve (diffusion operation AX=B).
!!
!! ARGUMENTS
!!   iVar       : Variable on which the diffusion operation is performed (e.g., TEMP_VAR)
!!   blockCount : The number of blocks in the list.   
!!   blockList  : The list of blocks on which the solution must be updated.  
!!
!!
!! SIDE EFFECTS
!!
!!  Also sets the HYPRE vector for B to 0.
!!  
!! NOTES
!!
!!   Uses HYPRE library.
!!   Solution is floored if gr_hypreUseFloor is true.
!!
!!***

!!REORDER(4): solnVec


subroutine gr_hypreUpdateSoln (iVar, blockCount, blockList)
  
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_getBlkRefineLevel
  use gr_hypreData,     ONLY : gr_hypreVecX, gr_hypreLower, &
       gr_hypreRefineMIN, gr_hypreUpper, gr_hypreVecB, gr_hypreUseFloor, gr_hypreFloor, &
       gr_hypreSolnIsDelta
  use gr_solversData,   ONLY : dbgContextSolvers => gr_solversDbgContext
  use RadTrans_interface,ONLY: RadTrans_dbgContext_t,RadTrans_getDbgContextPtr
  use Timers_interface, ONLY : Timers_start, Timers_stop    
  
  implicit none
  
#include "Flash.h"  
#include "constants.h"
#include "HYPREf.h"  
  
  integer,intent(IN) :: iVar
  integer,intent(IN) :: blockCount
  integer,intent(IN) :: blockList (blockCount)
  
  !! LOCAL VARIABLES
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits  
  real :: values(1)
  integer :: blockID,part,level,var
  integer :: i, j, k, lb, ierr, pos(NDIM)
  integer :: ii
  real, allocatable, dimension(:) :: BoxVal
  integer :: datasize(MDIM)
  integer :: component
  type(RadTrans_dbgContext_t),pointer :: dbgContextRt
  logical :: skipUpdate
    
  
  call Timers_start("gr_hypreUpdateSoln")    
  
  component = dbgContextSolvers%component
  if(component==3) then
     call RadTrans_getDbgContextPtr(dbgContextRt)
     skipUpdate = ( &
          dbgContextRt%flashErrCode .NE. 0 .AND. &
          dbgContextRt%retriable    .NE. 0 )
  else
     skipUpdate = .FALSE.
  end if

  !!-----------------------------------------------------------------------
  !!    Update solution: small negative numbers are floored to 
  !!    gr_hypreFloor if gr_hypreUseFloor is true.
  !!    Required in the context of MGD.
  !!-----------------------------------------------------------------------  
  
  call HYPRE_SStructVectorGather(gr_hypreVecX, ierr)  
  
  var = 0
  !print*,gr_hypreSolnIsDelta;stop !this is false
  do lb = 1, blockCount     
     blockID = blockList(lb)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     call Grid_getBlkPtr(blockID, solnVec)  
     call Grid_getBlkRefineLevel(blockID,level)                   
     
     part = level - gr_hypreRefineMIN
     
     datasize(1:MDIM)=blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM)+1
     
     allocate(BoxVal(product(dataSize(1:NDIM))))
     
     if (.NOT. skipUpdate) then
        BoxVal = 0.0

        !! Use GetBoxValues more efficient then GetValues.
        call HYPRE_SStructVectorGetBoxValues(gr_hypreVecX, part,gr_hypreLower(lb,1:NDIM), &
             gr_hypreUpper(lb,1:NDIM), var, BoxVal(:), ierr)          

        do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
           do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
              do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)                      

                 ii = (k-blkLimits(LOW,KAXIS)+1)                                +  &
                      (j-blkLimits(LOW,JAXIS))*dataSize(KAXIS)                  +  &
                      (i-blkLimits(LOW,IAXIS))*dataSize(KAXIS)*dataSize(JAXIS)

                 if (gr_hypreSolnIsDelta) then
                    solnVec(iVar,i,j,k) = BoxVal(ii) !! DEV: + ...
                 else
                    solnVec(iVar,i,j,k) = BoxVal(ii)
                 end if
                 if (gr_hypreUseFloor) then
                    if (solnVec(iVar,i,j,k) < gr_hypreFloor) solnVec(iVar,i,j,k) = gr_hypreFloor
                 end if

              end do
           end do
        end do
     end if

     ! Zero out gr_hypreVecB for the next round.
     BoxVal = 0.0
     call HYPRE_SStructVectorSetBoxValues(gr_hypreVecB, part,gr_hypreLower(lb,1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), var, BoxVal(:), ierr) 
     
     deallocate (BoxVal)
     
     call Grid_releaseBlkPtr(blockID, solnVec)
  end do
  
  call Timers_stop("gr_hypreUpdateSoln") 
  
  return
  
end subroutine gr_hypreUpdateSoln
