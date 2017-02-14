!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Async/ed_commHandleOffBlkRay
!!
!!  NAME     
!!   ed_commHandleOffBlkRay
!!
!!  SYNOPSIS
!!   ed_commHandleOffBlkRay(integer, intent(IN) :: ray)
!!
!!  DESCRIPTION 
!!    This subroutine moves the ray to the correct block.  If the ray exists 
!!    on a block on my MPI rank then it updates the block ID field in ed_rays.
!!    If the ray exists on a block on another MPI rank then the ray is sent
!!    to that MPI rank and ed_rays is set to NONEXISTENT.
!!    
!!  ARGUMENTS
!!    ray : index of the ray in ed_rays
!!
!!  SIDE EFFECTS
!!    ed_rays is updated.
!!
!!***

#include "constants.h"
#include "EnergyDeposition.h"

subroutine ed_commHandleOffBlkRay(ray)
  use EnergyDeposition_data, ONLY : ed_rays, ed_meshMe
  use ed_commData, ONLY : ed_commRayPosIndex
  !$ use ed_commData, ONLY : ed_commLock
  use UTPipeline, ONLY : UTPipeline_progressComm, UTPipeline_sendItem
  use Grid_interface, ONLY: Grid_getListOfBlocks, Grid_getBlkBoundBox, &
       Grid_outsideBoundBox, Grid_getBlkNeighBlkIDFromPos
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer, intent(IN) :: ray
  real, dimension(LOW:HIGH,1:MDIM) :: bndBox
  real, dimension(MDIM) :: pos
  integer, dimension(MDIM) :: negh
  integer, dimension(BLKNO:PROCNO) :: neghID
  integer :: oldBlockID, blockID, procID
  logical :: outside, isHandled

  pos(:) = ed_rays(ed_commRayPosIndex(:),ray)
  oldBlockID = int(ed_rays(RAY_BLCK,ray))

  call Grid_getBlkBoundBox(oldBlockID, bndBox)
  call Grid_outsideBoundBox(pos, bndBox, outside, negh)
  if (.not.outside) call Driver_abortFlash("Ray still inside the block???")
  call Grid_getBlkNeighBlkIDFromPos(oldBlockID, pos, negh, blockID, procID)

  ed_rays(RAY_BLCK,ray) = blockID
  if (ed_meshMe /= procID) then
     ed_rays(RAY_PROC,ray) = procID
     !$ call omp_set_lock(ed_commLock)
     call UTPipeline_sendItem(ed_rays(:,ray), procID, isHandled)
     do while (.not.isHandled)
        call UTPipeline_progressComm()
        call UTPipeline_sendItem(ed_rays(:,ray), procID, isHandled)
     end do
     !$ call omp_unset_lock(ed_commLock)
     ed_rays(:,ray) = NONEXISTENT
  end if
end subroutine ed_commHandleOffBlkRay
