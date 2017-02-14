!!****if* source/physics/TreeRay/TreeRayMain/OpticalDepth/tr_odInitFieldVar
!!
!! NAME
!!
!!  tr_odInitFieldVar
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!! RESULT
!!
!!
!!***

subroutine tr_odInitFieldVar()
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr &
  & , Grid_getListOfBlocks
  implicit none
#include "constants.h"
#include "Flash.h"
  integer :: blockID, blockCount, blockList(MAXBLOCKS)
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  
  call Grid_getListOfBlocks(LEAF,blockList,blockCount) ! get list of LEAF blocks
  do blockID = 1, blockCount
    call Grid_getBlkPtr(blockList(blockID),solnData,CENTER)
    solnData(CDTO_VAR, :, :, :) = 0.0
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
    solnData(CHID_VAR, :, :, :) = 0.0
    solnData(CDH2_VAR, :, :, :) = 0.0
#endif
#if CHEMISTRYNETWORK == 5
    solnData(CDCO_VAR, :, :, :) = 0.0
#endif
    call Grid_releaseBlkPtr(blockList(blockID),solnData)
  enddo
  return
end subroutine tr_odInitFieldVar



