#include "constants.h"
#include "Flash.h"

module flash_interfaces
  implicit none
  interface
     subroutine Grid_getBlkPtr(blockID, dataPtr, gridDataStruct)
       implicit none
       integer, intent(in) :: blockID
       real, dimension(:,:,:,:), pointer :: dataPtr
       integer, optional, intent(IN) :: gridDataStruct
     end subroutine Grid_getBlkPtr
  end interface

  interface
     subroutine Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, &
          gridDataStruct)
       implicit none
       integer,intent(IN) :: blockId
       integer, dimension(2,MDIM), intent(OUT) :: blkLimits, blkLimitsGC
       integer, optional, intent(IN) :: gridDataStruct
     end subroutine Grid_getBlkIndexLimits
  end interface
end module flash_interfaces
