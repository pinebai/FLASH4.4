!!****f* source/Grid/Grid_getBlkIDFromPosForListsOfBlocks
!!
!! NAME
!!  Grid_getBlkIDFromPosForListsOfBlocks
!!
!! SYNOPSIS
!!
!!  call Grid_getBlkIDFromPosForListsOfBlocks(real(IN) :: pos(:), 
!!                            integer(IN) :: blkList(:),
!!                            integer(IN) :: blkCount,
!!                            integer(OUT) :: ansBlockID,
!!                            integer(OUT) :: ansProcID,
!!                   optional,integer(IN)  :: comm,
!!                   optional,integer(IN)  :: blockType)
!!  
!!  call Grid_getBlkIDFromPos(real(IN)    :: pos(:), 
!!                            integer(IN) :: blkList(:),
!!                            integer(IN) :: blkCount,
!!                            integer(OUT) :: ansBlockID,
!!                            integer(OUT) :: ansProcID,
!!                   optional,integer(IN)  :: comm,
!!                   optional,integer(IN)  :: blockType)
!!
!! DESCRIPTION 
!! 
!!  Returns the processor and block ID
!!  containing the cell that overlaps with the 
!!  specified position co-ordinate.
!! 
!! 
!! ARGUMENTS 
!!
!!  pos        :: co-ordinates of the point
!!  blkList    :: the list of blocks to search
!!  blkCount   :: the count of blocks in the list
!!  ansBlockID    :: the local block ID of the block that contains the point;
!!                   or NONEXISTENT if no matching block was found by any task.
!!  ansProcID     :: the processor that owns the block that contains the point;
!!                   or NONEXISTENT if no matching block was found by any task.
!!  comm       :: if a communicator other than the default mesh communicator is
!!                desired, it should be specified here.
!!  blockType  :: if given with the value LEAF, indicates that all the blocks in
!!                the passed list are leaf bocks; an implementation of this API
!!                can use this assurance to implement some optimizations.
!!
!!
!! NOTES
!!
!!  Calls must be done collectively, i.e., by all MPI tasks in the communicator
!!  involved (whether given explicitly as comm or defaulted). In a collective call,
!!  all tasks are expected to call with the same values for the coordinate tuple
!!  in pos, but with different block lists.
!!  More succinctly, all intent(in) arguments should have the same values on all
!!  tasks, except for blkList and blkCount.
!!
!!  On return, ansBlockID and ansBlockID will contain the same values on all tasks.
!!
!!  BITTREE is not used in implementations of this interface.
!!
!! EXAMPLE
!!
!!***

#include "constants.h"

subroutine Grid_getBlkIDFromPosForListsOfBlocks(pos,blkList, blkCount,ansBlockID, ansProcID,comm,blockType)

  implicit none

  real, dimension(1:MDIM), intent(IN) :: pos
  integer,intent(IN)  :: blkCount
  integer,dimension(blkCount),intent(IN) :: blkList
  integer, intent(OUT) :: ansBlockID, ansProcID
  integer,optional,intent(IN) :: comm
  integer, optional, intent(IN) :: blockType
  
  ansProcID=NONEXISTENT
  ansBlockID=NONEXISTENT
  
end subroutine Grid_getBlkIDFromPosForListsOfBlocks
