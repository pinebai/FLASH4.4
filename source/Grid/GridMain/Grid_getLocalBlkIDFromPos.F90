!!****if* source/Grid/GridMain/Grid_getLocalBlkIDFromPos
!!
!! NAME
!!  Grid_getLocalBlkIDFromPos
!!
!! SYNOPSIS
!!
!!  call Grid_getLocalBlkIDFromPos(real(IN) :: pos(:), 
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
#include "Flash.h"


! A simple version that just searches through a block list until a block matching the position is found.
! Some shortcuts may be taken if the position is obviously outsise the domain or if we are only
! searching for LEAF nodes.
subroutine Grid_getLocalBlkIDFromPosSimple(pos,ansBlockID, ansProcID,blkList, blkCount,blockType)

  use Grid_data, ONLY : gr_minCellSizes, gr_globalDomain, gr_meshMe
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkCornerID
  implicit none

  real, dimension(1:MDIM), intent(IN) :: pos
  integer, intent(OUT) :: ansBlockID, ansProcID
  integer, OPTIONAL,intent(IN)  :: blkCount
  integer, OPTIONAL,dimension(:),intent(IN),target :: blkList
  integer, OPTIONAL, intent(IN) :: blockType

  integer  :: blkCountEff
  integer,POINTER,dimension(:) :: blkListEff
  integer :: blockTypeEff

  integer :: i !,j,blk, ierr, mycomm
  integer, dimension(MDIM) :: cornerID, stride, cornerIDHigh, posInd
  logical,dimension(MDIM) :: inBlk
  logical :: onUpperBoundary(NDIM), found
  real,dimension(MDIM) :: temp

  if (present(blockType)) then
     blockTypeEff = blockType
  else
     blockTypeEff = LEAF
  end if

  ansBlockID=NONEXISTENT
  ansProcID=NONEXISTENT
  found = .FALSE.
  nullify(blkListEff)


  if (present(blkCount)) then
     blkCountEff = blkCount
  else
     blkCountEff = MAXBLOCKS
     if (present(blkList)) blkCountEff = min(blkCountEff,size(blkList))
  end if
  if (blkCountEff==0) goto 3

  if (present(blkList)) then
     blkListEff => blkList
  else
     allocate(blkListEff(blkCountEff))
     call Grid_getListOfBlocks(LEAF,blkListEff,blkCountEff)
  end if
  if (blkCountEff==0) goto 2


  temp(1:NDIM)=(gr_globalDomain(HIGH,1:NDIM) - pos(1:NDIM))
  onUpperBoundary(:) = (temp(1:NDIM)==0.0)
  if(minval(temp).ge.0.0) then
     temp(1:NDIM)=(pos(1:NDIM)-gr_globalDomain(LOW,1:NDIM))/gr_minCellSizes(1:NDIM)
     if(minval(temp).ge.0.0) then
        posInd(1:NDIM)=int(temp(1:NDIM))+1
        where (onUpperBoundary(:)) posInd(1:NDIM)=posInd(1:NDIM)-1

        do i = 1,blkCountEff
           ansBlockID=blkListEff(i)
           call Grid_getBlkCornerID(ansBlockID,cornerID,stride,cornerIDHigh)
           inBlk=.true.
           inBlk(1:NDIM)=(posInd(1:NDIM).ge.cornerID(1:NDIM)).and.&
                         (posInd(1:NDIM).le.cornerIDHigh(1:NDIM))
           found=inBlk(IAXIS).and.inBlk(JAXIS).and.inBlk(KAXIS)
           if(found) then
              EXIT
           end if
        end do
     end if
  end if

  if (found) then
     ansProcID = gr_meshMe
  else
     ansBlockID = NONEXISTENT
  end if

2 if (.NOT. present(blkList)) then
     if (associated(blkListEff)) deallocate(blkListEff)
  end if

3 return

end subroutine Grid_getLocalBlkIDFromPosSimple
