!!****if* source/Grid/GridMain/Grid_getBlkIDFromPos
!!
!! NAME
!!  Grid_getBlkIDFromPos
!!
!! SYNOPSIS
!!
!!  call Grid_getBlkIDFromPos(real(IN)     :: pos(:), 
!!                            integer(OUT) :: ansBlockID,
!!                            integer(OUT) :: ansProcID,
!!                   optional,integer(IN)  :: comm)
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
!!  ansBlockID    :: the local blockID of the block that contains the point;
!!                   or NONEXISTENT if no matching block was found globally.
!!  ansProcID     :: the ID of the processor that contains the point;
!!                   or NONEXISTENT if no matching block was found anywhere.
!!  comm       :: if communication is necessary, a communicator must be
!!                specified here.
!!                Communication is necessary unless the setup uses BITTREE.
!!                The comm argument is ignored if BITTREE is used.
!!
!! NOTES
!!
!!  If a communicator comm is present and used, then this routine must be
!!  called collectively by all MPI tasks in the communicator. In a collective
!!  call, all tasks are expected to call with the same values for the
!!  coordinate tuple in pos. In other words, The intent(in) arguments should
!!  have the same values on all tasks.
!!
!!  On return, ansBlockID and ansBlockID will contain the same values on all tasks
!!  that provided the same input values for pos (and comm, if appropriate).
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_getBlkIDFromPos(pos,ansBlockID, ansProcID,comm)

#ifdef BITTREE
  use gr_interface, ONLY : gr_xyzToBlock
#else
  use Grid_data, ONLY : gr_minCellSizes, gr_globalDomain, gr_meshMe
  use Grid_interface, ONLY : Grid_getListOfBlocks,Grid_getBlkCornerID
  use Driver_interface, ONLY : Driver_abortFlash
#endif

  implicit none

  include "Flash_mpi.h"
  real, dimension(1:MDIM), intent(IN) :: pos
  integer, intent(OUT) :: ansBlockID, ansProcID
  integer, optional, intent(IN) :: comm

#ifndef BITTREE
  integer  :: blkCount
  integer,dimension(MAXBLOCKS) :: blkList
  integer :: i,j,blk, ierr, mycomm
  integer,parameter :: bufsize=2
  integer,dimension(bufsize) :: sbuf, rbuf
  integer, dimension(MDIM) :: cornerID, stride, cornerIDHigh, posInd
  logical,dimension(MDIM) :: inBlk
  logical :: onUpperBoundary(NDIM), found
  real,dimension(NDIM) :: temp
#endif

#ifdef BITTREE
  call gr_xyzToBlock(pos,ansProcID,ansBlockID)
  if (ansProcID == NONEXISTENT) ansBlockID = NONEXISTENT
#else
  if (.NOT.present(comm)) then
     call Driver_abortFlash('The specific routine Grid_getBlkIDFromPos requires a communicator argument!')
  end if
  call Grid_getListOfBlocks(LEAF,blkList,blkCount)
  ansBlockID=NONEXISTENT
  ansProcID=NONEXISTENT
  temp(1:NDIM)=(gr_globalDomain(HIGH,1:NDIM) - pos(1:NDIM))
  onUpperBoundary(:) = (temp(1:NDIM)==0)
  if(minval(temp).ge.0) then
     temp(1:NDIM)=(pos(1:NDIM)-gr_globalDomain(LOW,1:NDIM))/gr_minCellSizes(1:NDIM)
     if(minval(temp).ge.0) then
        posInd(1:NDIM)=int(temp(1:NDIM))+1
        where (onUpperBoundary(:)) posInd(1:NDIM)=posInd(1:NDIM)-1
        sbuf=0
        do i = 1,blkCount
           ansBlockID=blkList(i)
           call Grid_getBlkCornerID(ansBlockID,cornerID,stride,cornerIDHigh)
           inBlk=.true.
           inBlk(1:NDIM)=(posInd(1:NDIM).ge.cornerID(1:NDIM)).and.&
                (posInd(1:NDIM).le.cornerIDHigh(1:NDIM))
           found=inBlk(IAXIS).and.inBlk(JAXIS).and.inBlk(KAXIS)
           if(found) then
              sbuf(1)=gr_meshMe+1 !! to compensate for the fact that 0 is a valid procID
              sbuf(2)=ansBlockID
           end if
        end do
        mycomm=comm

        call MPI_AllReduce(sbuf,rbuf, 1, MPI_2INTEGER, MPI_MAXLOC,mycomm, ierr)
        ansProcID=rbuf(1)-1 !! remove the added 1
        ansBlockID=rbuf(2)
     end if
  end if
#endif
end subroutine Grid_getBlkIDFromPos
