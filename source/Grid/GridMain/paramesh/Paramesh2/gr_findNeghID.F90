!!****if* source/Grid/GridMain/paramesh/Paramesh2/gr_findNeghID
!!
!! NAME
!!
!!  gr_findNeghID
!!
!! SYNOPSIS
!!
!!  call gr_findNeghID(integer(IN)  :: blockid,
!!                     real(IN)     :: pos(MDIM),
!!                     integer(IN)  :: negh(MDIM),
!!                     integer(OUT) :: neghid(BLKNO:PROCNO))
!!
!! DESCRIPTION
!!
!!   Given the physical coordinates of a point outside the current
!!   block, this routine finds the processor number and blockID
!!   within that processor for the neighboring block that contains 
!!   the point.
!!
!! ARGUMENTS
!!
!!   blockid : ID of block in current processor
!!
!!   pos :     coordinates of the point of interest
!!
!!   negh :    the location of the neighbor with respect to the current block,
!!             in other words specification on which face/edge/point is common
!!             between the current block and neighbor of interest. For example
!!             Negh(1:MDIM)=LEFT_EDGE indicates that the lowest left hand corner
!!             of the current block is the same as the highest right end corner
!!             of the neighbor. Similarly Negh(IAXIS)=LEFT_EDGE, Negh(JAXIS:KAXIS)
!!             = CENTER implies that the left face of current block is common
!!             with the right face of the neighbor
!!
!!   neghid : identity of the neighbor, the second number is the processor
!!            number where the neighbor is located, and the first number 
!!            is the blocknumber within the processor
!!
!!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif
subroutine gr_findNeghID(blockID,pos,Negh,neghID)
#include "constants.h"
#include "Flash.h"

  use tree, ONLY : neigh,nodetype, child
  use Grid_data, ONLY : gr_globalDomain, gr_meshMe
  use gr_interface, ONLY:  gr_findWhichChild
  use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_getBlkBC

  implicit none
  
  integer,intent(IN) :: blockID
  real,dimension(MDIM), intent(IN) :: pos
  integer,dimension(MDIM),intent(IN) :: Negh
  integer,dimension(BLKNO:PROCNO),intent(OUT) :: neghID

  real,dimension(LOW:HIGH,MDIM) :: bndBox
  real, dimension(MDIM) :: wpos, deltaDomain
  integer, dimension(LOW:HIGH, MDIM) :: faces, ignoreMe
  integer :: face,i,k,tempProcID, tempBlockID, whichChild, eachAxis
  logical :: onProc

  k=0
  tempBlockID=blockID

  !!  Assume that the neigh exists and is on the same processor
  onProc = .true.
  do i = 1, NDIM
     if(onProc) then 
        if(Negh(i)/=CENTER) then !! if actually a neigh along this dimension

           !! Negh values can be 1 or 3, whereas negh are 1 or 2 along any axis
           !! k is there to offset the dimension; it is 0 for x, 2 for y and 4 for z
           face=k+min(HIGH,Negh(i))

           !! find out if the negh along this dimension exists and is on mype
           tempProcID=neigh(PROCNO,face,tempBlockID)

#ifdef DEBUG_GRID
           if(tempProcID >= 0 .AND. neigh(BLKNO,face,tempBlockID) .LE. 0)then

999           format ('gr_findNeghID:',I5,' neigh(:,',I1,',',I3,')=',3(I5))
              print 999,gr_meshMe,face,tempBlockID,neigh(:,face,tempBlockID)

           end if
#endif
           if(tempProcID==gr_meshMe .AND. neigh(BLKNO,face,tempBlockID) .GT. 0)then

              !! if neigh exists and is on mype then continue processing
              tempBlockID=neigh(BLKNO,face,tempBlockID)

           else

              !! if the neigh along this dimension either doesn't exist (not on the
              !! same refinement level) or off processor, we are done with processing
              !! Setting onProc to false ensure that the processing is bypassed for
              !! subsequent dimensions
              onProc=.false.
              tempBlockId=NONEXISTENT
           end if
        end if
        k=k+2
     end if
  end do

  if(tempBlockID==NONEXISTENT) then
     neghID(PROCNO)=UNKNOWN
     neghID(BLKNO)=NONEXISTENT
  else
     if(nodetype(tempBlockID)==LEAF) then
        neghID(PROCNO)=tempProcID
        neghID(BLKNO)=tempBlockID
     else
        !Neighbor block is a PARENT.

        !When we have a neighbor on a periodic boundary adjust the position 
        !"pos" to its wrapped around position "wpos".
        !Also move the box to the other side of the domain so that 
        !"midPoint" will be correct in gr_findWhichChild.
        !These changes are both temporary and are local to this subroutine.
        !--------------------------------------------------------------------------
        call Grid_getBlkBC(blockID, ignoreMe, faces)
        call Grid_getBlkBoundBox(blockID, bndBox)
        wpos = pos
        
        deltaDomain(1:MDIM) = &
             (gr_globalDomain(HIGH,1:MDIM) - gr_globalDomain(LOW,1:MDIM))
        
        do eachAxis = 1, NDIM
           if ( (negh(eachAxis) == LEFT_EDGE) .and. & 
                (faces(LOW,eachAxis) == PERIODIC) ) then
              
              wpos(eachAxis) = pos(eachAxis) + deltaDomain(eachAxis)                
              bndBox(LOW:HIGH,eachAxis) = &
                   bndBox(LOW:HIGH,eachAxis) + deltaDomain(eachAxis)
              
           else if ( (negh(eachAxis) == RIGHT_EDGE) .and. & 
                (faces(HIGH,eachAxis) == PERIODIC) ) then
              
              wpos(eachAxis) = pos(eachAxis) - deltaDomain(eachAxis)
              bndBox(LOW:HIGH,eachAxis) = &
                   bndBox(LOW:HIGH,eachAxis) - deltaDomain(eachAxis)
           end if
        end do
        !--------------------------------------------------------------------------

        !Pass the bounding box of the source block and the relative coordinates 
        !of the neighbor in "Negh".  We get back whichChild which is a child 
        !coordinate in the neighbor block.
        call gr_findWhichChild(wpos,bndBox,Negh,whichChild)

        !tempBlock ID holds the ID of the parent block (i.e. source block's 
        !neighbor).
        tempProcID=child(PROCNO,whichChild,tempBlockID)
        if(tempProcID==gr_meshMe) then
           neghID(PROCNO)=tempProcID
           neghID(BLKNO)=child(BLKNO,whichChild,tempBlockID)
        else
           neghID(PROCNO)=UNKNOWN
           neghID(BLKNO)=NONEXISTENT
        end if
     end if
  end if
  
end subroutine gr_findNeghID
