!!****if* source/Grid/GridParticles/GridParticlesMove/paramesh/VirtualParticles/gr_ptMoveOffBlk
!!
!! NAME
!!
!!  gr_ptMoveOffBlk
!!
!! SYNOPSIS
!!
!!  gr_ptMoveOffBlk(real(INOUT)   :: dataBuf(propCount,localCount),
!!                 integer(in)     :: propCount,
!!                 integer(INOUT) :: localCount,
!!                 real(INOUT)    :: destBuf(propCount,localCount),
!!                 integer(INOUT) :: numDest)
!!
!! DESCRIPTION
!!     
!!    This routine is used in moving the non stationaly data elements 
!!    associated with structures like particles and ray, when a data element
!!    moves off a block without re-gridding. Here every element currently 
!!    on the processor is examined to see if it still belongs to the same block.
!!    If it does not, it is further examimned to see if it has moved out of the physical boundary.
!!    If is out of physical boundary, it may either leave the domain, stay on the same block
!!    or be moved  to destBuf, which holds elements to be passed to the next processor, depending
!!    on the boundary conditions. If it is still in the physical domain, it may have
!!    moved to another block on the same processor, in which case only its BLK
!!    needs to change, otherwise it is moved to destBuf.
!!
!! ARGUMENTS
!!
!!     dataBuf -           A 2 dimensional array of elements and their property
!!                           values
!!     propCount - number of element attributes
!!     localCount -   The number of elements in the element array that
!!                           are stored on this processor.
!!     globalCount - The maximum number of elements in the element array that
!!                           can be stored on this processor.    
!!
!!     destBuf   -           temporary storage for elements that need to move off processor
!!     numDest   -           number of elements in destBuf
!! 
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_PARTICLES
#endif
subroutine gr_ptMoveOffBlk(dataBuf,propCount,localCount,globalCount,destBuf,numDest)
#include "constants.h"
#include "Flash.h"
  use gr_ptData, ONLY : gr_ptBlkList, gr_ptBlk, gr_ptProc,&
       gr_ptPosx,gr_ptPosy,gr_ptPosz, gr_ptKeepLostParticles
  use Grid_data, ONLY : gr_meshMe
  use Grid_interface, ONLY : Grid_getListOfBlocks
  use Driver_data,ONLY: dr_nstep
  use gr_interface, ONLY:  gr_findNeghID, gr_findBlock

  implicit none
  integer, intent(IN) :: propCount,globalCount
  integer,intent(INOUT)::localCount,numDest
  real,dimension(propCount,globalCount),intent(INOUT)::dataBuf,destBuf
  integer :: blockID,lostElements,blkCount,newBlockID
  real,dimension(LOW:HIGH,MDIM) :: bndBox
  real,dimension(MDIM) :: pos
  integer,dimension(MDIM) :: Negh
  integer,dimension(BLKNO:PROCNO) :: neghID
  logical :: moved
  integer i,k,j, destCount
  real,allocatable,dimension(:,:) :: destParticles
  logical :: newBlkID=.true.
  destCount=0
  numDest = 0
  lostElements=0

  i=1
  j=localCount

  allocate(destParticles(propCount,NDIM**3+1))
  
  do k=1,localCount
     blockID=int(dataBuf(gr_ptBlk,i))  !! Find the ID of the block the element

     !! belonged to before time advance
     if((blockID/=NONEXISTENT).and.(blockID/=LOST)) then
        call gr_ptVPGenerate(dataBuf(:,i),propCount,destCount,destParticles,&
             blockID,newBlkID)
        
        if((destCount == -1).and.(.not.gr_ptKeepLostParticles))then
           dataBuf(:,i)=dataBuf(:,j)
           j=j-1
        else
           i=i+1
           if(destCount>0) then
              destBuf(:,numDest+1:numDest+destCount)=destParticles(:,1:destCount)
              numDest=numDest+destCount
           end if
        end if
     end if
  end do
  deallocate(destParticles)
end subroutine gr_ptMoveOffBlk
