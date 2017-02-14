!!****if* source/Grid/GridParticles/GridParticlesMove/Grid_moveParticles
!!
!! NAME
!!  Grid_moveParticles
!!
!! SYNOPSIS
!!
!!  Grid_moveParticles(real(INOUT)    :: dataBuf(propCount,maxCount),
!!                     integer(IN)    :: propCount 
!!                     integer(IN)    :: maxCount,
!!                     integer(INOUT) :: localCount,
!!                     integer(IN)    :: index_list(indexCount),
!!                     integer(IN)    :: indexCount
!!                     logical(IN)    :: coords_in_blk)
!!  
!! DESCRIPTION 
!!  
!!  This routine deals with moving data associated with quantities
!! such as particles or rays, which do not have a fixed association
!! with a specific point of the domain all through the time evolution.
!! As they change their association from block to block, this routine moves
!! them to the correct block and processor. 
!! Depending upon whether the movement is due to re-gridding or not, this
!! routine calls appropriate function to move the data as needed.
!!  
!!
!!
!! ARGUMENTS 
!!
!!  dataBuf : the data structure containing the data of interest
!!              It is two dimensional real array, the first dimension
!!              represents properties associated with the data 
!!              structure, and 
!!              second dimension is index to individual elements in 
!!              the datastructure.
!!
!! propCount : number of properties for this datastructure 
!!
!!  maxCount : This is parameter determined at runtime, 
!!             and is the maximum count of elements 
!!             that a simulation expects to have. 
!!             All the arrays  are allocated based on this number
!!  localCount : While coming in it contains the current 
!!               number of elements in the data structure mapped to
!!               this processor. After all the data structure 
!!               movement, the number might change, 
!!               and the new value is put back into it
!!  index_list : The list of fields in the incoming dataBuf,
!!               the list are used to make the indices of tag, block
!!               processor and physical location of each individual
!!               data item to the GridParticles subunit so it can
!!               move the concerned data item appropriately
!!  indexCount : The count of fields included in index_list
!!
!!  coords_in_blk   : if true then this routine should not make assumptions 
!!                    about being able to determine the destination of a 
!!                    particle at the source processor. The matching 
!!                    of a particle to a block in this situation depends upon
!!                    verifying that the position coordinates of the particle
!!                    are within the bounding box of the block.
!!
!! NOTES
!!   
!!
!! SEE ALSO
!!
!!  gr_ptMoveSieve
!!
!!
!!
!!***

subroutine Grid_moveParticles(dataBuf, propCount, maxCount, localCount, &
     index_list, indexCount, coords_in_blk)
 
  use Grid_data, ONLY : gr_meshMe,gr_useParticles, gr_meshNumProcs, gr_useEnergyDeposition
  use gr_ptData, ONLY : gr_ptBlkList,gr_ptBlkCount,&
       gr_ptDestBuf,gr_ptSourceBuf,gr_ptBlk,gr_ptProc,&
       gr_ptPosx,gr_ptPosy,gr_ptPosz
  use gr_ptInterface, ONLY : gr_ptMoveSieve,gr_ptLocalMatch, &
       gr_ptMoveOffBlk, gr_ptVerifyBlock, &
       gr_ptSetIndices, gr_ptResetIndices
  use Grid_interface, ONLY : Grid_getListOfBlocks
  implicit none

#include "constants.h"
#include "Flash.h"


  integer,intent(IN) :: maxCount, propCount, indexCount
  integer,dimension(indexCount), intent(IN) :: index_list
  integer,intent(INOUT) :: localCount

  real, dimension(propCount, maxCount),intent(INOUT) :: dataBuf
  logical, intent(IN) :: coords_in_blk

  integer :: numDest,m,ii
  logical :: moveDone

  !if we have turned off particles then return

  !write(*,*) 'indexCount=',indexCount,index_list


  if(.not. (gr_useParticles .or. gr_useEnergyDeposition)) return

  !! Here get the list of blocks and populate it for the rest of the 
  !! functions that will be called, so they don't have to keep doing 
  !! it over and over again.

  call Grid_getListOfBlocks(LEAF, gr_ptBlkList, gr_ptBlkCount)

  !! Now get all the indices into the data structure setup right
  !! for the rest of the unit

  call gr_ptSetIndices(index_list,indexCount)


  dataBuf(gr_ptBlk,localCount+1:maxCount)=NONEXISTENT
  
  if(coords_in_blk) then

     numDest=0
     dataBuf(gr_ptBlk,1:localCount)=UNKNOWN

     gr_ptSourceBuf(:,1:localCount)=dataBuf(:,1:localCount)
     m=0

     call gr_ptLocalMatch(dataBuf,m,propCount,maxCount,gr_ptSourceBuf,&
          localCount, gr_ptDestBuf,numDest)
     localCount=m
     dataBuf(gr_ptBlk,localCount+1:maxCount)=NONEXISTENT
     if (gr_meshNumProcs > 1) then
#ifdef BITTREE
        call gr_ptMovePttoPt(databuf,propCount,maxCount,&
             localCount, numDest)
#else
        call gr_ptMoveSieve(dataBuf,localCount,propCount,&
             maxCount, numDest)
#endif
     end if
  else
     numDest=0
     call gr_ensureValidNeighborInfo(0)
     call gr_ptMoveOffBlk(databuf,propCount,localCount,maxCount,gr_ptDestBuf,numDest)
     call gr_ptMovePttoPt(databuf,propCount,maxCount,&
          localCount, numDest)
  end if

#ifdef DEBUG_GRIDPARTICLES
  !Check to see whether each particle has been placed on the correct block.
  call gr_ptVerifyBlock(databuf,propCount,localCount,maxCount)
#endif

  call gr_ptResetIndices(index_list,indexCount)
  
end subroutine Grid_moveParticles


