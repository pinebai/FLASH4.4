!!****if* source/Grid/GridParticles/GridParticlesMove/Sieve/gr_ptMoveSieve
!!
!! NAME
!!  gr_ptMoveSieve
!!
!! SYNOPSIS
!!
!!  gr_ptMoveSieve(real(INOUT)    :: dataBuf(:,:),
!!                 integer(INOUT) :: localCount,
!!                 integer(IN)    :: propCount,
!!                 integer(IN)    :: maxCount,
!!                 integer(INOUT) :: numDest,
!!       optional, integet(IN)    :: matchProp)
!!  
!! DESCRIPTION 
!!  
!!  This routine moves non stationary Grid data elements to the 
!!  correct block and processor. The data elements such as particles 
!!  may have moved because of time advancement, or 
!!  after the grid has refined or derefined. Other data elements such
!!  as those associated with rays follow the  path of the ray.
!!
!!
!!  Overview of algorithm
!!
!!  * Sort the elements in the order of their associated blocks.  This makes for
!!  faster searching and processing of particles
!!
!!  * Find all the elements that belong on one of the local blocks
!!  and make sure that the appropriate block number is associated 
!!  with them.
!!
!!  * send elements in buffer to the first right neighbor. (send a
!!  dummy if no elements need to move off proc), and recv them from
!!  immediate left neighbor and add them to the elements data
!!  structure
!!
!!  * Also send gr_globalID and oldBlksProcessed arrays to neighbor
!!
!!  * repeat this process sending uprocessed elements to immediate
!!  left neighbor, then right neighbor + 2, then left neighbor + 2 and
!!  so on until all the elements have been processed.
!!
!!
!! ARGUMENTS 
!!
!!  dataBuf : List of data elements. 
!!            It is two dimensional real array, the first dimension
!!            represents each particle's properties, and second dimension is index to
!!              elements.
!!
!!  localCount : While coming in it contains the current number of elements mapped to
!!                      this processor. After all the data structure movement, the number
!!                      of local elements might change, and the new value is put back into it
!!  propCount : number of particle attributes
!!  maxCount : This is parameter determined at runtime, and is the maximum number of loca
!!             elements that a simulation expects to have. All the arrays in the elements
!!              unit are allocated based on this number
!!  numDest  : the count of data elements in the sieve at any iteration
!!  matchProp : property to be matched to find the destination of the data element
!!
!! NOTES
!!   
!!
!! SEE ALSO
!!
!!  gr_ptLocalMatch
!!
!!
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_ptMoveSieve(dataBuf,localCount,propCount,&
     maxCount, numDest, matchProp)

  use Grid_data, ONLY : gr_globalMe, gr_meshNumProcs,gr_useParticles, gr_useEnergyDeposition
  use gr_ptData, ONLY : gr_ptBlkList,gr_ptBlkCount,&
       gr_ptDestBuf,gr_ptSourceBuf,gr_ptSieveCheckFreq, gr_ptSieveFreq,&
       gr_ptBlk, gr_ptTag, gr_ptToggle, gr_ptProcDist
  use gr_ptInterface, ONLY : gr_ptLocalMatch
  use gr_ptSieveInterface, ONLY : gr_ptResetProcPair, gr_ptNextProcPair
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage
#ifdef FLASH_GRID_PARAMESH  ! If timesInLoop >= gr_meshNumProcs then
  ! something has gone wrong because the data should already be at the destination


  use tree, ONLY : lnblocks
#endif

  implicit none
  integer,intent(INOUT) :: localCount
  integer,intent(IN) :: propCount
  integer,intent(IN) :: maxCount
  integer,intent(INOUT) :: numDest
  integer, optional, intent(IN) :: matchProp
  include "Flash_mpi.h"

  real, dimension(propCount, maxCount),intent(INOUT) :: dataBuf

  

  integer,dimension(MPI_STATUS_SIZE) :: status

  logical :: stillProcessing, mustCommunicate
  real :: nothingReal
  integer :: timesInLoop,  ierr, ierr2, errLen
  integer :: numSource,i,j,maxLocalCount, src, dest
  integer :: sendCount, recvCount,sendTag,recvTag

  integer :: propToMatch

  character(len=MAX_STRING_LENGTH) :: strBuff
  character(len=10) :: numToStr
  character (len=MPI_MAX_ERROR_STRING) :: errStr

  !if we have turned off particles then return

  if(.not. (gr_useParticles .or. gr_useEnergyDeposition)) return

  if(present(matchProp)) then
     propToMatch=matchProp
  else
     propToMatch=gr_ptBlk
  endif
  

  nothingReal = -1.0  !dummy nothing particle tag to send when no elements

  gr_ptSieveFreq=gr_ptSieveCheckFreq

  stillProcessing = (numDest>0)

  !write(*,*) 'numDest=',numDest,maxCount
  !write(*,*) 'Xpos=',databuf(POSX_PART_PROP,2),'Ypos=',databuf(POSY_PART_PROP,2),'Zpos=',databuf(POSZ_PART_PROP,2) 
  if(.not.stillProcessing)gr_ptDestBuf(1,1)=nothingReal
  


  call gr_ptResetProcPair(stillProcessing,mustCommunicate)
  timesInLoop = 0 ! keep track of how many times we are in this loop.  


  
  ! If we get >= gr_meshNumProcs then
  ! something has gone wrong because the elements 
  ! have to be on one of the procs!
  do while(mustCommunicate)
     
     timesInLoop = timesInLoop + 1

     sendTag = 1
     recvTag = 1
     sendCount = max(1,numDest*propCount)
     recvCount = maxCount*propCount
     call gr_ptNextProcPair(timesInLoop, src, dest, stillProcessing, mustCommunicate)
     if(mustCommunicate) then 
        call MPI_SENDRECV(gr_ptDestBuf(1,1), sendCount, FLASH_REAL, &
             dest, sendTag, gr_ptSourceBuf(1,1), &
             recvCount, FLASH_REAL, src, & 
             recvTag, FLASH_COMM, status, ierr)
        
#ifdef DEBUG_GRIDPARTICLES
        if(ierr /= MPI_SUCCESS) then
           call MPI_ERROR_STRING(ierr,errStr,errLen,ierr2)
           call Driver_abortFlash("Grid_moveParticles SEND_RECV: " // errStr(1:errLen))
        end if
#endif
        
        !get the count of how much was returned from MPI_Recv status
        call MPI_GET_COUNT(status, FLASH_REAL, numSource, ierr)
        
#ifdef DEBUG_GRIDPARTICLES
        if(ierr /= MPI_SUCCESS) then
           call MPI_ERROR_STRING(ierr,errStr,errLen,ierr2)
           call Driver_abortFlash("Grid_moveParticles GET_COUNT: " // errStr(1:errLen))
        end if
#endif
        
        if(numSource==1)then
           numDest=0
           gr_ptDestBuf(1,1)=nothingReal
           stillProcessing=.false.
        else
           numSource=numSource/propCount
           call gr_ptLocalMatch(dataBuf,localCount,propCount,maxCount,&
                gr_ptSourceBuf,numSource,gr_ptDestBuf,numDest)
           stillProcessing=(numDest>0)
        endif
        
#ifdef DEBUG_GRIDPARTICLES
        if(ierr /= MPI_SUCCESS) then
           call MPI_ERROR_STRING(ierr,errStr,errLen,ierr2)
           call Driver_abortFlash("Grid_moveParticles MPI_ALLREDUCE: " // errStr(1:errLen))
        end if
#endif 
     end if
     
  end do !loop over distance
  
  
#ifdef DEBUG_GRIDPARTICLES
  !if(ierr /= MPI_SUCCESS) then
  !   call MPI_ERROR_STRING(ierr,errStr,errLen,ierr2)
  !   call Driver_abortFlash("Grid_moveParticles MPI_ALLREDUCE: " // errStr(1:errLen))
  !end if
  
  !if(gr_meshMe==MASTER_PE) then
  
  !   write (numToStr, '(i10.10)') maxLocalCount
  !   strBuff="maximum particles in a block"//numToStr
  !   call Logfile_stampMessage( strBuff)
  !   print*,'The maximum number of particles per processor is ',maxLocalCount
  !end if
  
#ifdef FLASH_GRID_PARAMESH  
  do i=1,localCount
     if(dataBuf(gr_ptBlk,i) > lnblocks*1.0) then
        print *, "gr_ptBlk, lnblocks ", dataBuf(gr_ptBlk,i), lnblocks
        call Driver_abortFlash("gr_ptBlk > localNumBlocks Grid_moveParticles end")
     end if
  end do
#endif
  
#endif
  
  
end subroutine gr_ptMoveSieve
