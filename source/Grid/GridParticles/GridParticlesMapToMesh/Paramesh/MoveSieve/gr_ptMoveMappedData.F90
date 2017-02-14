!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/MoveSieve/gr_ptMoveMappedData
!!
!! NAME
!!  gr_ptMoveMappedData
!!
!! SYNOPSIS
!!
!!  gr_ptMoveMappedData(integer,intent(IN) :: varGrid, &
!!                       integer,intent(IN) :: bufferSize, &
!!                       real,dimension(bufferSize),intent(INOUT) :: sendBuf, &
!!                       integer,intent(INOUT) :: sendCount, &
!!                       real,dimension(bufferSize),intent(INOUT) :: recvBuf)
!!
!! DESCRIPTION
!!
!! Routine which manages the communication of smeared grid cells between processors. 
!! The smeared grid cells are stored in a send buffer, which is passed between
!! processors.  A general communication approach is taken because we may not know 
!! the destination processor until the metadata corner ID is analysed.  As such, 
!! an approach is taken which (on average) minimises the number of times the send buffer 
!! is passed beween processors in order for the smeared grid cells reach the correct processor.
!! In general, if a nearest neighbor block exists on another processor, then that processor 
!! will be a close neighbor of the current processor.  This is because the blocks 
!! are ordered according to a Morton-Space filling curve, and as such, there is a 
!! high degree of spatial locality.  
!!
!! Each processor passes its send buffer back and forth, starting with its nearest neighbors, and 
!! then gradually working to its furthest away neighbors.  This exploits the spatial 
!! locality, and means, in general, the size of send buffer decreases rapidly with each 
!! passing round.
!! 
!! ARGUMENTS
!!               varGrid:   Index of gridded variable to receive interpolated
!!                              quantity
!!               bufferSize:  The size of the sendBuf and recvBuf arrays
!!               sendBuf:  An array used to store data intended for another processor
!!               sendCount:  The number of data elementes to be sent to another
!!                           processor
!!               recvBuf:  An array containing the data just receieved from another 
!!                         processor
!!
!! PARAMETERS
!! 
!!***

subroutine gr_ptMoveMappedData(varGrid,bufferSize,sendBuf,sendCount,recvBuf)

  use Grid_data, ONLY : gr_globalMe, gr_meshNumProcs, gr_meshComm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY: Logfile_stampMessage
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_ptInterface, ONLY : gr_ptPackUnpackData, gr_ptDumpState

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "gr_ptMapToMesh.h"

  integer,intent(IN) :: varGrid
  integer,intent(IN) :: bufferSize
  real,dimension(bufferSize),intent(INOUT) :: sendBuf
  integer,intent(INOUT) :: sendCount
  real,dimension(bufferSize),intent(INOUT) :: recvBuf

  logical :: mustCommunicate
  integer :: procDist, toggle, timesInLoop, dest, src, ierr
  integer :: recvCount,sendTag,recvTag, ierr2, errLen
  integer,dimension(MPI_STATUS_SIZE) :: status
  logical :: stillProcessing
  character (len=MPI_MAX_ERROR_STRING) :: errStr

  if(gr_meshNumProcs==1) then
     if(sendCount/=0) call Driver_abortFlash("gr_ptMoveMappedData: sendSize>0 for 1 proc")
  else

     !All processors must participate!!!
     if(sendCount == 0) then
        stillProcessing=.false.
        sendBuf(1) = real(NONEXISTENT)
        sendCount = 1
     else
        stillProcessing=.true.
     end if
     !print *, "Processor", gr_meshMe, "sendBuf contains", sendCount, "elements"


     mustCommunicate = .false.
     call MPI_ALLREDUCE(stillProcessing, mustCommunicate, 1, MPI_LOGICAL, MPI_LOR, gr_meshComm, ierr)  

#ifdef DEBUG_GRIDMAPPARTICLES
     if(ierr /= MPI_SUCCESS) then
        call MPI_ERROR_STRING(ierr,errStr,errLen,ierr2)
        call Driver_abortFlash("gr_ptMoveMappedData MPI_ALLREDUCE: " // errStr(1:errLen))
     end if
#endif

#ifdef DEBUG_GRIDMAPPARTICLES_VERBOSE
     print *, "Processor", gr_meshMe, "round=", 0, "still processing=", stillprocessing, & 
     ", anyone processing=", mustCommunicate
#endif


     procDist = 1
     toggle = 1
     timesInLoop = 0 ! keep track of how many times we are in this loop.  
     ! If we get >= gr_meshNumProcs then
     ! something has gone wrong because the particles 
     ! have to be on one of the procs!
     do while(mustCommunicate)

        timesInLoop = timesInLoop + 1

        if(timesInLoop >= gr_meshNumProcs) then
           write(*,*)"[gr_ptMoveMappedData]: timesInLoop >= gr_meshNumProcs"

           !sendBuf contains grid points that we will pass onto the next processor.  
           !However, as we have already passed around these grid points to all 
           !processors, we know they will never be extracted from sendBuf.
           !Therefore, print the metadata describing the grid points to file 
           !to figure out what went wrong.
           call gr_ptDumpState(bufferSize, sendBuf, sendCount)

           !Wait for all processors to finish dumping their state 
           !to file.  We can use a barrier safely here because all processors 
           !have the same value in "mustCommunicate" variable.
           call MPI_Barrier(gr_meshComm, ierr)
           call Logfile_stampMessage("[gr_ptMoveMappedData]: timesInLoop >= gr_meshNumProcs")
           call Driver_abortFlash &
                ("[gr_ptMoveMappedData]: particles should have been processed by now (see log files)")
        end if


        dest = gr_globalMe + procDist*toggle
        src = gr_globalMe - procDist*toggle


        !catch ending conditions
        if(dest >= gr_meshNumProcs) then
           dest = dest - gr_meshNumProcs
        end if

        if(src < 0) then
           src = src + gr_meshNumProcs
        end if

        if(src >= gr_meshNumProcs) then
           src = src - gr_meshNumProcs
        end if

        if(dest < 0) then
           dest = dest + gr_meshNumProcs
        end if

#ifdef DEBUG_GRIDMAPPARTICLES
        if(dest < 0) then
           call Driver_abortFlash("gr_ptMoveMappedData: dest < 0")
        end if

        if(dest >= gr_meshNumProcs) then
           call Driver_abortFlash("gr_ptMoveMappedData: dest >= gr_meshNumProcs")
        end if

        if(src < 0) then
           call Driver_abortFlash("gr_ptMoveMappedData: src < 0")
        end if

        if(src >= gr_meshNumProcs) then
           call Driver_abortFlash("gr_ptMoveMappedData: scr >= gr_meshNumProcs")
        end if

        if(src == gr_globalMe) then
           call Driver_abortFlash("src == gr_globalMe, and we still have unmatched particles, bad")
        end if

        if(dest == gr_globalMe) then
           call Driver_abortFlash("dest == gr_globalMe, and we still have unmatched particles, bad")
        end if
#endif

        sendTag = 1
        recvTag = 1
        call Timers_start("sendrecv")

        !Use the tried and tested communication technique once more, e.g.
        !print *, "Processor", gr_meshMe, "sending", sendCount, "REALS to processor", dest
        call MPI_SENDRECV(sendBuf, sendCount, FLASH_REAL, &
             dest, sendTag, recvBuf, bufferSize, FLASH_REAL, src, &
             recvTag, gr_meshComm, status, ierr)

        call Timers_stop("sendrecv")

#ifdef DEBUG_GRIDMAPPARTICLES
        if(ierr /= MPI_SUCCESS) then
           call MPI_ERROR_STRING(ierr,errStr,errLen,ierr2)
           call Driver_abortFlash("gr_ptMoveMappedData SEND_RECV: " // errStr(1:errLen))
        end if
#endif

        !get the count of how much was returned from MPI_Recv status
        call MPI_GET_COUNT(status, FLASH_REAL, recvCount, ierr)

#ifdef DEBUG_GRIDMAPPARTICLES
        if(ierr /= MPI_SUCCESS) then
           call MPI_ERROR_STRING(ierr,errStr,errLen,ierr2)
           call Driver_abortFlash("gr_ptMoveMappedData GET_COUNT: " // errStr(1:errLen))
        end if
#endif

        !print *, "Processor", gr_meshMe, "received", recvCount, "REALS from processor", src

        !Extract the data from recvBuf, and repack other processor's data into sendBuf.
        call gr_ptPackUnpackData(varGrid, bufferSize, sendBuf, sendCount, recvBuf, recvCount)

        !All processors must participate!!!
        !If we receive nothing, then send nothing!
        if((sendCount == 0).or.(recvCount == 1)) then
           stillProcessing=.false.
           sendBuf(1) = real(NONEXISTENT)
           sendCount = 1
        else
           stillProcessing=.true.
        end if

        mustCommunicate = .false.
        call MPI_ALLREDUCE(stillProcessing, mustCommunicate, 1, MPI_LOGICAL, &
             MPI_LOR, gr_meshComm, ierr) 


#ifdef DEBUG_GRIDMAPPARTICLES
        if(ierr /= MPI_SUCCESS) then
           call MPI_ERROR_STRING(ierr,errStr,errLen,ierr2)
           call Driver_abortFlash("gr_ptMoveMappedData MPI_ALLREDUCE: " // errStr(1:errLen))
        end if
#endif

#ifdef DEBUG_GRIDMAPPARTICLES_VERBOSE
        print *, "Processor", gr_meshMe, "round=", timesInLoop, "still processing=", stillprocessing, & 
             ", anyone processing=", mustCommunicate        
#endif

        !set dest and src for next time through loop
        toggle = toggle * (-1)     
        procDist = procDist + 1
     end do !loop over distance

  end if

end subroutine gr_ptMoveMappedData
