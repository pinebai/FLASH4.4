!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/gr_ptPackUnpackData
!!
!! NAME
!!  gr_ptPackUnpackData
!!
!! SYNOPSIS
!!
!!  gr_ptPackUnpackData(integer,intent(IN) :: varGrid, &
!!                      integer,intent(IN) :: bufferSize, &
!!                      real,dimension(bufferSize),intent(INOUT) :: sendBuf, &
!!                      integer,intent(INOUT) :: sendCount, &
!!                      real,dimension(bufferSize),intent(INOUT) :: recvBuf, &
!!                      integer,intent(IN) :: recvCount)
!!
!! DESCRIPTION
!!
!!  This routine determines whether the data in the receive buffer (recvBuf) is 
!!  meant to reside on this processor.  The receive buffer contains data and 
!!  metadata.  The metadata may indicate the destination processor directly, or
!!  contain the corner ID of the destination block.  If the metadata processor ID 
!!  contains NONEXISTENT, then this processor will search each of its blocks 
!!  and attempt to match against the metadata corner ID, and if found, update the 
!!  processor ID.  The data is unpacked if the processor IDs match.  The 
!!  relative refinement of the source block guard cell region to the destination block 
!!  does not matter, as the data was packed at the refinement level of the destination
!!  block.  Therefore, data is always unpacked the same way.  If the data is meant 
!!  for another processor, then it is copied directly into the send buffer (sendbuf).
!!
!! ARGUMENTS
!!               varGrid:   Index of gridded variable to receive interpolated
!!                              quantity
!!               bufferSize:    Size of the recvBuf and sendBuf arrays
!!               sendBuf:       The array containing data to be sent to another 
!!                              processor
!!               sendCount:     The number of data elements we have repacked 
!!                              into sendBuf
!!               recvBuf:       The array containing data we just received from
!!                              another processor
!!               recvCount:     The number of data elements we have received 
!!                              from the incoming message, which is placed 
!!                              in recvBuf 
!!
!! 
!!***

!!REORDER(4): solnData

subroutine gr_ptPackUnpackData(varGrid, bufferSize, sendBuf, sendCount, recvBuf, recvCount)

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkIndexLimits
  use gr_ptInterface, ONLY : gr_ptSearchBlk, gr_ptParseMetadata
  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs,gr_meshComm
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "gr_ptMapToMesh.h"

  integer,intent(IN) :: varGrid
  integer,intent(IN) :: bufferSize
  real,dimension(bufferSize),intent(OUT) :: sendBuf
  integer,intent(OUT) :: sendCount
  real,dimension(bufferSize),intent(IN) :: recvBuf
  integer,intent(IN) :: recvCount

  integer,dimension(MDIM) :: neghCornerID
  integer :: iCoord, jCoord, kCoord, iCell, iCopy
  integer :: headerPtr, dataPtr, sendBufPtr
  integer :: numbElements, i, j, k    
  real :: value
  integer, dimension(LOW:HIGH,MDIM) :: sectionCoords
  integer, dimension(MDIM) :: negh
  real,dimension(:,:,:,:),pointer :: solnData


  sendCount = 0

  if(recvCount > 1) then
     !If we receive some data, we must scan through it and pick out the relevant elements.
     headerPtr = 1
     sendBufPtr = 1

     EachHeader: do

        !Extract the data description.
        call gr_ptParseMetadata(bufferSize, recvBuf, headerPtr, negh, &
             neghCornerID, sectionCoords, numbElements)

        if(negh(BLKPROC) == NONEXISTENT) then
           call gr_ptSearchBlk(neghCornerID,negh)           
        end if


        if(negh(BLKPROC) == gr_meshMe) then

           !print *, "Processor", gr_meshMe, "data match... extracting data to block", negh(BLKID)
           call Grid_getBlkPtr(negh(BLKID), solnData, CENTER)

           dataPtr = headerPtr + SIZE_HEADER        

           do k = 0, (sectionCoords(HIGH,KAXIS) - sectionCoords(LOW,KAXIS)), 1
              do j = 0, (sectionCoords(HIGH,JAXIS) - sectionCoords(LOW,JAXIS)), 1
                 do i = 0, (sectionCoords(HIGH,IAXIS) - sectionCoords(LOW,IAXIS)), 1

                    value = recvBuf(dataPtr)
                    iCoord = sectionCoords(LOW,IAXIS)+i
                    jCoord = sectionCoords(LOW,JAXIS)+j
                    kCoord = sectionCoords(LOW,KAXIS)+k


#ifdef DEBUG_GRIDMAPPARTICLES_VERBOSE
                    if(value/=0.0) then
                       if(negh(REFLEVELDIF)==1) then
                          write(*,"(a,i6,a,ES12.5,a,1X,i3,1X,i3,1X,i3,a,i6)") &
                               " Processor", gr_meshMe, " received prolonged data=", value, &
                               " and now storing at", iCoord, jCoord, kCoord, " in block", negh(BLKID)
                       
                       else if(negh(REFLEVELDIF)==-1) then
                          write(*,"(a,i6,a,ES12.5,a,1X,i3,1X,i3,1X,i3,a,i6)") &
                               " Processor", gr_meshMe, " received restricted data=", value, &
                               " and now storing at", iCoord, jCoord, kCoord, " in block", negh(BLKID)
                       end if
                    end if
#endif


                    !Assert that the i, j, k coordinate maps to a central section.
                    if( (iCoord < NGUARD+1) .or. (iCoord > NXB + NGUARD) ) then
                       print *, "[gr_ptPackUnpackData]: iCoord is out of range = ", iCoord
                       call Driver_abortFlash("[gr_ptPackUnpackData]: Error, i not mapping to central region")
                    else if ( (K2D == 1) .and. ((jCoord < NGUARD+1) .or. (jCoord > NYB + NGUARD)) ) then
                       print *, "[gr_ptPackUnpackData]: jCoord is out of range = ", jCoord
                       call Driver_abortFlash("[gr_ptPackUnpackData]: Error, j not mapping to central region")
                    else if ( (K3D == 1) .and. ((kCoord < NGUARD+1) .or. (kCoord > NZB + NGUARD)) ) then
                       print *, "[gr_ptPackUnpackData]: kCoord is out of range = ", kCoord
                       call Driver_abortFlash("[gr_ptPackUnpackData]: Error, k not mapping to central region")
                    end if


                    solnData(varGrid,iCoord,jCoord,kCoord) = solnData(varGrid,iCoord,jCoord,kCoord) + value

                    dataPtr = dataPtr + 1

                 end do
              end do
           end do

           call Grid_releaseBlkPtr(negh(BLKID),solnData,CENTER)


        else
           !Not my data, so repack it, and pass it onto the next processor.
           !print *, "Processor", gr_meshMe, "non data match... repacking data"
           do iCopy = 0, SIZE_HEADER + numbElements - 1, 1
              sendBuf(sendBufPtr) = recvBuf(headerPtr+iCopy)
              sendBufPtr = sendBufPtr + 1
           end do

           !sendCount is total data repacked.
           sendCount = sendCount + SIZE_HEADER + numbElements

        end if

        !Skip over the data and move onto the next header description.
        headerPtr = headerPtr + SIZE_HEADER + numbElements

        !All data has been checked, so break from the extraction loop.
        if((headerPtr-1) == recvCount) then
           exit 
        else if((headerPtr-1) > recvCount) then
           call Driver_abortFlash("[gr_ptPackUnpackData] Metadata loop will not exit cleanly!") 
        end if

     end do EachHeader

  end if

end subroutine gr_ptPackUnpackData
