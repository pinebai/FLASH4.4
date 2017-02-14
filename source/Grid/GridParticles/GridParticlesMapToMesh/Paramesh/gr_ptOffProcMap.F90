!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/gr_ptOffProcMap
!!
!! NAME
!!  gr_ptOffProcMap
!!
!! SYNOPSIS
!!
!!  gr_ptOffProcMap(integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: srcCoords, &
!!                  integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: destCoords, &
!!                  integer,intent(IN) :: bufferSize, &
!!                  real,dimension(bufferSize),intent(INOUT) :: sendBuf, &
!!                  integer,intent(INOUT) :: sendSize, &
!!                  integer,intent(INOUT) :: sendBufPtr, &
!!                  integer,dimension(BLKID:REFLEVELDIF), intent(IN):: negh, &
!!                  integer,dimension(MDIM), intent(IN) :: neghCornerID)
!!
!! DESCRIPTION
!!
!! This routine copies the guard cell region into the one-dimensional send buffer.
!! Data is stored as metadata + smeared guard cells.  The smeared guard cells are 
!! stored in increasing memory order (rather than contiguous order) in the send buffer.
!! The coordinates inside the guard cell region are converted into an appropriate 1D 
!! offset.  Therefore, the end destination of the data, can be obtained from the position 
!! of the guard cell in the send buffer and the metadata containing the coordinates of 
!! the actual guard cell region.
!!
!! ARGUMENTS
!!               srcCoords:  The guard cell region from the source block
!!               destCoords:  The guard cell region for the destination block
!!               bufferSize:  The size of the sendBuf and recvBuf arrays
!!               sendBuf: An array to be sent as an MPI message to another processor
!!               sendSize:  The number of data elements we will send in the message
!!               sendBufPtr:  A pointer to the next free location in sendBuf
!!               negh:  An array containing information about the destination block
!!               neghCornerID:  An array containing the corner ID of the destination block
!!
!! PARAMETERS
!! 
!!***

subroutine gr_ptOffProcMap(srcCoords, destCoords, bufferSize, sendBuf, &
     sendSize, sendBufPtr, negh, neghCornerID)

  use Grid_data, ONLY : gr_meshNumProcs, gr_meshMe
  use gr_ptData, ONLY : gr_ptBuf
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_ptInterface, ONLY : gr_ptProlongSmear
  use ut_conversionInterface, ONLY : ut_ConvertToMemoryOffset

  implicit none

#include "constants.h"
#include "Flash.h"
#include "gr_ptMapToMesh.h"

  integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: srcCoords, destCoords
  integer,intent(IN) :: bufferSize
  real,dimension(bufferSize),intent(INOUT) :: sendBuf
  integer,intent(INOUT) :: sendSize
  integer,intent(INOUT) :: sendBufPtr
  integer,dimension(BLKID:REFLEVELDIF), intent(IN):: negh
  integer,dimension(MDIM), intent(IN) :: neghCornerID

  integer,dimension(MDIM) :: destStart,destEnd,srcStart,srcEnd
  integer, dimension(MDIM) :: upperLimit,elementCoord
  real, dimension(1:2, 1:2, 1:2) :: prolongedSection
  real :: val, srcToDestRefRatio
  integer :: i,j,k, i1,j1,k1,i2,j2,k2, incSrc,incDest
  integer :: numbElements, eachDataElement, ii, jj, kk, memoryOffset


  if (gr_meshNumProcs==1) then
     call Driver_abortFlash("[gr_ptOffProcMap]: Shouldn't be here when using 1 processor!")
  end if


  !-------------------------------------------------------------------------
  !                          PACK THE METADATA
  !-------------------------------------------------------------------------
  !The metadata describing the data in sendBuf has the same format
  !irrespective of the refinement level of the off-processor block.
  sendBuf(sendBufPtr+BLKID-1) = real(negh(BLKID))
  sendBuf(sendBufPtr+BLKPROC-1) = real(negh(BLKPROC))
  sendBuf(sendBufPtr+REFLEVELDIF-1) = real(negh(REFLEVELDIF))
  sendBuf(sendBufPtr+CORNERID-1) = real(neghCornerID(IAXIS))
  sendBuf(sendBufPtr+CORNERID) = real(neghCornerID(JAXIS))
  sendBuf(sendBufPtr+CORNERID+1) = real(neghCornerID(KAXIS))

  sendBuf(sendBufPtr+COORDSID-1) = real(destCoords(LOW,IAXIS))
  sendBuf(sendBufPtr+COORDSID) = real(destCoords(HIGH,IAXIS))
  sendBuf(sendBufPtr+COORDSID+1) = real(destCoords(LOW,JAXIS))
  sendBuf(sendBufPtr+COORDSID+2) = real(destCoords(HIGH,JAXIS))
  sendBuf(sendBufPtr+COORDSID+3) = real(destCoords(LOW,KAXIS))
  sendBuf(sendBufPtr+COORDSID+4) = real(destCoords(HIGH,KAXIS))

  sendBufPtr = sendBufPtr + SIZE_HEADER


  !-------------------------------------------------------------------------
  !                             PACK THE DATA
  !-------------------------------------------------------------------------
  if(negh(REFLEVELDIF)==0) then
     !write(*,*)"Neighbor at the same refinement, and on another processor."
     do k = 0, (srcCoords(HIGH,KAXIS) - srcCoords(LOW,KAXIS)), 1
        do j = 0, (srcCoords(HIGH,JAXIS) - srcCoords(LOW,JAXIS)), 1
           do i = 0, (srcCoords(HIGH,IAXIS) - srcCoords(LOW,IAXIS)), 1

              sendBuf(sendBufPtr) = gr_ptBuf(srcCoords(LOW,IAXIS)+i, &
                   srcCoords(LOW,JAXIS)+j, srcCoords(LOW,KAXIS)+k)
              sendBufPtr = sendBufPtr + 1

           end do
        end do
     end do

  else

     !First copy the coordinates into smaller arrays which
     !are easier to work with.
     srcStart(:) = srcCoords(LOW,:)
     srcEnd(:) = srcCoords(HIGH,:)
     destStart(:) = destCoords(LOW,:)
     destEnd(:) = destCoords(HIGH,:)

     srcToDestRefRatio = 2.0**(-negh(REFLEVELDIF))
     incSrc=int(max(1.0,srcToDestRefRatio))
     incDest=int(max(1.0,(1.0/srcToDestRefRatio)))


     !Determine the appropriate upper limit:
     if(negh(REFLEVELDIF) == 1) then
        upperLimit(:) = srcEnd(1:MDIM)-srcStart(1:MDIM)+1

        numbElements = (destCoords(HIGH,IAXIS)-destCoords(LOW,IAXIS)+1) * &
             (destCoords(HIGH,JAXIS)-destCoords(LOW,JAXIS)+1) * &
             (destCoords(HIGH,KAXIS)-destCoords(LOW,KAXIS)+1)

        !Zero the sendBuf because the inner loop over ii, jj, kk will only 
        !write to certain memory offsets.
        do eachDataElement = sendBufPtr, sendBufPtr + numbElements - 1, 1
           sendBuf(eachDataElement) = 0.0
        end do

     else if(negh(REFLEVELDIF) == -1) then
        upperLimit(:) = destEnd(1:MDIM)-destStart(1:MDIM)+1
     else
        call Driver_abortFlash("[gr_ptOffProcMap]: Unrecognised refinement difference between source and dest blocks")
     end if



     do k = 1, upperLimit(KAXIS)

        k1 = srcStart(KAXIS) + ((k-1)*incSrc)
        k2 = destStart(KAXIS) + ((k-1)*incDest)

        do j = 1, upperLimit(JAXIS)

           j1 = srcStart(JAXIS) + ((j-1)*incSrc)
           j2 = destStart(JAXIS) + ((j-1)*incDest)

           do i = 1, upperLimit(IAXIS)

              i1 = srcStart(IAXIS) + ((i-1)*incSrc)
              i2 = destStart(IAXIS) + ((i-1)*incDest)


              if(negh(REFLEVELDIF) == 1) then

                 val = gr_ptBuf(i1,j1,k1)
                 if(val /= 0.0) then

#ifdef DEBUG_GRIDMAPPARTICLES_VERBOSE
                    write(*,"(a,i6,a,1X,i3,1X,i3,1X,i3,a,1X,i3,a,1X,i3,1X,i3,a,1X,i3,1X,i3,a,1X,i3)") &
                         " Processor", gr_meshMe, " prolonging from ", i1, j1, k1, &
                         " into range ",i2,":",i2+1, j2,":",j2+K2D, k2,":",k2+K3D
#endif

                    call gr_ptProlongSmear(i1,j1,k1,prolongedSection)

                    !Copy prolonged section into sendBuf.
                    do kk = k2, k2+K3D
                       do jj = j2, j2+K2D
                          do ii = i2, i2+1

                             elementCoord(IAXIS) = ii
                             elementCoord(JAXIS) = jj
                             elementCoord(KAXIS) = kk

                             !Returns the memory offset into the array at element ii,jj,kk
                             call ut_ConvertToMemoryOffset(MDIM, elementCoord, destStart, &
                                  destEnd, memoryOffset)

                             sendBuf(sendBufPtr + memoryOffset) = & 
                                  prolongedSection( (ii-i2+1), (jj-j2+1), (kk-k2+1) )

#ifdef DEBUG_GRIDMAPPARTICLES_VERBOSE
                             write(*,"(a,i6,a,ES12.5,a,1X,i3,1X,i3,1X,i3,a,i6)") &
                                  " Processor", gr_meshMe, " sending prolonged data=", &
                                  prolongedSection( (ii-i2+1), (jj-j2+1), (kk-k2+1) ), &
                                  " destined for element", ii, jj, kk, " in block", negh(BLKID)
#endif

                          end do
                       end do
                    end do

                 end if


              else if(negh(REFLEVELDIF) == -1) then

#ifdef DEBUG_GRIDMAPPARTICLES_VERBOSE
                 write(*,"(a,i6,a,1X,i3,a,1X,i3,1X,i3,a,1X,i3,1X,i3,a,1X,i3,a,1X,i3,1X,i3,1X,i3)") &
                      " Processor", gr_meshMe, " restricting from range ", &
                      i1,":",i1+incSrc-1,j1,":",j1+((incSrc-1)*K2D),k1,":",k1+((incSrc-1)*K3D), &
                      " into ",i2,j2,k2
#endif

                 !Divide by 2.0**NDIM because the density from the smaller source cell 
                 !is distributed into a larger cell space.
                 sendBuf(sendBufPtr) = &
                      (sum( gr_ptBuf(i1:i1+incSrc-1,j1:j1+((incSrc-1)*K2D), k1:k1+((incSrc-1)*K3D)) ) / &
                      2.0**NDIM)

#ifdef DEBUG_GRIDMAPPARTICLES_VERBOSE
                 write(*,"(a,i6,a,ES12.5,a,1X,i3,1X,i3,1X,i3,a,i6)") &
                      " Processor", gr_meshMe, " sending restricted data=", &
                      sendBuf(sendBufPtr), " destined for element", i2, j2, k2, " in block", negh(BLKID)
#endif

                 sendBufPtr = sendBufPtr + 1

              end if

           end do
        end do
     end do


     !We must update sendBufPtr.
     if(negh(REFLEVELDIF) == 1) then
        sendBufPtr = sendBufPtr + numbElements
     end if

  end if

  sendSize = sendBufPtr - 1  !This is no. of REALS to be sent in the message.

  return

end subroutine gr_ptOffProcMap
