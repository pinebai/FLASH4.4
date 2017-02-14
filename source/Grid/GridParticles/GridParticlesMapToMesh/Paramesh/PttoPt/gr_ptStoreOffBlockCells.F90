!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/PttoPt/gr_ptStoreOffBlockCells
!!
!! NAME
!!  gr_ptStoreOffBlockCells
!!
!! SYNOPSIS
!!
!!   gr_ptStoreOffBlockCells (integer,dimension(MAXBLOCKS), intent(IN) :: particlesPerBlk
!!                            integer,dimension(MAXBLOCKS), intent(IN) :: blockList
!!                            integer,intent(IN) :: blockCount
!!                            integer,dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
!!                            integer,dimension(MDIM), intent(IN) :: blkSize
!!                            integer,dimension(MDIM), intent(IN) :: guard
!!                            integer, intent(OUT) :: BufferSize)
!!
!! DESCRIPTION
!!
!! This subroutine has two purposes.  The first is to calculate the maximum 
!! size of the send / receive buffer on each processor, and the second is to 
!! store useful information in a user defined type (UDT).  The first part is non-trivial 
!! because the send / receive buffer consists of metadata and smeared
!! mesh points assigned to off-processor neighboring blocks.  Unfortunately, the only way to 
!! determine the number of off-processor mesh points is to loop over each block's guard cell 
!! regions and count each mesh point.  As this is rather time-consuming, we store useful
!! information about the computational domain in a UDT (gr_ptDomain) whilst counting.
!! 
!! An initial estimate of the size of the send / receive buffer is determined locally.  
!! This value is compared across all processors, and the maximum size is used for each 
!! each processor's send / receive buffer. 
!!
!! ARGUMENTS
!!
!!                          particlesPerBlk: Number of particles residing on each block.
!!                          blockList: List of all leaf blocks existing on this processor.
!!                          blockCount: Number of leaf blocks existing on this processor.
!!                          blkLimitsGC: Upper and lower indicies of source block including guard cells.
!!                          blkSize: Size of the source block (same for each block).
!!                          guard: Number of guard cells for the source block.
!!                          BufferSize: The size of the send / receive buffer.
!! 
!!***

subroutine gr_ptStoreOffBlockCells(particlesPerBlk, blockList, blockCount, blkLimitsGC, blkSize, guard, BufferSize)

  use gr_ptMapData, ONLY : gr_ptDomain, gr_ptRecvSpecifier, gr_ptRecvSpecifierTmp, &
       gr_ptRecvTotal, gr_ptRecvTotalTmp, gr_ptNumMessagesToSend
  use Grid_interface, ONLY : Grid_getBlkCornerID, Grid_getBlkBoundBox, Grid_getDeltas
  use gr_ptInterface, ONLY : gr_ptFindNegh, gr_ptGetSrcDestCoords
  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs, gr_meshComm
  use Logfile_interface, ONLY : Logfile_open, Logfile_close
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "gr_ptMapToMesh.h"

  integer,dimension(MAXBLOCKS), intent(IN) :: particlesPerBlk, blockList
  integer,intent(IN) :: blockCount
  integer,dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
  integer,dimension(MDIM), intent(IN) :: blkSize, guard
  integer, intent(OUT) :: BufferSize

  real, dimension(LOW:HIGH,MDIM) :: bndBlk
  real, dimension(MDIM) :: cellSpacing
  integer,dimension(LOW:HIGH,MDIM) :: guardCoords,srcCoords,destCoords
  integer,dimension(MDIM) :: guardCellID,srcCornerID,srcStride
  integer,dimension(BLKID:REFLEVELDIF,ABSMAXNEGH):: negh
  integer,dimension(MDIM,ABSMAXNEGH) :: neghCornerID
  integer :: blkNo, blockID
  integer :: regionIter, allCenters
  integer :: i, j, k, ib, jb, kb, ie, je, ke
  integer :: numNegh, ierr
  integer :: cellsToSend, eachDim, n
  integer :: maxRecvBufferSize, localSendSize, maxSendBufferSize
#ifdef DEBUG_GRIDMAPPARTICLES
  integer :: logUnit
  integer, save :: commIteration = 0
  logical, parameter :: logUnitLocal = .true.
  commIteration = commIteration + 1
#endif

  !MPI_IN_PLACE is an MPI-2 construct so we use 2 arrays, e.g.
  !gr_ptRecvTotalTmp and gr_ptRecvTotal.
  gr_ptRecvSpecifier = 0
  gr_ptRecvSpecifierTmp = 0
  gr_ptRecvTotalTmp = 0
  gr_ptRecvTotal = 0


  do blkNo = 1, blockCount
     blockID = blockList(blkNo)

     gr_ptDomain(blkNo) % blockID = NONEXISTENT
     if(particlesPerBlk(blockID)>0) then

        gr_ptDomain(blkNo) % blockID = blockID
        call Grid_getBlkBoundBox(blockID, bndBlk)
        gr_ptDomain(blkNo) % bndBlk(LOW:HIGH,1:NDIM) = bndBlk(LOW:HIGH,1:NDIM)
        call Grid_getDeltas(blockID, cellSpacing)
        gr_ptDomain(blkNo) % cellSpacing(1:NDIM) = cellSpacing(1:NDIM)

        call Grid_getBlkCornerID(blockID, srcCornerID, srcStride)

        regionIter = 0 
        allCenters = 2**NDIM
        kb=blkLimitsGC(LOW,KAXIS); ke=kb+K3D*(guard(KAXIS)-1)
        do k=LEFT_EDGE,LEFT_EDGE+K3D*(RIGHT_EDGE-1)
           jb=blkLimitsGC(LOW,JAXIS); je=jb+K2D*(guard(JAXIS)-1)
           do j=LEFT_EDGE,LEFT_EDGE+K2D*(RIGHT_EDGE-1)
              ib=blkLimitsGC(LOW,IAXIS); ie=ib+guard(IAXIS)-1
              do i=LEFT_EDGE,RIGHT_EDGE
                 if((i*j*k)/=allCenters) then !! make sure it is gCell region

                    regionIter = regionIter + 1
                    guardCellID(IAXIS)=i; guardCellID(JAXIS)=j; guardCellID(KAXIS)=k

                    !Find the neighboring blocks for this guard cell region.
                    call gr_ptFindNegh(blockID,guardCellID,negh,neghCornerID,numNegh)

                    gr_ptDomain(blkNo) % haloRegion(regionIter) % numNegh = numNegh

                    if(numNegh > 0) then

                       !Coordinates of the guard cell region of interest.
                       guardCoords(LOW,IAXIS)=ib; guardCoords(HIGH,IAXIS)=ie
                       guardCoords(LOW,JAXIS)=jb; guardCoords(HIGH,JAXIS)=je
                       guardCoords(LOW,KAXIS)=kb; guardCoords(HIGH,KAXIS)=ke


                       do n = 1,numNegh

                          call gr_ptGetSrcDestCoords(blkSize,guard,guardCellID,srcCornerID,srcStride, &
                               neghCornerID(:,n),guardCoords,negh(:,n),srcCoords,destCoords)


                          !If the neighbor belongs on another processor make space for its data.
                          if (negh(BLKPROC,n) /= gr_meshMe) then

                             cellsToSend = 1
                             do eachDim = 1, NDIM
                                cellsToSend = cellsToSend * (destCoords(HIGH,eachDim) - destCoords(LOW,eachDim) + 1)
                             end do

                             !Count the number of messages we are sending to each process.
                             gr_ptRecvSpecifierTmp(negh(BLKPROC,n)) = & 
                                  gr_ptRecvSpecifierTmp(negh(BLKPROC,n)) + 1
                             !Count the number of doubles we are sending to each process.
                             gr_ptRecvTotalTmp(negh(BLKPROC,n)) = &
                                  gr_ptRecvTotalTmp(negh(BLKPROC,n)) + SIZE_HEADER + cellsToSend

                          end if


                          gr_ptDomain(blkNo) % haloRegion(regionIter) % neighbor(n) % negh(:) = negh(:,n)
                          gr_ptDomain(blkNo) % haloRegion(regionIter) % neighbor(n) % cornerID(1:NDIM) = &
                               neghCornerID(1:NDIM,n)
                          gr_ptDomain(blkNo) % haloRegion(regionIter) % neighbor(n) % srcCoords(LOW:HIGH,1:NDIM) = &
                               srcCoords(LOW:HIGH,1:NDIM)
                          gr_ptDomain(blkNo) % haloRegion(regionIter) % neighbor(n) % destCoords(LOW:HIGH,1:NDIM) = &
                               destCoords(LOW:HIGH,1:NDIM)


                       end do !End of loop over neighbors

                    end if !End of if more than zero neighbors

                 end if !End of if all centers
                 ib=ie+1;ie=ib+mod(i,CENTER)*blkSize(IAXIS)+mod(i+1,CENTER)*guard(IAXIS)-1
              end do !End of do i
              jb=je+1;je=jb+mod(j,CENTER)*blkSize(JAXIS)+mod(j+1,CENTER)*guard(JAXIS)-1
           end do !End of do j
           kb=ke+1;ke=kb+mod(k,CENTER)*blkSize(KAXIS)+mod(k+1,CENTER)*guard(KAXIS)-1
        end do !End of do k

     end if !End of if a particle on this block
  end do !End of loop over blocks


  !Anything in tmp arrays are local calculations and describe what
  !we are sending to each process.  The subsequent reduced arrays contain
  !global information that is used to determine messages to be received.
  gr_ptNumMessagesToSend = sum(gr_ptRecvSpecifierTmp)
  localSendSize = sum(gr_ptRecvTotalTmp)

  call MPI_AllReduce(localSendSize, maxSendBufferSize, 1, MPI_INTEGER, &
       MPI_MAX, gr_meshComm, ierr)
  call MPI_Allreduce(gr_ptRecvSpecifierTmp, gr_ptRecvSpecifier, gr_meshNumProcs, &
       FLASH_INTEGER, MPI_SUM, gr_meshComm, ierr)
  call MPI_Allreduce(gr_ptRecvTotalTmp, gr_ptRecvTotal, gr_meshNumProcs, &
       FLASH_INTEGER, MPI_SUM, gr_meshComm, ierr)

  maxRecvBufferSize = maxval(gr_ptRecvTotal)
  bufferSize = max(maxRecvBufferSize, maxSendBufferSize)

#ifdef DEBUG_GRIDMAPPARTICLES
  if (gr_meshNumProcs <= 8) then
     call Logfile_open(logUnit,logUnitLocal)
     write(logUnit,'(a, i10)') "Communication iteration", commIteration
     do i = 0, gr_meshNumProcs-1
        write(logUnit,'(a,i10,a,i10,a,i10)') "Send", gr_ptRecvSpecifierTmp(i), &
             " messages of combined size", gr_ptRecvTotalTmp(i), &
             " to process", i
     end do
     write(logUnit,'(a,i10,a,i10)') "Send", sum(gr_ptRecvSpecifierTmp), &
          " messages of combined size", sum(gr_ptRecvTotalTmp)
     write(logUnit,'(a,i10,a,i10)') "Recv", gr_ptRecvSpecifier(gr_meshMe), &
          " messages of combined size", gr_ptRecvTotal(gr_meshMe)
     write(logUnit,'(a,i10)') "Global max send size", maxSendBufferSize
     write(logUnit,'(a,i10)') "Global max recv size", maxRecvBufferSize
     write(logUnit,'(a,i10)') "Global max send/recv size", bufferSize
     write(logUnit,'(a)') ""
     call Logfile_close(logUnitLocal)
  end if
#endif

end subroutine gr_ptStoreOffBlockCells
