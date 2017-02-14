!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/MoveSieve/gr_ptStoreOffBlockCells
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

!#ifndef DEBUG_GRIDMAPPARTICLES
!#define DEBUG_GRIDMAPPARTICLES
!#endif

subroutine gr_ptStoreOffBlockCells(particlesPerBlk, blockList, blockCount, blkLimitsGC, blkSize, guard, BufferSize)

  use gr_ptMapData, ONLY : gr_ptDomain
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkCornerID, Grid_getBlkBoundBox, Grid_getDeltas
  use gr_ptInterface, ONLY : gr_ptFindNegh, gr_ptGetSrcDestCoords
  use Grid_data, ONLY : gr_meshMe, gr_meshComm
  use gr_interfaceTypeDecl
  use gr_interface, ONLY : gr_checkGridState, gr_findAllNeghID


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
  integer :: totalHeader, totalData, blkNo, blockID
  integer :: regionIter, allCenters
  integer :: i, j, k, ib, jb, kb, ie, je, ke
  integer :: numNegh, localDerivedSize, ierr
  integer :: cellsToSend, eachDim, length, n
  type (AllBlockRegions_t) :: surrBlksSummary
  logical, dimension(2**(NDIM-1)) :: matchArray
  integer,dimension(BLKID:REFLEVELDIF) :: neghA, neghB
  integer :: a, b


#ifdef DEBUG_GRIDMAPPARTICLES
  !Useful in debug mode to run through all blocks to see if we can find valid data.
  do blkNo = 1, blockCount
     blockID = blockList(blkNo)
     call gr_findAllNeghID(blockID,surrBlksSummary)
  end do
  call MPI_Barrier(gr_meshComm, ierr)
#endif !DEBUG_GRIDMAPPARTICLES


  totalHeader = 0
  totalData = 0
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


        call gr_findAllNeghID(blockID,surrBlksSummary)

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



                    !The following is temporary code used to validate gr_ptFindNegh
                    !against gr_findAllNeghID.  First check the number of neighbors 
                    !that they return is the same.
                    !----------------------------- START ----------------------------
                    if (numNegh /= surrBlksSummary % regionInfo(i,j,k) % numNegh) then
                       call Driver_abortFlash("Number of neighbors mismatch")
                    end if

                    if (numNegh > 0) then
                       !Now we need to check that the actual neighbors are the same.

                       matchArray(1:numNegh) = .false.

                       A_gr_ptFindNegh: do a = 1, numNegh
                          neghA(BLKNO:TYPENO) = negh(BLKNO:TYPENO,a)

                          B_gr_findAllNeghID: do b = 1, numNegh 
                             neghB(BLKNO:TYPENO) = surrBlksSummary % regionInfo(i,j,k) % &
                                  details(BLKNO:TYPENO,b)

                             !DEV: Type is known to be different.  I need to convert it
                             !from format in surrBlksSummary to a type which
                             !mapping routines understand.
                             if (all(neghA(BLKNO:PROCNO) == neghB(BLKNO:PROCNO))) then
                                matchArray(a) = .true.
                                exit  !exits loop B, but stays in loop A.
                             end if

                          end do B_gr_findAllNeghID
                       end do A_gr_ptFindNegh

                       if (any(matchArray(1:numNegh) .eqv. .false.)) then
                          call Driver_abortFlash("Actual neighbor mismatch!")
                       end if

                    end if
                    !----------------------------- END   ----------------------------



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

                             !We have a single header for each neighbor belonging on another block.
                             totalHeader = totalHeader + SIZE_HEADER

                             cellsToSend = 1
                             do eachDim = 1, NDIM
                                cellsToSend = cellsToSend * (destCoords(HIGH,eachDim) - destCoords(LOW,eachDim) + 1)
                             end do

                             totalData = totalData + cellsToSend
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


  localDerivedSize = totalHeader + totalData
  call MPI_ALLREDUCE (localDerivedSize, BufferSize, 1, MPI_INTEGER, MPI_MAX, gr_meshComm, ierr)

  !print *, "Processor", gr_meshMe, "will set its send/receive buffer to size", BufferSize


#ifdef DEBUG_GRIDMAPPARTICLES
  !The following is temporary code used to check the grid state.
  !Moved to the bottom to have the aborts come AFTER diagnosis
  !----------------------------- START ----------------------------
  call gr_checkGridState()
  print *, "done with gr_checkGridState"
  !----------------------------- END   ----------------------------
#endif

end subroutine gr_ptStoreOffBlockCells
