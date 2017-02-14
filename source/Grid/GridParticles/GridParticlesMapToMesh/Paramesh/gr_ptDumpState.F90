!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/gr_ptDumpState
!!
!! NAME
!!  gr_ptDumpState
!!
!! SYNOPSIS
!!
!!  gr_ptDumpState(integer(IN) :: bufferSize
!!                real(IN), dimension(bufferSize) :: dataBuffer
!!                integer(IN) :: bufferContentSize)
!!
!! DESCRIPTION
!!
!! This subroutine is only invoked when the data buffer containing the
!! off-proc grid points has been passed around all processors.  This only
!! occurs when we have an error, as all data should have been extracted by now.
!! It does the following:
!!
!! A). Prints metadata about the destination of the unmatched grid points.
!! B). Prints metadata about the blocks actually on this processor.
!! C). Prints metadata about the neighboring blocks (found in gr_ptFindNegh()).
!!
!! All of this is printed to file in distinct (A,B,C) sections for each 
!! processor.  The cause of the unmatched grid points can then hopefully 
!! resolved by grep(ing) through the different sections of each data file.
!!
!!
!! ARGUMENTS
!!               bufferSize: Size of the communication buffer.
!!               dataBuffer: The actual communication buffer.
!!               bufferContentSize: The amount of actual data in the 
!!                                  the communication buffer.
!!
!! PARAMETERS
!! 
!!***

subroutine gr_ptDumpState(bufferSize, dataBuffer, bufferContentSize)

  use gr_ptMapData, ONLY : gr_ptSmearLen, gr_ptDomain, NUMBGUARDREGIONS
  use gr_ptData, ONLY : gr_ptBlkCount, gr_ptBlkList
  use gr_ptInterface, ONLY : gr_ptParseMetadata
  use Grid_data, ONLY : gr_meshMe
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkCornerID
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_open, Logfile_close
  use tree, ONLY : lrefine

#include "Flash.h"
#include "constants.h"
#include "gr_ptMapToMesh.h"

  implicit none
  integer,intent(IN) :: bufferSize
  real,dimension(bufferSize),intent(IN) :: dataBuffer
  integer,intent(IN) :: bufferContentSize

  integer, dimension(LOW:HIGH,MDIM) :: sectionCoords
  integer, dimension(MAXBLOCKS) :: listofBlocks
  integer, dimension(MDIM) :: cornerID, stride, negh, neghCornerID
  integer :: headerPtr, numbElements, blk, regionIter, n, numNegh
  integer :: blkCount, blockID, logUnit
  logical, parameter :: logUnitLocal = .true.

  !We leave 12 spaces for integers, so that we can capture 
  !any 4-byte integer value with space to spare.
  character(len=*), parameter  :: FMT1D = &
       "(1X, A, 1X, I12, 1X, I12, 1X, I12, 1X, I12)"
  character(len=*), parameter  :: FMT2D = &
       "(1X, A, 1X, I12, 1X, I12, 1X, I12, 1X, I12, 1X, I12)"
  character(len=*), parameter  :: FMT3D = &
       "(1X, A, 1X, I12, 1X, I12, 1X, I12, 1X, I12, 1X, I12, 1X, I12)"

#if NDIM==1
  character(len=*), parameter :: FMTString = FMT1D
#elif NDIM==2
  character(len=*), parameter :: FMTString = FMT2D
#elif NDIM==3
  character(len=*), parameter :: FMTString = FMT3D
#endif


  call Logfile_open(logUnit,logUnitLocal)

  if (gr_ptSmearLen <= 0) then

     write(logUnit,*) &
          "[gr_ptDumpState]: gr_ptSmearLen <= 0: Shouldn't be communicating!"

  else

     write(logUnit,*) "NOTE: level indicates refinement."
     write(logUnit,*) "It is relative for A & C, where: " // &
          "Neighbor is more refined (+1), same (0), less refined (-1)"
     write(logUnit,*) &
          "ID            block        proc         level      " // &
          "cornerX      cornerY      cornerZ"

     !Dump metadata describing unmatched grid points.
     !----------------------------------------------------------
     if (bufferContentSize > 1) then
        !Only scan though the metadata if there is metadata to scan.
        !gr_ptMoveMappedData sets sendCount (what we call bufferContentSize) 
        !to 1 if this process has no data to send.
        headerPtr = 1
        EachHeader: do

           !Extract the data description.
           call gr_ptParseMetadata(bufferSize, dataBuffer, headerPtr, negh, &
                neghCornerID, sectionCoords, numbElements)
           
           write(logUnit, FMTString) "A:", negh(BLKID), negh(BLKPROC), &
                negh(REFLEVELDIF), neghCornerID(1:NDIM)
           
           !Skip over the data and move onto the next header description.
           headerPtr = headerPtr + SIZE_HEADER + numbElements
           
           !All data has been checked, so break from the extraction loop.
           if((headerPtr-1) == bufferContentSize) then
              exit 
           else if((headerPtr-1) > bufferContentSize) then
              call Driver_abortFlash("[gr_ptDumpState]: Metadata loop will not exit cleanly!") 
           end if
           
        end do EachHeader
     end if
     !----------------------------------------------------------



     !Dump the block ID and corner ID of every block on this processor.
     !----------------------------------------------------------
     call Grid_getListOfBlocks(ALL_BLKS, listofBlocks, blkCount)
     do blk = 1, blkCount
        blockID = listofBlocks(blk)

        call Grid_getBlkCornerID(blockID, cornerID, stride)

        write(logUnit, FMTString) "B:", blockID, gr_meshMe, &
             lrefine(blockID), cornerID(1:NDIM)

     end do
     !----------------------------------------------------------



     !Dump the contents of the neighbor information data structure.
     !This contains all neighbors for each LEAF block that has particles.  
     !Mistakes may have crept in here because finding all neighbors is non-trivial.
     !----------------------------------------------------------
     do blk = 1, gr_ptBlkCount   !Count of number of LEAF blocks.

        !Neighbors were calclulated if the block contained particles.
        if (gr_ptDomain(blk) % blockID /= NONEXISTENT) then

           do regionIter = 1, NUMBGUARDREGIONS
              numNegh = gr_ptDomain(blk) % haloRegion(regionIter) % numNegh

              do n = 1, numNegh              
                 negh(:) = gr_ptDomain(blk) % &
                      haloRegion(regionIter) % neighbor(n) % negh(:)

                 neghCornerID(1:NDIM) = gr_ptDomain(blk) % &
                      haloRegion(regionIter) % neighbor(n) % cornerID(1:NDIM)

                 write(logUnit, FMTString) "C:", negh(BLKID), negh(BLKPROC), &
                      negh(REFLEVELDIF), neghCornerID(1:NDIM)
              end do

           end do  !End of guard cell region loop.
        end if  !End of if neighbors calculated for this block.
     end do  !End of LEAF block loop.
     !----------------------------------------------------------

  end if
  
  call Logfile_close(logUnitLocal)
  
end subroutine gr_ptDumpState
