!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftGenMap
!!
!! NAME
!!
!!  gr_pfftGenMap
!!
!! SYNOPSIS
!!  
!!  gr_pfftGenMap()
!!
!! DESCRIPTION
!!
!! This routine does not send the actual data, but uses communication
!! to generate the information about the communication pattern to be expected
!! when actual data movement happens.
!!
!!***

subroutine gr_pfftGenMap()
#include "Flash.h"
#include "constants.h"
#include "Pfft.h"

  use Logfile_interface, ONLY : Logfile_open, Logfile_close, Logfile_stamp
  use gr_pfftData, ONLY : pfft_procGrid, pfft_ndim, pfft_myPE, &
       pfft_me, pfft_commWithTopology, pfft_numProcs, pfft_inLen, pfft_regionBndBox,&
       pfft_inRegion
  use gr_pfftReconfigData, ONLY : pfft_pencilSize, pfft_procLookup
  use gr_pfftinterface, ONLY : gr_pfftGetDestPfftCoords, &
       gr_pfftCreateSendNode, gr_pfftCommunicateNodeMetaData, gr_pfftGridPointTable
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkCornerID,&
       Grid_getBlkIndexLimits
  use Driver_interface, ONLY : Driver_abortFlash
#ifdef FLASH_GRID_PARAMESH
  use Grid_data, ONLY : gr_oneRefLev
#endif

  implicit none
  integer,dimension(IAXIS:KAXIS) :: stride,cornerID, blkSize
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(MAXBLOCKS) :: blkList
  integer, dimension(1:MDIM) :: numProcs, firstPfftProc, guard, &
       pfftProcCoords, pfftProcOffset
  !Prefix g means global, prefix l means local.
  integer, dimension(1:MDIM) :: gBlockStart, gBlockEnd, &
       gBlockFragStart, gBlockFragEnd, gPencilStart, &
       gPencilEnd, lBlockFragStart, lBlockFragEnd, &
       lPencilFragStart, lPencilFragEnd
  integer :: axis, blkCount, blockID, fragmentSize, i, iProc, jProc, kProc, &
       overallBlkSize, sumFragmentSize, pencilSize, commBufMemEstimate, &
       ierr, pfftProc, dims, procPos, pencilStart, pencilEnd
  logical :: zeroSizePencil

  if(pfft_ndim == 1) then
     return
  end if

  !This is a real "back of the cigarette packet" calculation, but the 
  !average memory / process devoted to our persistent communication buffers is:
  commBufMemEstimate = product(pfft_inLen(1:pfft_ndim)) * 2 * 8 !8-byte reals.
  !(This is reasonably accurate if each process has the same number of FLASH 
  !grid points to send.  Note, the max memory / process for just the receive 
  !buffer is: product(pfft_inLen(1:pfft_ndim)) * 8.)
  call Logfile_stamp( commBufMemEstimate, &
    "[gr_pfftGenMap] Estimated comm. buffer allocation (bytes)")

  !Stores values in pfft_procLookup.
  call gr_pfftGridPointTable(pfft_inLen)

  !Initialise any variables, arrays, datatypes.
  call gr_pfftInitMapData()
  

  !Calculate the number of grid points this processors will own in pencil space.
  !----------------------------------------------------------------------------
  do axis = 1, pfft_ndim
     gPencilStart(axis) = &
          pfft_procLookup(axis) % procInfo(pfft_me(axis)) % globalStartGridPoint
     gPencilEnd(axis) = &
          pfft_procLookup(axis) % procInfo(pfft_me(axis)) % globalEndGridPoint
  end do
  
  if ( (any(gPencilStart(1:pfft_ndim) == NONEXISTENT)) .or. &
       (any(gPencilEnd(1:pfft_ndim) == NONEXISTENT)) ) then
     !A NONEXISTENT means this processor will own no points in pencil space.
     !However, we still need this processor because it has FLASH grid points 
     !which must be sent to a pencil space processor.
     pfft_pencilSize = 0
     print *, "(INFO) Processor:", pfft_myPE, "has no pencil grid points."
  else
     pfft_pencilSize = product(&
          (gPencilEnd(1:pfft_ndim) - gPencilStart(1:pfft_ndim) + 1))
  end if
  !----------------------------------------------------------------------------


#if defined(FLASH_GRID_UG)
  call Grid_getListOfBlocks(LEAF,blkList,blkCount)
  
#elif defined(FLASH_GRID_PARAMESH)
  if(pfft_inRegion) then
     if(gr_oneRefLev==NONEXISTENT) then
        call Grid_getListOfBlocks(INREGION,blkList,blkCount,region_bndBox=pfft_regionBndBox)        
     else
        call Grid_getListOfBlocks(INREGION,blkList,blkCount,gr_oneRefLev,pfft_regionBndBox)        
     end if
  else
     call Grid_getListOfBlocks(REFINEMENT,blkList,blkCount,gr_oneRefLev)
  end if
#endif


  blockLoop: do i = 1, blkCount

     blockID = blkList(i)
     call Grid_getBlkCornerID(blockID,cornerID,stride,inRegion=pfft_inRegion)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     blkSize(1:MDIM) = blkLimits(HIGH,1:MDIM) - blkLimits(LOW,1:MDIM) + 1
     guard(1:MDIM) = blkLimits(LOW,1:MDIM) - blkLimitsGC(LOW,1:MDIM)

     !Starting position of the block if this was assumed to be the 
     !UG at least refined level
     gBlockStart(1:MDIM) = 1  !For <3D sims
     gBlockStart(1:NDIM) = ceiling(real(cornerID(1:NDIM)) / real(stride(1:NDIM)))
     gBlockEnd(1:MDIM) = 1    !For <3D sims
     gBlockEnd(1:NDIM) = gBlockStart(1:NDIM) + blkSize(1:NDIM) - 1
     overallBlkSize = product(blkSize(1:MDIM))

     
     !We first need to determine how many fragments to create 
     !from this one FLASH block.  We use the numProcs array to 
     !store how many PFFT processors the block is distributed over.
     numProcs(1:MDIM) = 1
     call gr_pfftGetDestPfftCoords(gBlockStart, firstPfftProc)
     
     do axis = 1, pfft_ndim        
        procPos = firstPfftProc(axis)

        fragmentLoop: do 

           if (procPos >= pfft_procGrid(axis)) then 
              call Driver_abortFlash &
                   ("[gr_pfftGenMap]: Going to overrun data structure")
           end if

           gPencilEnd(axis) = & 
                pfft_procLookup(axis) % procInfo(procPos) % globalEndGridPoint

           !If the end position of the FLASH block is beyond the 
           !PFFT processor end position, it means that some of the block must
           !reside on the next PFFT processor.
           if (gPencilEnd(axis) < gBlockEnd(axis)) then          
              procPos = procPos + 1
              numProcs(axis) = numProcs(axis) + 1
           else
              exit
           end if

        end do fragmentLoop
     end do


     sumFragmentSize = 0
     !Now we loop over each processor and calculate the fragment size and
     !position.  We then use this information to construct a fragment object 
     !which is then appended to a list.

     !In <3D simulations we must ensure unused dimensions hold the value 1.
     gBlockFragStart = 1; gBlockFragEnd = 1
     lBlockFragStart = 1; lBlockFragEnd = 1
     lPencilFragStart = 1; lPencilFragEnd = 1


     kAxisLoop: do kProc = 0, numProcs(KAXIS)-1
        pfftProcCoords(KAXIS) = firstPfftProc(KAXIS) + kProc
        if (kProc == 0) then
           gBlockFragStart(KAXIS) = gBlockStart(KAXIS)
        else
           !There are no "holes" in the pencil grid, so we can calculate 
           !the global fragment start position by adding 1 to the previous 
           !end position.
           gBlockFragStart(KAXIS) = gBlockFragEnd(KAXIS) + 1
        end if


        jAxisLoop: do jProc = 0, numProcs(JAXIS)-1
           pfftProcCoords(JAXIS) = firstPfftProc(JAXIS) + jProc
           if (jProc == 0) then
              gBlockFragStart(JAXIS) = gBlockStart(JAXIS)
           else
              gBlockFragStart(JAXIS) = gBlockFragEnd(JAXIS) + 1
           end if


           iAxisLoop: do iProc = 0, numProcs(IAXIS)-1
              pfftProcCoords(IAXIS) = firstPfftProc(IAXIS) + iProc
              if (iProc == 0) then
                 gBlockFragStart(IAXIS) = gBlockStart(IAXIS)
              else
                 gBlockFragStart(IAXIS) = gBlockFragEnd(IAXIS) + 1
              end if


              do axis = 1, pfft_ndim                  
                   
                 gPencilStart(axis) = pfft_procLookup(axis) &
                      % procInfo(pfftProcCoords(axis)) % globalStartGridPoint
                 gPencilEnd(axis) = pfft_procLookup(axis) &
                      % procInfo(pfftProcCoords(axis)) % globalEndGridPoint
                 

                 !The fragment end position is either:
                 !  * The pencil end position or/ * The block end position.
                 if (gBlockEnd(axis) > gPencilEnd(axis)) then
                    gBlockFragEnd(axis) = gPencilEnd(axis)
                 else
                    gBlockFragEnd(axis) = gBlockEnd(axis)
                 end if


                 !These describe the position of the grid points relative to the
                 !actual FLASH block of interest.
                 lBlockFragStart(axis) = gBlockFragStart(axis) - &
                      gBlockStart(axis) + guard(axis) + 1  
                 lBlockFragEnd(axis) = gBlockFragEnd(axis) - &
                      gBlockStart(axis) + guard(axis) + 1


                 !These describe the position of the grid points relative to the
                 !actual pencil block of interest.
                 lPencilFragStart(axis) = gBlockFragStart(axis) - &
                      gPencilStart(axis) + 1
                 lPencilFragEnd(axis) = lPencilFragStart(axis) + &
                      (gBlockFragEnd(axis) - gBlockFragStart(axis))

              end do
              
              call MPI_Cart_rank(pfft_commWithTopology, pfftProcCoords(1), &
                   pfftProc, ierr)

              fragmentSize = product( & 
                   (gBlockFragEnd(1:NDIM) - gBlockFragStart(1:NDIM) + 1) )

              sumFragmentSize = sumFragmentSize + fragmentSize
              if (sumFragmentSize > overallBlkSize) then
                 call Driver_abortFlash ("[gr_pfftGenMap]:" // &
                      "Block fragment calculation error (1)")
              end if
              
#ifdef DEBUG_PFFT
              print *, "Processor:", pfft_myPE, "sending a fragment to:", &
                   pfftProc, "Start:", gBlockFragStart, &
                   ", end:", gBlockFragEnd
#endif DEBUG_PFFT

              call gr_pfftCreateSendNode(pfft_myPE, blockID,  &
                   lBlockFragStart, lBlockFragEnd, &
                   pfftProc, lPencilFragStart, lPencilFragEnd)

           end do iAxisLoop
        end do jAxisLoop
     end do kAxisLoop
     !-----------------------------------------------------------------------

     !Check that the sum of fragments sizes is equal to the total block size.
     if (sumFragmentSize /= overallBlkSize) then  
        call Driver_abortFlash("[gr_pfftGenMap]: " // &
             "Block fragment calculation error (2)")
     end if

  end do blockLoop

  call gr_pfftCommunicateNodeMetaData()

end subroutine gr_pfftGenMap
