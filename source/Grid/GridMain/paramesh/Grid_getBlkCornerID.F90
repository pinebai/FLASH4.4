!!****if* source/Grid/GridMain/paramesh/Grid_getBlkCornerID
!!
!! NAME
!!  Grid_getBlkCornerID
!!
!! SYNOPSIS
!!
!!  call Grid_getBlkCornerID(integer(IN)  :: blockId,
!!                           integer(OUT) :: cornerID(MDIM),
!!                           integer(OUT) :: stride(MDIM),
!!                  optional,integer(OUT) :: cornerIDHigh(MDIM))
!!  
!! DESCRIPTION 
!! 
!!  Returns the global integer indices of the left most interior zone
!!  of the block and the stride of indices along each dimension.
!!  Together the cornerID and the stride make a unique identifier for
!!  each block on the grid.
!!
!!  The block's cornerID is determined by calculating the global
!!  integer index of each cell on the global domain as if the grid
!!  were fully refined.  (The uniform grid is like a fully refined
!!  grid when calculating the cornerID.)  Another way to put it is
!!  that the cornerID is just the left most interior cell global index
!!  of a block in a dimension, as if the grid were fully refined.
!! 
!!  stride is defined as = 2**(lrefine_max - lrefine(blockId)) and
!!  indicates the spacing factor between one cell and one directly to
!!  its right when calculating the cornerID.
!!
!!  This means that if a block is at the highest refinement level,
!!  then it will always have a stride = 1 and the cornerID of the
!!  block to the current block's right will always be cornerID plus
!!  the index size of the block(NXB with fixed block sizes).  (The
!!  uniform grid always has a stride = 1.)
!!
!!  CornerID counting starts at 1 and only the interior cells (no
!!  guardcells) are used to calculate the cornerID.
!! 
!! 
!! ARGUMENTS 
!!
!!  blockId :: the local blockID
!!  cornerID :: global integer indices of start of the interior zone
!!              of the block
!!     
!!  stride  :: spacing factor between indices. In UG, stride is always = 1.
!!             For PARAMESH, stride may be more than 1, depending
!!             on how far down you are in the tree.
!!
!!  cornerIDHigh :: global integer indices of the last interior zone
!!              of the block
!!
!!  inRegion :: if present and true, cornerID is computed relative to region
!!              specified in Grid scope variable gr_region
!!
!! EXAMPLE
!!
!! Example 1:
!!  In a 1 dimensional UG case with 2 blocks and nxb=8
!!  The cornerID for block 1 = 1 and the cornerID for block 2 = 9 
!!
!!
!! Example 2:
!!  In a 1 dimensional PARAMESH case with lrefine_max = 4, nxb=8, the
!!  global indices blkLimits from 1:64.  If the grid is fully refined then
!!  there are 8 blocks in the x direction and their cornerIDs are
!!  1,9,17,25,33,41,49,57 respectively. stride=1
!!  
!!  If the entire grid is at a refinement level = 3 then there are 4
!!  blocks in the x direction and their cornerIDs are 1,17,33,49
!!  respectively.  stride = 2
!!
!!  If the entire grid is at a refinement level = 2 then there are 2
!!  blocks in the x direction and their cornerIDs are 1,33
!!  respectively.  stride = 4 (meaning it takes stride*nxb to get to
!!  the next block's cornerID)
!!
!!  And so on.
!!
!!  Multiple blocks can have the same cornerID as in the above
!!  example, but they can not have the same cornerID AND stride.
!! 
!!***

subroutine Grid_getBlkCornerID(blockId, cornerID, stride,cornerIDHigh, inRegion)
  use Grid_data, ONLY : gr_oneBlock,gr_globalDomain, gr_meshMe, &
       gr_nBlockX, gr_nBlockY, gr_nBlockZ,gr_region, gr_delta
  use tree, ONLY : lrefine, lrefine_max,bnd_box
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_open, Logfile_close
  implicit none

#include "constants.h"
#include "Flash.h"
  integer,intent(IN) :: blockId
  integer,dimension(MDIM), intent(OUT) :: cornerID, stride
  integer,dimension(MDIM),optional,intent(OUT) :: cornerIDHigh
  logical,optional,intent(IN) :: inRegion
  integer,dimension(MDIM) :: cID
  real, dimension(MDIM) :: delta, halfDelta
  real :: factor
  integer :: logUnit
  logical, parameter :: logUnitLocal = .true.
  real,dimension(NDIM) :: half_delta
  integer:: offset,i
  logical :: doRegionCalc

  cornerID = gr_oneBlock(blockId)%cornerID
  stride(:) = 2**(lrefine_max - lrefine(blockId))

  if (present(inRegion)) then
     doRegionCalc = inRegion
  else
     doRegionCalc = .false.
  end if

  if (doRegionCalc) then
     half_delta = gr_delta(1:NDIM,lrefine_max)/2.0
     do i = 1,NDIM
        if(gr_region(LOW,i)>gr_globalDomain(LOW,i)) then
           offset= (gr_region(LOW,i)-gr_globalDomain(LOW,i)+half_Delta(i))/&
                                  gr_delta(i,lrefine_max)
           cornerID(i)=cornerID(i)-offset
        end if
     end do
     if (present(cornerIDHigh)) then
        call Driver_abortFlash("Not clear what Anshu wants to do in this case")
     end if
  else

  !! This section of the API is to verify the calculation.

     factor = 2**(lrefine_max-1)
     
     delta(IAXIS) = (gr_globalDomain(HIGH,IAXIS) - gr_globalDomain(LOW,IAXIS)) & 
          / (1.0 * NXB * gr_nBlockX * factor)
#if NDIM >= 2
     delta(JAXIS) = (gr_globalDomain(HIGH,JAXIS) - gr_globalDomain(LOW,JAXIS)) & 
          / (1.0 * NYB * gr_nBlockY * factor)
#endif
#if NDIM == 3
     delta(KAXIS) = (gr_globalDomain(HIGH,KAXIS) - gr_globalDomain(LOW,KAXIS)) & 
          / (1.0 * NZB * gr_nBlockZ * factor)
#endif
     
     
     halfDelta(1:NDIM) = delta(1:NDIM)/2.0
     cID(1:NDIM) = ((bnd_box(LOW,1:NDIM,blockID) - gr_globalDomain(LOW,1:NDIM) + halfDelta(1:NDIM)) & 
          / delta(1:NDIM)) + 1
     
     if (any(cornerID(1:NDIM) /= cID(1:NDIM))) then
        print *, "Processor:", gr_meshMe, "block:", blockID, &
             "refinement:", lrefine(blockID), &
             "Original:", cornerID(1:NDIM), "verified:", cID(1:NDIM), &
             "bnd_box:", bnd_box(LOW:HIGH,1:NDIM,blockID)
        
        call Logfile_open(logUnit,logUnitLocal)
        write(logUnit,*) "**** Corner ID mismatch ****"
        write(logUnit,*) "block:", blockID, &
             "refinement:", lrefine(blockID), &
             "Original:", cornerID(1:NDIM), "verified:", cID(1:NDIM), &
             "bnd_box:", bnd_box(LOW:HIGH,1:NDIM,blockID)
        call Logfile_close(logUnitLocal)
        call Driver_abortFlash("corner ID calculations inconsistent")    
     end if
     if(present(cornerIDHigh)) then
        cornerIDHigh(IAXIS)=cornerID(IAXIS)+stride(IAXIS)*NXB-1
        cornerIDHigh(JAXIS)=cornerID(JAXIS)+stride(JAXIS)*NYB-1
        cornerIDHigh(KAXIS)=cornerID(KAXIS)+stride(KAXIS)*NZB-1
     end if
  end if
end subroutine Grid_getBlkCornerID

