!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/gr_ptGetSrcDestCoords
!!
!! NAME
!!  gr_ptGetSrcDestCoords
!!
!! SYNOPSIS
!!
!!  gr_ptGetSrcDestCoords(integer,dimension(MDIM), intent(IN)  :: blkSize
!!                        integer,dimension(MDIM), intent(IN)  :: guard
!!                        integer,dimension(MDIM), intent(IN)  :: guardCellID
!!                        integer,dimension(MDIM), intent(IN)  :: srcCornerID
!!                        integer,dimension(MDIM), intent(IN)  :: destCornerID
!!                        integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: guardCoords
!!                        integer,dimension(BLKID:REFLEVELDIF), intent(IN):: negh
!!                        integer,dimension(LOW:HIGH,MDIM), intent(OUT)  :: srcCoords
!!                        integer,dimension(LOW:HIGH,MDIM), intent(OUT)  :: destCoords
!!
!! DESCRIPTION
!!
!! This procedure determines two key pieces of information:
!!   1.  The coordinate range of the guard cells in the source block (srcCoords). 
!!   2.  The coordinate range in the destination block in which we will copy the 
!!       source block's guard cells (destCoords).
!!
!! The corner IDs of the source block (srcCornerID) and destination block (destCornerID) 
!! are crucial, as it allows us to find srcCoords and destCoords, even when the 
!! destination block exists on another processor and its block ID is unknown.
!!
!! The procedure takes into account the resolution difference of abutting blocks. 
!!
!! ARGUMENTS
!!               blkSize:  An array containing the size of the source block.
!!               guard:  An array containing the number of guard cells for the source block.
!!               guardCellID:  An array containing the relative coordinates of the guard
!!                             cell region with respect to the source block.
!!               srcCornerID:  An array containing the source block corner ID.
!!               destCornerID:  An array containing the destination block corner ID.
!!               guardCoords:  An array containing the coordinates of the guard cell region
!!                             with respect to the source block.
!!               negh:  An array containing information about the destination block.
!!               srcCoords:  An array containing the coordinates of the guard cell region 
!!                           that need to be copied to the destination block.
!!               destCoords:  An array containing the coordinates in the destination block
!!                            in which to copy the source block guard cells.
!!
!! NOTES
!!
!! guardCoords and srcCoords can differ if the neighboring block is more refined.  
!! This only happens when the source block's guard cell region has several 
!! neighboring destination blocks (i.e. when there are 2 or 4 neighboring blocks).
!!
!!***

subroutine gr_ptGetSrcDestCoords(blkSize, guard, guardCellID, srcCornerID, srcStride, destCornerID, &
     guardCoords, negh, srcCoords, destCoords)     

  use Driver_interface, ONLY : Driver_abortFlash
  use gr_ptMapData, ONLY : gr_ptSmearLen
  implicit none

#include "constants.h"
#include "Flash.h"
#include "gr_ptMapToMesh.h"

  integer,dimension(MDIM), intent(IN) :: blkSize, guard, guardCellID, srcCornerID, srcStride, destCornerID
  integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: guardCoords
  integer,dimension(BLKID:REFLEVELDIF), intent(IN):: negh
  integer,dimension(LOW:HIGH,MDIM), intent(OUT)  :: srcCoords, destCoords

  integer,dimension(MDIM) :: guardCellZeroBasedID
  integer :: eachAxis, value, idSpan, oddValue, nBlksFromDomCorner, numbSrcCells, numbDestCells
  integer :: smearLength

  !Perform neccessary initialisations. We may later adjust these values.
  !It is also useful for simulations in which NDIM < MDIM.
  srcCoords(LOW:HIGH,1:MDIM) = guardCoords(LOW:HIGH,1:MDIM)
  destCoords(LOW:HIGH,1:MDIM) = guardCoords(LOW:HIGH,1:MDIM)

  smearLength = gr_ptSmearLen
  !Reduce the quantity of cells that need to be exchanged by considering
  !the smear length.  Here, smearLength is zero for NGP and one for CIC and TSC.
  !i.e. exhange only cells that can possibly receive mapping.
  !(OPTIMISATION, WHICH HAS NOT YET BEEN INCORPORATED!)
  !numbSrcCells = max(1, smearLength)  !We shouldn't be here if smearLength==0, but just in case.
  !do eachAxis = 1, NDIM
  !   if(guardCellID(eachAxis) == LEFT_EDGE) then
  !      srcCoords(LOW,eachAxis) = srcCoords(HIGH,eachAxis) - (numbSrcCells - 1)
  !   else if(guardCellID(eachAxis) == RIGHT_EDGE) then
  !      srcCoords(HIGH,eachAxis) = srcCoords(LOW,eachAxis) + (numbSrcCells - 1)
  !   end if
  !end do


  if(negh(REFLEVELDIF)==0) then
     !Neighbor block is at the same refinement.

     !Gives a number between -1 and +1, which identifies the guard 
     !cell region for each dimension in a block.
     guardCellZeroBasedID(1:NDIM)=guardCellID(1:NDIM)-2  

     !e.g. blkSize = 8, guard = 4
     !---------------------------
     !
     !guardCoords for IAXIS: (ib=1, ie=4): on dest block, we want: (iib=9, iie=12)
     !guardCellZeroBasedID = -1 (block to left)
     !destCoords for IAXIS: 1-(-1*8)=9, 4-(-1*8)=12
     !
     !guardCoords for IAXIS: (ib=13, ie=16): on dest block, we want: (iib=5, iie=8)
     !guardCellZeroBasedID = +1 (block to right)
     !destCoords for IAXIS: 13-(1*8)=5, 16-(1*8)=8
     destCoords(LOW,1:NDIM) = srcCoords(LOW,1:NDIM) - (guardCellZeroBasedID(1:NDIM) * blkSize(1:NDIM))
     destCoords(HIGH,1:NDIM) = srcCoords(HIGH,1:NDIM) - (guardCellZeroBasedID(1:NDIM) * blkSize(1:NDIM))


  else if(negh(REFLEVELDIF)==-1) then
     !Neighbor block is less refined.

     do eachAxis = 1, NDIM

        idSpan = srcStride(eachAxis) * blkSize(eachAxis)
        nBlksFromDomCorner = (srcCornerID(eachAxis)-1) / idSpan

        !The number of blocks will be even or odd
        oddValue = mod(nBlksFromDomCorner,2)  !0 is even, 1 is odd.
        numbDestCells = max(1, guard(eachAxis)/2)   !numbDestCells = max(1, smearLength/2)   !eventually.

        if(guardCellID(eachAxis) == LEFT_EDGE) then

           destCoords(HIGH,eachAxis) = 1 + guard(eachAxis) + (blkSize(eachAxis)/2) + & 
                ((1-oddValue)*(blkSize(eachAxis)/2)) - 1
           destCoords(LOW,eachAxis) = destCoords(HIGH,eachAxis) - (numbDestCells - 1)


        else if(guardCellID(eachAxis) == RIGHT_EDGE) then

           destCoords(LOW,eachAxis) = 1 + guard(eachAxis) + ((1-oddValue) * (blkSize(eachAxis)/2))
           destCoords(HIGH,eachAxis) = destCoords(LOW,eachAxis) + (numbDestCells - 1)


        else if(guardCellID(eachAxis) == CENTER) then

           !Destination covers the entire size of the source block in this dimension.
           destCoords(LOW,eachAxis) = 1 + guard(eachAxis) + (oddValue * blkSize(eachAxis)/2)
           destCoords(HIGH,eachAxis) = destCoords(LOW,eachAxis) + (blkSize(eachAxis)/2) - 1

        else

           call Driver_abortFlash("[gr_ptGetSrcDestCoords]: Unrecognised guard cell region identifier.")

        end if
     end do


  else if(negh(REFLEVELDIF)==1) then
     !Neighbor block is more refined.

     do eachAxis = 1, NDIM

        numbDestCells = max(1, guard(eachAxis)*2)  !numbDestCells = max(1,smearLength*2)   !eventually.

        if(guardCellID(eachAxis) == LEFT_EDGE) then

           destCoords(HIGH,eachAxis) = 1 + guard(eachAxis) + blkSize(eachAxis) - 1
           destCoords(LOW,eachAxis) = destCoords(HIGH,eachAxis) - (numbDestCells - 1)


        else if(guardCellID(eachAxis) == RIGHT_EDGE) then

           destCoords(LOW,eachAxis) = 1 + guard(eachAxis)
           destCoords(HIGH,eachAxis) = destCoords(LOW,eachAxis) + (numbDestCells - 1)


        else if(guardCellID(eachAxis) == CENTER) then
           !Neighbor is level in this dimension.
           !Therefore the destination is the complete neighbor section.
           destCoords(LOW,eachAxis) = 1 + guard(eachAxis)
           destCoords(HIGH,eachAxis) = guard(eachAxis) + blkSize(eachAxis)


           !But which subsection in the source block do we need to copy?
           !Check the relative location of the destination block.
           idSpan = (srcStride(eachAxis)/2) * blkSize(eachAxis)
           nBlksFromDomCorner = (destCornerID(eachAxis)-1) / idSpan
           oddValue = mod(nBlksFromDomCorner,2)  !0 is even, 1 is odd.


           !Copy from a region half the size of source block in this dimension.
           srcCoords(LOW,eachAxis) = 1 + guard(eachAxis) + (oddValue * blkSize(eachAxis)/2)
           srcCoords(HIGH,eachAxis) = srcCoords(LOW,eachAxis) + (blkSize(eachAxis)/2) - 1

        else

           call Driver_abortFlash("[gr_ptGetSrcDestCoords]: Unrecognised guard cell region identifier.")

        end if
     end do

  else

     call Driver_abortFlash("[gr_ptGetSrcDestCoords]: Unexpected refinement.")

  end if


  return

end subroutine gr_ptGetSrcDestCoords
