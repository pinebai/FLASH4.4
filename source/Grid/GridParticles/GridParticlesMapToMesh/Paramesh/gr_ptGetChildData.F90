!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/gr_ptGetChildData
!!
!! NAME
!!  gr_ptGetChildData
!!
!! SYNOPSIS
!!
!!  gr_ptGetChildData(integer(IN)  :: guardCellID(MDIM), &
!!                    integer(IN)  :: blkSize(MDIM), &
!!                    integer(IN)  :: parentCornerID(MDIM), &
!!                    integer(IN)  :: srcStride(MDIM), &
!!                    integer(IN)  :: neghBlkID, &
!!                    integer(IN)  :: neghProcID, &
!!                    integer(OUT) :: negh(BLKID:REFLEVELDIF,ABSMAXNEGH), &
!!                    integer(OUT) :: neghCornerID(MDIM,ABSMAXNEGH), &
!!                    integer(OUT) :: numNegh
!!
!!
!! DESCRIPTION
!!
!!  Helper routine for gr_ptFindNegh.  Obtains the block and process ID of 
!!  the neighbor LEAF blocks at relative position specified by guardCellID.  These 
!!  neighbor blocks are at a finer refinement level than our source block.
!!  Our 'real' neighbor is a PARENT which is identified by neghBlkID, neghProcID.
!!
!!
!! ARGUMENTS
!!
!!  guardCellID:    The guard cell region of interest relative to source block
!!  blkSize:        The size of our source block.
!!  parentCornerID: The corner ID of our 'real' neighbor (has block Type=PARENT)
!!  srcStride:      The stride at source block level.
!!  neghBlkID:      The block ID of our 'real' neighbor (has block Type=PARENT)
!!  neghProcID:     The proc ID of our 'real' neighbor (has block Type=PARENT)
!!  negh:           Data structure for storing blkID/procID of all LEAF neighbors
!!  neghCornerID:   Data structure for storing corner IDs of all LEAF neighbors
!!  numNegh:        Number of LEAF neighbors.
!!
!!
!! NOTES
!!
!!***

!Useful to turn on additional assertations in this file by default.
#ifndef DEBUG_GRIDMAPPARTICLES
#define DEBUG_GRIDMAPPARTICLES
#endif

#include "constants.h"
#include "Flash.h"
#include "gr_ptMapToMesh.h"

subroutine gr_ptGetChildData(guardCellID, blkSize, parentCornerID, srcStride, &
     neghBlkID, neghProcID, negh, neghCornerID, numNegh)

  use gr_ptInterface, ONLY : gr_ptSearchBlk
  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs
  use gr_interface, ONLY : gr_getBlkHandle
  use Driver_interface, ONLY : Driver_abortFlash
  use tree, ONLY : child

  implicit none
  interface
     subroutine gr_ptFindValidChildren(maxDims, maxChildren, blockChild, &
          guardCellID, dims, actualChildren, validChild)
       implicit none
       integer, intent(IN) :: maxDims, maxChildren
       integer, dimension(maxDims,maxChildren), intent(IN) :: blockChild
       integer, dimension(maxDims), intent(IN) :: guardCellID
       integer, intent(IN) :: dims, actualChildren
       logical, dimension(actualChildren), intent(OUT) :: validChild
     end subroutine gr_ptFindValidChildren
  end interface


  integer, dimension(MDIM), intent(IN) :: guardCellID, blkSize, &
       parentCornerID, srcStride
  integer, intent(IN) :: neghBlkID, neghProcID
  integer, dimension(BLKNO:TYPENO,ABSMAXNEGH), intent(OUT) :: negh
  integer, dimension(MDIM,ABSMAXNEGH), intent(OUT) :: neghCornerID
  integer, intent(OUT) :: numNegh

  integer :: validChildren, eachChild, childCounter, blkHandle

  integer, parameter :: TOTAL_CHILDREN = 2**MDIM
  integer, parameter :: USABLE_CHILDREN = 2**NDIM

  !This is how the children are laid out in a parent block.
  integer, parameter, dimension(MDIM,TOTAL_CHILDREN) :: blockChild = &
       reshape (source = (/&
       LOW,  LOW,  LOW , &
       HIGH, LOW,  LOW , &
       LOW,  HIGH, LOW , &
       HIGH, HIGH, LOW , &
       LOW,  LOW,  HIGH, &
       HIGH, LOW,  HIGH, &
       LOW,  HIGH, HIGH, &
       HIGH, HIGH, HIGH  &
       /), shape = (/MDIM,TOTAL_CHILDREN/))
  logical, dimension(USABLE_CHILDREN) :: validChild
  integer, dimension(MDIM) :: switchCornerID

  !We dimension "negh" argument according to BLKNO:TYPENO and not 
  !BLKID:REFLEVELDIF.  This is so that the function prototype can be
  !placed in gr_ptInterface without needing to include gr_ptMapToMesh.h.  
  !We set REFLEVELDIF equal to TYPENO and BLKID to BLKNO in 
  !gr_ptMapToMesh.h, so the following is just a double-check.
  if ( (TYPENO /= REFLEVELDIF) .or. (BLKNO /= BLKID) ) then
     call Driver_abortFlash("Preprocessor definitions must be the same")
  end if

  negh = -200   !Random initial value for spotting errors.
  neghCornerID = 1  !For when NDIM < MDIM.
  numNegh = NONEXISTENT
  validChildren = 0
  switchCornerID = 0

  !Update a logical array named validChild in which .true. elements
  !indicate neighboring child blocks touching the guard cell region.
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, USABLE_CHILDREN, validChild)


  !Check that the number of child blocks we calculate is correct.
  numNegh = 2 ** (count(guardCellID(1:NDIM) == CENTER))
  validChildren = count(validChild(1:USABLE_CHILDREN) .eqv. .true.)
  if (validChildren /= numNegh) then
     print *, "Predicted:", numNegh, ", Calculated:", validChildren
     call Driver_abortFlash("[gr_ptGetChildData]: Number of children mismatch")
  end if


  !-------------------------------------------------------------------
  !Corner IDs
  !-------------------------------------------------------------------
  childCounter = 0
  do eachChild = 1, USABLE_CHILDREN
     if (validChild(eachChild) .eqv. .true.) then      
        
        childCounter = childCounter + 1
        
        !LOW = 1, HIGH = 2.  Therefore, we keep parent corner ID if we have a LOW.
        switchCornerID(1:NDIM) = blockChild(1:NDIM,eachChild) - 1 !0 or 1.

        !The LEAF neighbor is more refined so by definition the source block 
        !stride cannot be 1 (must be 2,4,...).  srcStride/2 is always valid.
        neghCornerID(1:NDIM,childCounter) = parentCornerID(1:NDIM) + &
             (switchCornerID(1:NDIM) * (srcStride(1:NDIM)/2) * blkSize(1:NDIM))

     end if
  end do

  if (childCounter /= numNegh) then
     call Driver_abortFlash &
          ("[gr_ptGetChildData]: Counter not equal to number of neighbors (1)")
  end if


  !-------------------------------------------------------------------
  !Block / process IDs
  !-------------------------------------------------------------------
  !The neighboring blocks are more refined and may exist on 
  !our processor or on a remote processor.  We store the block 
  !handle of our neighboring block which exists at the same refinement 
  !(i.e. the parent).  Then we use the "child" PARAMESH data 
  !structure to find the block and process IDs of our LEAF
  !block neighbors.


  negh(REFLEVELDIF,1:numNegh) = 1  !We always know the refinement level.
  blkHandle = NONEXISTENT  !Initialize, so we can detect gr_getBlkHandle failure

  !Recall that neghBlkID and neghProcID give our neighbor which is a 
  !parent block, but we need to know the LEAF neighbors.
  call gr_getBlkHandle(neghBlkID, neghProcID, blkHandle)

  if (blkHandle /= NONEXISTENT) then
     
     !The left most neighbor is the first neighbor.  Therefore, we 
     !can obtain all other neighbors by moving to the right in each dimension.
     childCounter = 0
     do eachChild = 1, USABLE_CHILDREN
        

        if (validChild(eachChild) .eqv. .true.) then 
     
           childCounter = childCounter + 1
           negh(BLKID,childCounter) = child(BLKID,eachChild,blkHandle)
           negh(BLKPROC,childCounter) = child(BLKPROC,eachChild,blkHandle)

#ifdef DEBUG_GRIDMAPPARTICLES
           !Check child for junk.  Crash if we find junk.
           if ( (negh(BLKID,childCounter) <= 0) .or. &
                (negh(BLKID,childCounter) > MAXBLOCKS) .or. &
                (negh(BLKPROC,childCounter) < 0) .or. & 
                (negh(BLKPROC,childCounter) >= gr_meshNumProcs) ) then
              
              print *, "Processor:", gr_meshMe, &
                   "encountered junk in PARAMESH child data structure, BlockID:", &
                   negh(BLKID,childCounter), "ProcID:", negh(BLKPROC,childCounter), &
                   "... reverting to cornerID."
              
              call Driver_abortFlash("Invalid neighbor block/proc data 2")
              !Another option is to rely on the cornerID (commented out).
              !negh(BLKID,childCounter) = NONEXISTENT
              !negh(BLKPROC,childCounter) = NONEXISTENT
              !call gr_ptSearchBlk(neghCornerID(:,childCounter),negh(:,childCounter))
           end if
#endif !DEBUG_GRIDMAPPARTICLES
           
           
        end if
     end do

     if (childCounter /= numNegh) then
        call Driver_abortFlash &
             ("[gr_ptGetChildData]: Counter not equal to number of neighbors (2)")
     end if

  else
     !If we don't find a usable block handle, we can still rely on corner ID.
     print *, "Processor:", gr_meshMe, & 
          "block handle not found (2) ... shouldn't happen."
  end if

end subroutine gr_ptGetChildData


!Pass all state so that I can test this subroutine independently.
subroutine gr_ptFindValidChildren(maxDims, maxChildren, blockChild, &
     guardCellID, dims, actualChildren, validChild)

  implicit none

  integer, intent(IN) :: maxDims, maxChildren
  integer, dimension(maxDims,maxChildren), intent(IN) :: blockChild
  integer, dimension(maxDims), intent(IN) :: guardCellID
  integer, intent(IN) :: dims, actualChildren
  logical, dimension(actualChildren), intent(OUT) :: validChild

  integer :: eachChild, eachAxis
  
  validChild(:) = .true.
  do eachChild = 1, actualChildren

     do eachAxis = 1, dims

        !We can only eliminate children when they are not in the center.
        if ( (guardCellID(eachAxis) == LEFT_EDGE .and. & 
             blockChild(eachAxis,eachChild) == LOW) .or. &
             (guardCellID(eachAxis) == RIGHT_EDGE .and. &
             blockChild(eachAxis,eachChild) == HIGH) ) then
           validChild(eachChild) = .false.
        end if

     end do
  end do

end subroutine gr_ptFindValidChildren
