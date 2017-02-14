!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/gr_ptFindNegh
!!
!! NAME
!!  gr_ptFindNegh
!!
!! SYNOPSIS
!!
!!  gr_ptFindNegh(integer(IN)  :: srcBlkID,
!!                integer(IN)  :: guardCellID(MDIM),
!!                integer(OUT) :: negh(BLKID:REFLEVELDIF,4),
!!                integer(OUT) :: neghCornerID(MDIM,4),
!!                integer(OUT) :: numNegh)
!!
!! DESCRIPTION
!!
!!  This routine is designed to find information about the neighboring blocks 
!!  to a particular guard cell region.  The caller describes the guard cell region 
!!  by providing the input arguments: srcBlkID and guardCellID.
!!  The guardCellID array must be populated with one of 1,2,3 in each dimension,
!!  which corresponds to the LEFT_EDGE(1), CENTER(2), RIGHT_EDGE(3) regions.  
!!
!!  The neighboring blocks are described by the output arguments negh, neghCornerID 
!!  and numNegh.  The number of neighboring blocks (numNegh) range from 0 up to 1 
!!  (1D simulations), 2(2D simulations), 4(3D simulations).  A 0 is obtained when 
!!  the guard cell region touches a Paramesh physical boundary.  The  
!!  negh and neghCornerID data arrays will only contain valid neigbor information for 
!!  the number of neighbors returned, e.g. when the number of neighbors (numNegh) is 1, 
!!  negh(BLKID,2) will contain invalid data.
!!
!!  The algorithm works by: 
!!
!!      Use the Paramesh surr_blks data structure to find the neighboring 
!!      block of the source block guard cell region.  A single value is returned:
!!         A:  Valid block ID.
!!         B:  Invalid block ID (NONEXISTANT).
!!         C:  Invalid block ID (PARAMESH_PHYSICAL_BOUNDARY).
!!
!!      A NONEXISTANT value indicates there is no valid neighbor at the source 
!!      block's resolution.  Therefore, the neighbor is at a lower resolution, 
!!      and cannot be found directly.  The data can be accessed from surr_blks at
!!      index = block handle, where block handle is specific to the particular 
!!      remote block.
!!
!!      A PARAMESH_PHYSICAL_BOUNDARY value indicates the source block is touching 
!!      an external boundary.
!!      
!!      If the neighbor block ID is valid then one of two situations are possible.
!!         A:  Neighboring block at the same resolution.
!!         B:  Neighboring block is a parent, in which case its children are 
!!             at a higher resolution than the source block.
!!
!!      When we have obtained as much information from the surr_blks structure 
!!      as possible, we will query PARAMESH cached information through a block 
!!      handle.
!!
!!
!! ARGUMENTS
!!
!!               srcBlkID:  The Paramesh leaf block identifier of the source block
!!               guardCellID:  The position of guard cell region relative to 
!!                             the source block.
!!               negh:  An array containing information about the destination block(s).
!!               neghCornerID:  An array containing the corner ID of the 
!!                              destination block(s).
!!               numNegh:  The number of neighbors to a guard cell region
!!
!!***

!Useful to turn on additional assertations in this file by default.
#ifndef DEBUG_GRIDMAPPARTICLES
#define DEBUG_GRIDMAPPARTICLES
#endif

!A debugging mode in which we delete any information about 
!neighboring blocks & process IDs.  Thia means we must rely on
!the corner ID and refinement level of the neighboring blocks during
!the data exchange.
!#define RELY_ON_CORNERID

#include "constants.h"
#include "Flash.h"
#include "gr_ptMapToMesh.h"

subroutine gr_ptFindNegh(srcBlkID,guardCellID,negh,neghCornerID,numNegh)

  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs, gr_globalDomain
  use tree, ONLY : surr_blks, parent, lrefine, which_child, nodetype, bnd_box, &
       lrefine_max
  use Grid_interface, ONLY : Grid_getBlkCornerID, &
       Grid_getBlkIndexLimits, Grid_getBlkBC
  use gr_ptInterface, ONLY : gr_ptSearchBlk, gr_ptGetChildData
  use gr_interface, ONLY : gr_getBlkHandle
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  interface 
     subroutine gr_ptApparentCornerID(cornerID, guardCellID, faces, apparentCornerID)
       implicit none
       integer, dimension(1:MDIM), intent(IN) :: cornerID, guardCellID
       integer, dimension(LOW:HIGH, MDIM), intent(IN) :: faces
       integer, dimension(1:MDIM), intent(OUT) :: apparentCornerID
     end subroutine gr_ptApparentCornerID
  end interface
  interface
     subroutine gr_ptCheckNeighborsSensible(numNegh, negh, neghCornerID)
       implicit none
       integer, intent(IN) :: numNegh
       integer, dimension(BLKID:REFLEVELDIF,ABSMAXNEGH), intent(IN) :: negh
       integer, dimension(MDIM,ABSMAXNEGH), intent(IN) :: neghCornerID
     end subroutine gr_ptCheckNeighborsSensible
  end interface
  interface
     subroutine gr_ptCheckGuardCellAdjust2(srcBlkID, guardCellID, parentGuardCellID)
       implicit none
       integer, intent(IN) :: srcBlkID
       integer, dimension(MDIM), intent(IN) :: guardCellID, parentGuardCellID
     end subroutine gr_ptCheckGuardCellAdjust2
  end interface

  integer, intent(IN) :: srcBlkID
  integer, dimension(MDIM), intent(IN) :: guardCellID
  integer, dimension(BLKID:REFLEVELDIF,ABSMAXNEGH),intent(OUT):: negh
  integer, dimension(MDIM,ABSMAXNEGH),intent(OUT) :: neghCornerID
  integer, intent(OUT) :: numNegh

  integer, dimension(LOW:HIGH, MDIM) :: blkLimits, blkLimitsGC, faces, ignoreMe
  integer, dimension(MDIM) :: srcCornerID, srcApparentCornerID, &
       parentCornerID, parentApparentCornerID, parentGuardCellID, &
       srcStride, blkSize, neghAtSameRefineCornerID, &
       blockPosition, mystride, srcCornerIDcalculated
  integer :: neghBlkID, neghProcID, neghBlkType, gCell, n
  integer :: blkHandle, nBlksFromDomCorner, eachAxis, idSpan, oddValue
  logical :: changeToCenter


#ifdef DEBUG_GRIDMAPPARTICLES
  !Check that this block actually exists and is a LEAF.
  if (nodetype(srcBlkID) /= LEAF) then
     print *, "Processor:", gr_meshMe, "ERROR! source block is type:", &
          nodetype(srcBlkID)
     call Driver_abortFlash("[gr_ptFindNegh]: Source block is not a LEAF")
  end if
#endif !DEBUG_GRIDMAPPARTICLES



  !Obtain information about this particular block.
  call Grid_getBlkCornerID(srcBlkID, srcCornerID, srcStride)
  call Grid_getBlkBC(srcBlkID, ignoreMe, faces)
  call Grid_getBlkIndexLimits(srcBlkID, blkLimits, blkLimitsGC)
  blkSize = blkLimits(HIGH,:) - blkLimits(LOW,:) + 1


  !! First determine the ID of the neighboring block
  !! If it is at same or higher level of resolution, we get a valid 
  !! value in neghBlkID, if the neighbor is at lower resolution then
  !! the block has no valid neighbor, its parent has one, and the 
  !! returned value is NONEXISTENT. For the purpose of active particles
  !! the parent's neighbor is a valid neighbor for the child if the child
  !! doesn't have a neighbor of its own.
  neghBlkID = surr_blks(BLKID,guardCellID(IAXIS),&
       1+(guardCellID(JAXIS)-1)*K2D,&
       1+(guardCellID(KAXIS)-1)*K3D,srcBlkID)
  neghProcID = surr_blks(BLKPROC,guardCellID(IAXIS),1+(guardCellID(JAXIS)-1)*K2D,&
       1+(guardCellID(KAXIS)-1)*K3D,srcBlkID)
  neghBlkType = surr_blks(BLKTYPE,guardCellID(IAXIS),1+(guardCellID(JAXIS)-1)*K2D,&
       1+(guardCellID(KAXIS)-1)*K3D,srcBlkID)

  negh = -200  !Random initial value for spotting errors.
  neghCornerID = -200



  if (neghBlkID <= PARAMESH_PHYSICAL_BOUNDARY) then
     !This is the case for e.g. an external reflecting boundary.

     numNegh = 0
     return

  else if(neghBlkID == NONEXISTENT) then 

     numNegh = 1     !! There can be only one block that shares the corner, edge or face

     !--------------------------------------------------------------------
     ! Determine the corner ID of the source block's parent.
     !--------------------------------------------------------------------
     parentCornerID = 1 !Necessary for the case where NDIM < MDIM.
     parentGuardCellID = 1 !Necessary for surr_blks call when NDIM < MDIM.

     do eachAxis = 1, NDIM

        !Find the relative position of the source block in the global domain.
        !First, calculate the distance from this block's corner ID to the next block's
        !corner ID.  Next, count how many blocks it takes to reach the 
        !source block from the computational domain corner.
        idSpan = srcStride(eachAxis) * blkSize(eachAxis)
        nBlksFromDomCorner = (srcCornerID(eachAxis)-1) / idSpan

        !The number of blocks will be even or odd, which thus tells us 
        !the location of the source block relative to its parent.
        oddValue = mod(nBlksFromDomCorner,2)  !0 is even, 1 is odd.

        parentCornerID(eachAxis) = srcCornerID(eachAxis) - &
             (oddValue * srcStride(eachAxis) * blkSize(eachAxis))

     end do


     !--------------------------------------------------------------------
     ! Adjust the guard cell region, so it is relative to the parent.
     !-------------------------------------------------------------------- 
     do eachAxis = 1, NDIM

        gCell = guardCellID(eachAxis)        
        changeToCenter = &
             (((gCell == LEFT_EDGE) .and. (srcCornerID(eachAxis) > parentCornerID(eachAxis))) &
             .or. &
             ((gCell == RIGHT_EDGE) .and. (srcCornerID(eachAxis) == parentCornerID(eachAxis))))

        if (changeToCenter .eqv. .true.) then
           parentGuardCellID(eachAxis) = CENTER
        else
           parentGuardCellID(eachAxis) = gCell
        end if
     end do


#ifdef DEBUG_GRIDMAPPARTICLES
     !Check the validity of the above expression.
     call gr_ptCheckGuardCellAdjust2(srcBlkID, guardCellID, parentGuardCellID)
#endif !DEBUG_GRIDMAPPARTICLES


     !--------------------------------------------------------------------
     !Determine the corner ID for the neighbor of the parent.
     !-------------------------------------------------------------------- 
     !The parent of the source block may exist on another processor.  
     !Therefore, we pass the boundary conditions of the source block (faces) 
     !into gr_ptApparentCornerID.
     !-------------------------------------------------------------------- 
     call gr_ptApparentCornerID(parentCornerID, parentGuardCellID, faces, &
          parentApparentCornerID)


     !parentGuardCellID(eachAxis)-2) is the direction: (-1,0,+1)
     neghCornerID(1:NDIM,1) = parentApparentCornerID(1:NDIM) + &
          ((2*srcStride(1:NDIM))*blkSize(1:NDIM)*(parentGuardCellID(1:NDIM)-2))


     !--------------------------------------------------------------------
     ! Now that we know where the neighbor block touches our parent block 
     ! we can find the neighbor's block and process ID.
     !-------------------------------------------------------------------- 
     blkHandle = NONEXISTENT  !Initialize, so we can detect gr_getBlkHandle failure
     call gr_getBlkHandle(parent(BLKNO,srcBlkID), parent(PROCNO,srcBlkID), blkHandle)

     negh(REFLEVELDIF,1) = -1  !We always know the refinement level.
     if (blkHandle /= NONEXISTENT) then

        negh(BLKID,1) = &
             surr_blks(BLKID,parentGuardCellID(IAXIS),parentGuardCellID(JAXIS),&
             parentGuardCellID(KAXIS),blkHandle)
        negh(BLKPROC,1) = &
             surr_blks(BLKPROC,parentGuardCellID(IAXIS),parentGuardCellID(JAXIS),&
             parentGuardCellID(KAXIS),blkHandle)


#ifdef DEBUG_GRIDMAPPARTICLES
           !Check surr_blks for junk, and revert to cornerID if we find junk.
           if ( (negh(BLKID,1) <= 0) .or. &
                (negh(BLKID,1) > MAXBLOCKS) .or. &
                (negh(BLKPROC,1) < 0) .or. & 
                (negh(BLKPROC,1) >= gr_meshNumProcs) ) then

              !We had this error occur when the data in surr_blks was out of date.
              !The new call to gr_ensureValidNeighborInfo in Grid_mapParticlesToMesh should 
              !ensure this is now fixed.

              blockPosition(1:NDIM) = &
                   nint( (bnd_box(LOW,1:NDIM,srcBlkID) - gr_globalDomain(LOW,1:NDIM)) / &
                         (bnd_box(HIGH,1:NDIM,srcBlkID) - bnd_box(LOW,1:NDIM,srcBlkID)) )
              myStride(1:NDIM) = 2**(lrefine_max - lrefine(srcBlkID))
              srcCornerIDcalculated(1:NDIM) = &
                   blockPosition(1:NDIM) * (myStride(1:NDIM) * blkSize(1:NDIM)) + 1
              
              print *, "Processor:", gr_meshMe, &
                   "encountered following junk in PARAMESH surr_blks data structure, BlockID:", &
                   negh(BLKID,1), "ProcID:", negh(BLKPROC,1), &
                   "for source block:", srcBlkID, " at refine:", lrefine(srcBlkID), &
                   ".  Some info: parent (BLK, PROC):",  &
                   parent(BLKNO,srcBlkID), parent(PROCNO,srcBlkID), &
                   "src corner ID:", srcCornerID(1:NDIM), ", parent corner ID:", &
                   parentCornerID(1:NDIM), ", negh corner ID:", neghCornerID(1:NDIM,1), &
                   ", guardCellID:", guardCellID(1:NDIM), ", parent guardCellID:", &
                   parentGuardCellID(1:NDIM), ", srcStride:", srcStride(1:NDIM), & 
                   ", blkSize:", blkSize(1:NDIM), ", which_child:", which_child(srcBlkID), &
                   ", predicted src corner ID:", srcCornerIDCalculated(1:NDIM)        

              call Driver_abortFlash("Invalid neighbor block/proc data 1")
              !Another option is to rely on the cornerID (commented out).
              !negh(BLKID,1) = NONEXISTENT
              !negh(BLKPROC,1) = NONEXISTENT
              !call gr_ptSearchBlk(neghCornerID(:,1),negh(:,1))

           end if
#endif !DEBUG_GRIDMAPPARTICLES


     else
        !If we don't find a usable block handle, we can still rely on corner ID.
        print *, "Processor:", gr_meshMe, & 
             "block handle not found (1) ... shouldn't happen."
     end if

  else

     !The neighbor has a known identifier.

     !--------------------------------------------------------------------
     !We must operate on the source block corner ID to uphold periodicity.
     !--------------------------------------------------------------------
     call gr_ptApparentCornerID(srcCornerID, guardCellID, faces, &
          srcApparentCornerID)

     neghAtSameRefineCornerID(1:NDIM) = srcApparentCornerID(1:NDIM) + &
             ((guardCellID(1:NDIM)-2)*srcStride(1:NDIM)*blkSize(1:NDIM))


     if (neghBlkType == LEAF) then
        !! If the neighbor is a leaf block, then clearly it is at 
        !! the same refinement level, so the difference in level
        !! from the source block is 0, and there is only one neighbor
        !! along the current corner, edge or face.
        numNegh = 1
        negh(BLKID,1) = neghBlkID
        negh(BLKPROC,1) = neghProcID
        negh(REFLEVELDIF,1) = 0
        neghCornerID(1:NDIM,1) = neghAtSameRefineCornerID(1:NDIM)

     else  !! If the neighbor is a parent then, it has been refined.
        !! the number of LEAF neighbors then is 1, 2 or 4 depeding upon
        !! whether it is a corner, edge or a face. A 1D problem can only
        !! have a corner, a 2D corner and edge, and a 3D problem can have
        !! all three. If the Parent neighbor is on the current processor, then
        !! all of its children can be determinisitically found, otherwise they
        !! have to be searched using cornerID. numNegh=1 implies a corner, so
        !! one neighbor, numNegh=2 implies an edge, so two neighbors, and numNegh=4
        !! implies a face, so 4 neighbors

        !This subroutine updates negh, neghCornerID and numNegh.
        call gr_ptGetChildData(guardCellID, blkSize, neghAtSameRefineCornerID, srcStride, &
             neghBlkID, neghProcID, negh, neghCornerID, numNegh)

     end if

  end if

#ifdef DEBUG_GRIDMAPPARTICLES
  call gr_ptCheckNeighborsSensible(numNegh, negh, neghCornerID)
#endif !DEBUG_GRIDMAPPARTICLES


  !This is a mode for debugging.  It removes whatever we have stored 
  !for the block & process IDs.  This means we have to rely on the cornerID.
  !NOTE: This is only an option if we are using the Sieve implementation.
#ifdef RELY_ON_CORNERID
  !We want to retain the refinement level information in negh.
  negh(BLKID,:) = NONEXISTENT
  negh(BLKPROC,:) = NONEXISTENT

  do n = 1, numNegh
     if (negh(BLKPROC,n) == NONEXISTENT) then
        !It is possible that the neighbor exists on our processor after all.
        !We will do a linear search for the corner ID.
        call gr_ptSearchBlk(neghCornerID(:,n),negh(:,n))
     end if
  end do
#endif !RELY_ON_CORNERID

end subroutine gr_ptFindNegh


!-----------------------------------------------------------------------
! Helper routines (START)
!-----------------------------------------------------------------------
subroutine gr_ptApparentCornerID(cornerID, guardCellID, faces, apparentCornerID)

  use Grid_interface, ONLY : Grid_getGlobalIndexLimits

  implicit none
  integer, dimension(1:MDIM), intent(IN) :: cornerID, guardCellID
  integer, dimension(LOW:HIGH, MDIM), intent(IN) :: faces
  integer, dimension(1:MDIM), intent(OUT) :: apparentCornerID
  integer, dimension(1:MDIM) :: globalIndexLimits
  integer :: eachAxis

  apparentCornerID = 1
  call Grid_getGlobalIndexLimits(globalIndexLimits)

  do eachAxis = 1, NDIM

     !If there are PERIODIC boundary conditions, we must operate on the 
     !corner ID of the block to uphold periodicity.
     if ( ((guardCellID(eachAxis) == LEFT_EDGE ) .and. (faces(LOW ,eachAxis) == PERIODIC))&
          .or.&
          ((guardCellID(eachAxis) == RIGHT_EDGE) .and. (faces(HIGH,eachAxis) == PERIODIC)) ) then

        apparentCornerID(eachAxis) = cornerID(eachAxis) - &
             ((guardCellID(eachAxis)-2) * globalIndexLimits(eachAxis))

     else
        apparentCornerID(eachAxis) = cornerID(eachAxis)
     end if

  end do

end subroutine gr_ptApparentCornerID
!-----------------------------------------------------------------------
! Helper routines (FINISH)
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! Checking routines (START)   (Only take input args!)
!-----------------------------------------------------------------------
subroutine gr_ptCheckNeighborsSensible(numNegh, negh, neghCornerID)

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkCornerID
  use Grid_data, ONLY : gr_meshMe

  implicit none

  integer, intent(IN) :: numNegh
  integer, dimension(BLKID:REFLEVELDIF,ABSMAXNEGH), intent(IN) :: negh
  integer, dimension(MDIM,ABSMAXNEGH), intent(IN) :: neghCornerID

  integer, dimension(MDIM) :: neghCornerIDCHECK, dummyStride
  integer :: n


  !Check that our calculated corner IDs are sensible: 
  do n = 1, numNegh

     if ( negh(BLKPROC,n) == gr_meshMe ) then

        !This verifies that our technique for finding corner IDs is correct.
        !We can only check if the neighboring blocks are on the same processor.
        call Grid_getBlkCornerID(negh(BLKID,n), neghCornerIDCHECK, dummyStride)

        if (any(neghCornerID(1:NDIM,n) /= neghCornerIDCHECK(1:NDIM))) then
           print *, "Processor:", gr_meshMe, ", CornerID: predicted=", &
                neghCornerID(1:NDIM,n), ", actual=", neghCornerIDCHECK(1:NDIM)            
           call Driver_abortFlash &
                ("[gr_ptCheckNeighborsSensible]: Incorrect Corner ID prediction (1)")
        end if

     else

        !We can only check that the corner ID is not a crazy value
        if (any(neghCornerID(1:NDIM,n) < 0)) then
           print *, "Processor:", gr_meshMe, ", CornerID: predicted=", &
                neghCornerID(1:NDIM,n)
           call Driver_abortFlash &
                ("[gr_ptCheckNeighborsSensible]: Incorrect Corner ID prediction (2)")
        end if

     end if

     !Possible refinement levels are -1, 0, +1.
     if ( (negh(REFLEVELDIF,n) < -1) .or. (negh(REFLEVELDIF,n) > 1) ) then
        print *, "Refinement level not possible:", negh(REFLEVELDIF,n)
        call Driver_abortFlash &
             ("[gr_ptCheckNeighborsSensible]: Incorrect refinement level")
     end if

  end do

end subroutine gr_ptCheckNeighborsSensible

!!*********************************************************************************

subroutine gr_ptCheckGuardCellAdjust2(srcBlkID, guardCellID, parentGuardCellID)

  use tree, only : which_child
  use Grid_data, only : gr_meshMe

  implicit none

  integer, intent(IN) :: srcBlkID
  integer, dimension(MDIM), intent(IN) :: guardCellID, parentGuardCellID

  integer, parameter :: TOTAL_CHILDREN = 2**MDIM
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
  integer, dimension(MDIM) :: parentGuardCellIDTmp
  integer :: eachAxis, child


  parentGuardCellIDTmp = NONEXISTENT
  child = which_child(srcBlkID)

  do eachAxis = 1, NDIM

     if ( (blockChild(eachAxis,child) == LOW) .and. &
          (guardCellID(eachAxis) == LEFT_EDGE) ) then
        parentGuardCellIDTmp(eachAxis) = LEFT_EDGE

     else if ( (blockChild(eachAxis,child) == HIGH) .and. &
          (guardCellID(eachAxis) == RIGHT_EDGE) ) then
        parentGuardCellIDTmp(eachAxis) = RIGHT_EDGE

     else
        parentGuardCellIDTmp(eachAxis) = CENTER

     end if

  end do


  if (any(parentGuardCellIDTmp(1:NDIM) /= parentGuardCellID(1:NDIM))) then
     print *, "Processor:", gr_meshMe, "gets conflicting answers for block:", &
          srcBlkID, "parentGuardCellID - Orig:", parentGuardCellID, &
          "New:", parentGuardCellIDTmp             
  end if


end subroutine gr_ptCheckGuardCellAdjust2
!-----------------------------------------------------------------------
! Checking routines (FINISH)
!-----------------------------------------------------------------------
