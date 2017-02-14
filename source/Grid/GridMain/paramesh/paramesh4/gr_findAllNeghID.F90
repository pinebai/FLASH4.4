!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_findAllNeghID
!!
!! NAME
!!
!!  gr_findAllNeghID
!!
!! SYNOPSIS
!!
!!  call gr_findAllNeghID(     integer(IN)  :: blockID,
!!             type(AllBlockRegions_t)(OUT) :: surrBlksSummary)
!!
!! DESCRIPTION
!!
!!  For a given block, find the processor ID and blockID for all of its
!!  geometric neighbors.
!!
!!  It is assumed that blockid indicates a leaf block.
!!  Geometric neighbors are leaf blocks adjacent to the block given by blockid
!!  in all directions, including diagonal (edge and corner) neighbors.
!!
!!
!! ARGUMENTS
!!
!!   blockID : ID of block in current processor
!!   surrBlksSummary : Derived data type holding information about the block's surrounding
!!                     blocks.
!!
!! NOTES
!!
!!  The concept of geometric neighbor is different from the concept of neighbor used
!!  internally within PARAMESH.
!!    *  Geometric neighbors are leaf blocks, i.e., always have the finest resolution at
!!       given locations.
!!    *  PARAMESH neighbors are neighboring blocks that exist at the same refinement
!!       level in the mesh hierarchy as the referencing block.  
!!       PARAMESH arrays like surr_blks and neigh refer to these same-level neighbors.
!!       PARAMESH same-level neighbors may have different nodetypes.
!!  Translation from the concept of same-level neighbor to that of geometric neighbor
!!  is the purpose of this routine.
!!
!!  The derived type (AllBlockRegions_t) is declared in module gr_interfaceTypeDecl.
!!
!!  The structure of the surrBlksSummary argument can be considered equivalent to
!!  a combination of
!!           integer(OUT) :: surrBlksNumNegh(3, 1:1+2*K2D, 1:1+2*K3D)
!!           integer(OUT) :: surrBlksDetails(BLKNO:TYPENO,2^(NDIM-1), 3, 1:1+2*K2D, 1:1+2*K3D)
!!
!!   BLKNO,TYPENO,NDIM,K2D,K3D are defined in Flash.h or constants.h.
!!
!! SEE ALSO
!!
!!  gr_interfaceTypeDecl
!!***

subroutine gr_findAllNeghID(blockID, surrBlksSummary)

#include "constants.h"
#include "Flash.h"

  use Grid_data, ONLY : gr_meshNumProcs, gr_meshMe
  use Driver_interface, ONLY : Driver_abortFlash
  use tree, ONLY : surr_blks, parent, child, maxblocks_tr, bnd_box
  use paramesh_dimensions, only : maxblocks_alloc
  use gr_interface, ONLY : gr_findWhichChildren
  use gr_interfaceTypeDecl, ONLY: AllBlockRegions_t !for the definition of surrBlksSummary.

  implicit none
  interface
     subroutine gr_getParentsGuardCellView(srcBlkID, srcGuardCellID, parentGuardCellID)
       integer, intent(IN) :: srcBlkID
       integer, dimension(MDIM), intent(IN) :: srcGuardCellID
       integer, dimension(MDIM), intent(OUT) :: parentGuardCellID
     end subroutine gr_getParentsGuardCellView
  end interface
  interface
     subroutine print_info(blk, proc, blockID, myPE, i, j, k)
       integer, intent(IN) :: blk, proc, blockID, myPE, i, j, k
     end subroutine print_info
  end interface


  integer, intent(IN) :: blockID
  type (AllBlockRegions_t), intent(OUT) :: surrBlksSummary

  integer,dimension(BLKNO:TYPENO) :: negh_prop
  integer :: blkHandle, proc, blk
  integer :: i,j,k, ENDK, ENDJ, ENDI, neghCount, numNegh, allCenters
  integer, dimension(2**(NDIM-1)) :: childID
  integer, dimension(MDIM) :: parentGuardCellID, guardCellID

  if (maxblocks_tr /= maxblocks_alloc) then
     if (gr_meshMe == 0) print *, "Paramesh datastructures inconsistent size!"
  end if
  

  !! First determine the correct endpoints of all the loops.
  ENDI = 3
  ENDJ = max(1,K2D*3)
  ENDK = max(1,K3D*3)

  !! And initialize all the neghbour id's
  childID = NONEXISTENT
  numNegh = NONEXISTENT

  !! Loop over all possible neighbor configurations
  allCenters = 2**NDIM
  kAxisNegh: do k = 1, ENDK
     jAxisNegh: do j = 1, ENDJ
        iAxisNegh: do i = 1, ENDI

           !Initialize all neighbor details in region (i,j,k) to NONEXISTENT.
           surrBlksSummary % regionInfo(i,j,k) % numNegh = NONEXISTENT
           surrBlksSummary % regionInfo(i,j,k) &
                % details(BLKNO:TYPENO, 1:2**(NDIM-1)) = NONEXISTENT


           guardCellID(IAXIS) = i
           guardCellID(JAXIS) = j
           guardCellID(KAXIS) = k


           !! if none of the dimensions has LEFT_EDGE or RIGHT_EDGE
           !! the region is interior of the block. This routine is
           !! not interested in that region
           if((i*j*k) /= allCenters) then 

              !May not be a block at this position.  We may be at 
              !the edge of the domain, and so surr_blks at i,j,k contains 
              !the external boundary conditions.
              negh_prop(:) = surr_blks(:,i,j,k,blockID)



              !First check for an external boundary.
              if (negh_prop(BLKNO) <= PARAMESH_PHYSICAL_BOUNDARY) then

                 !This is the case for e.g. an external reflecting boundary.
                 !...Just copy the external boundary conditions.
                 numNegh = 0
                 surrBlksSummary % regionInfo(i,j,k) &
                      % details(BLKNO:TYPENO,1) = negh_prop(BLKNO:TYPENO)



              !Now check for wacky values in surr_blks.
              else if ( &
                   (negh_prop(BLKNO) < -1) .or. &
                   (negh_prop(BLKNO) == 0) .or. &
                   (negh_prop(BLKNO) > MAXBLOCKS) .or. &
                   (negh_prop(PROCNO) < -1) .or. &
                   (negh_prop(PROCNO) >= gr_meshNumProcs) &
                   ) then

                 print *, gr_meshMe, negh_prop(:)
                 call Driver_abortFlash("Unexpected surr_blks values")
                 

                 
              else if (negh_prop(BLKNO) == NONEXISTENT) then
                 !(Assumes NONEXISTENT is equal to -1).

                 !! If this query returns NONEXISTENT, that means that the
                 !! the neighbor is at a lower resolutioon. Now we need to
                 !! locate parent's neighbor. And here there will be only one
                 !! neighbor along this face/edge/corner

                 numNegh = 1
                 proc = parent(PROCNO,blockID)
                 blk = parent(BLKNO,blockID)



                 !! find the parent ID from metadata/cached metadata
                 blkHandle = -1
                 call gr_getBlkHandle(blk, proc, blkHandle)
                 if (blkHandle == -1) then
                    call print_info(blk, proc, blockID, gr_meshMe, i, j, k)
                    print *, "(1) No block handle.... increase maxblocks_alloc"
                    call Driver_abortFlash("Damn 1")
                 end if

                 if ((blkHandle < 1) .or. (blkHandle > maxblocks_alloc)) then                    
                    call print_info(blk, proc, blockID, gr_meshMe, i, j, k)
                    print *, "(1) Garbage block handle:", blkHandle
                    call Driver_abortFlash("Damn 2")
                 end if


                 !We need to translate the guard cell ID from the source blocks 
                 !perspective to that of its parent.
                 call gr_getParentsGuardCellView(blockID, guardCellID, &
                      parentGuardCellID)


                 surrBlksSummary % regionInfo(i,j,k) & 
                      % details(BLKNO:PROCNO,1) = &
                      surr_blks(BLKNO:PROCNO,parentGuardCellID(IAXIS), &
                      parentGuardCellID(JAXIS),parentGuardCellID(KAXIS), &
                      blkHandle)


                 !! save the information about relative resolution
                 surrBlksSummary % regionInfo(i,j,k) &
                      % details(TYPENO,1) = negh_prop(TYPENO)


              else if (negh_prop(TYPENO) == PARENT_BLK) then

                 !! PARENT_BLK indicates that there may be more than one
                 !! block in the neighborhood here, and they will be 
                 !! at a higher resolution.

                 numNegh = 2 ** (count(guardCellID(1:NDIM) == CENTER))
                 proc = negh_prop(PROCNO)
                 blk = negh_prop(BLKNO)

                 blkHandle = -1
                 call gr_getBlkHandle(blk, proc, blkHandle)
                 if (blkHandle == -1) then
                    call print_info(blk, proc, blockID, gr_meshMe, i, j, k)
                    print *, "(2) No block handle.... increase maxblocks_tr"
                    call Driver_abortFlash("Damn 3")
                 end if

                 if ((blkHandle < 1) .or. (blkHandle > maxblocks_tr)) then                    
                    call print_info(blk, proc, blockID, gr_meshMe, i, j, k)
                    print *, "(2) Garbage block handle:", blkHandle
                    call Driver_abortFlash("Damn 4")
                 end if


                 !! Get the places of all children of interest 
                 call gr_findWhichChildren(numNegh, guardCellID, childID)


                 !! Now for each child, find its ID and place it in 
                 !! the data structure to return
                 do neghCount = 1, numNegh
                    surrBlksSummary % regionInfo(i,j,k) &
                         % details(BLKNO:PROCNO,neghCount) = &
                         child(BLKNO:PROCNO,childID(neghCount),blkHandle)
                 end do

                 !! save the information about relative resolution
                 surrBlksSummary % regionInfo(i,j,k) &
                      % details(TYPENO,1:numNegh) = negh_prop(TYPENO)

              else

                 !! This is the situation when the neighbor is at the
                 !! same level of resulution. There is only one neighbor
                 !! very simply found.
                 numNegh = 1
                 surrBlksSummary % regionInfo(i,j,k) &
                      % details(BLKNO:TYPENO,1) = negh_prop(BLKNO:TYPENO)

              end if
           end if !! end of the condition to check if this is interior


           !Important!  We need to know how many neighbors are available.
           !This tells us what data to access in our surrBlksSummary datatype.
           surrBlksSummary % regionInfo(i,j,k) % numNegh = numNegh

        end do iAxisNegh
     end do jAxisNegh
  end do kAxisNegh

end subroutine gr_findAllNeghID



subroutine gr_getParentsGuardCellView(srcBlkID, srcGuardCellID, parentGuardCellID)

  use tree, only : which_child

  implicit none
  integer, intent(IN) :: srcBlkID
  integer, dimension(MDIM), intent(IN) :: srcGuardCellID
  integer, dimension(MDIM), intent(OUT) :: parentGuardCellID

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
  integer :: eachAxis, child


  parentGuardCellID = 1 !Required for when NDIM < MDIM.
  child = which_child(srcBlkID)

  !Obtain the source blocks position in the parent block.
  !e.g. child=5 is at location: LOW,LOW,HIGH.
  !We use this location to figure out how the source block's 
  !neighbor appears with respect to the parent.

  do eachAxis = 1, NDIM
     if ( (blockChild(eachAxis,child) == LOW) .and. &
          (srcGuardCellID(eachAxis) == LEFT_EDGE) ) then
        parentGuardCellID(eachAxis) = LEFT_EDGE

     else if ( (blockChild(eachAxis,child) == HIGH) .and. &
          (srcGuardCellID(eachAxis) == RIGHT_EDGE) ) then
        parentGuardCellID(eachAxis) = RIGHT_EDGE

     else
        parentGuardCellID(eachAxis) = CENTER

     end if
  end do

end subroutine gr_getParentsGuardCellView


subroutine print_info(blk, proc, blockID, myPE, i, j, k)

  use tree, only : bnd_box
  implicit none
  integer, intent(IN) :: blk, proc, blockID, myPE, i, j, k
  print *, "Block handle error for target block:", &
       blk, ", proc:", proc, ". My block is:", blockID, &
       "and proc is:", myPE, &
       " and we were trying to find neighbors to guard cell region:", &
       i, j, k, "and my global space is: IAXIS=", &
       bnd_box(LOW:HIGH,IAXIS,blockID), "JAXIS=", &
       bnd_box(LOW:HIGH,JAXIS,blockID), "KAXIS=", &
       bnd_box(LOW:HIGH,KAXIS,blockID)

end subroutine print_info
