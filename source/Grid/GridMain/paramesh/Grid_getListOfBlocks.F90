!!****if* source/Grid/GridMain/paramesh/Grid_getListOfBlocks
!!
!! NAME
!!  Grid_getListOfBlocks
!!
!! SYNOPSIS
!!
!!  Grid_getListOfBlocks(integer(IN)          :: blockType,
!!                       integer(OUT)         :: listofBlocks(MAXBLOCKS), 
!!                       integer(OUT)         :: count,
!!                       integer(IN,optional) :: refinementLevel,
!!                       real(IN,optional)    :: region_bndBox(LOW:HIGH,MDIM)
!!                       logical(IN,optional) :: includePartialBlocks)
!!  
!! DESCRIPTION 
!!  Returns a list and the number of blocks of a specified type on the local processor
!!  This routine can also be used to find blocks that are on a particular boundary.
!!  
!!
!! ARGUMENTS
!!
!!  blockType - specification of block type
!!              For all Grid implementations, valid values are :
!!
!!              ALL_BLKS    all local blocks on a processor
!!
!!              IBDRY_BLKS  blocks that are on physical boundary along IAXIS
!!              JBDRY_BLKS  blocks that are on physical boundary along JAXIS
!!              KBDRY_BLKS  blocks that are on physical boundary along KAXIS
!!              ANY_BDRY_BLKS  blocks that have any of their faces 
!!              on a physical boundary.
!!              ACTIVE_BLKS all currently active blocks, in paramesh
!!              context that means parent and leaf blocks
!!              
!!              values that have meaning only for paramesh are :
!!              LEAF, PARENT_BLK or ANCESTOR  representing
!!              the type of node in the Oct-tree managing the blocks.
!!              REFINEMENT the refinement level
!!              INREGION  All blocks within the region defined by
!!                          the accompanying optional argument
!!                          region_bndBox. If the optional argument
!!                          refinementLevel is also present then the blocks
!!                          are also checked to see if they are at the specified
!!                          refinement level, and only those that are get included
!!                          in the list
!!            
!!              All of these constants are defined in constants.h
!!
!!  listofBlocks - returned array holding the integer block number of all the blocks
!!                 on a local processor of type 'blockType'
!!
!!  count - number of blocks returned in listofBlocks
!!
!!  refinementLevel - requested refinement level, only valid with blockType = REFINEMENT
!!                    of INREGION
!!
!!  region_bndBox - when blocktype is specified as INREGION this argument defines the
!!                  bounding box of the region
!!
!!  includePartialBlocks - this argument is valid only when region_bndBox is present
!!                  when present and true, the blocks that are partially in the specified region
!!                  are included in the returned list, otherwise they are ignored.
!!
!! EXAMPLE
!!   
!!   Consider a 2 dimensional problem with 16 blocks, 4 blocks along IAXIS and 
!!   4 along JAXIS, and they are numbered in lexicographic order as follows, all
!!   on the same processor.
!!
!!    --- --- --- ---
!!   | 1 | 2 | 3 | 4 | 
!!    --- --- --- ---
!!   | 5 | 6 | 7 | 8 | 
!!    --- --- --- ---
!!   | 9 |10 |11 |12 | 
!!    --- --- --- ---
!!   |13 |14 |15 |16 | 
!!    --- --- --- ---
!!    
!!   call Grid_getListOfBlocks(JBDRY_BLKS, listOfBlocks, count)
!!      returns count = 8, and listOfBlocks = <1 2 3 4 13 14 15 16>
!!   
!!   call Grid_getListOfBlocks(ANY_BDRY_BLKS, listOfBlocks, count)
!!     returns count = 12 and listOfBlocks = < 1 2 3 4 5 8 9 12 13 14 15 16 >
!!
!!
!!***


subroutine Grid_getListOfBlocks(blockType, listOfBlocks,count,refinementLevel,&
     region_bndBox, includePartialBlocks)

  use tree, ONLY : nodetype,lnblocks,neigh,lrefine,bnd_box
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_oneBlock

#include "constants.h"
#include "Flash.h"  

  implicit none

  integer, intent(in) :: blockType
  integer,dimension(MAXBLOCKS),intent(out) :: listOfBlocks
  integer,intent(out) :: count
  integer,intent(IN), optional :: refinementLevel
  real, dimension(LOW:HIGH,MDIM), intent(IN), optional :: region_bndBox
  logical, intent(IN), optional :: includePartialBlocks

  integer :: i,j
  logical :: isBnd=.false.
  logical :: doOverlapTest, full_overlap, no_overlap, partialBlocks

  count = 0


  select case (blockType)
  case (IBDRY_BLKS)
     do i = 1,lnblocks 
        if((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK)) then
           isBnd = (neigh(1,1,i) .LE. PARAMESH_PHYSICAL_BOUNDARY)
           isBnd = isBnd.or.(neigh(1,2,i).LE.PARAMESH_PHYSICAL_BOUNDARY)
           if(isBnd) then
              count = count+1
              listOfBlocks(count)=i
           end if
        end if
     end do
  case (JBDRY_BLKS)
     do i = 1,lnblocks
        if((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK)) then
           isBnd = (neigh(1,3,i) .LE.PARAMESH_PHYSICAL_BOUNDARY)
           isBnd = isBnd.or.(neigh(1,4,i) .LE.PARAMESH_PHYSICAL_BOUNDARY)
           if(isBnd) then
              count = count+1
              listOfBlocks(count)=i
           end if
        end if
     end do
  case(KBDRY_BLKS)
     do i = 1,lnblocks
        if((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK)) then
           isBnd = (neigh(1,5,i) .LE.PARAMESH_PHYSICAL_BOUNDARY)
           isBnd = isBnd.or.(neigh(1,6,i) .LE.PARAMESH_PHYSICAL_BOUNDARY)
           if(isBnd) then
              count = count+1
              listOfBlocks(count)=i
           end if
        end if
     end do
  case(ANY_BDRY_BLKS)
     do i = 1,lnblocks
        if((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK)) then
           isBnd = (neigh(1,1,i) .LE.PARAMESH_PHYSICAL_BOUNDARY)
           isBnd = isBnd.or.(neigh(1,2,i) .LE.PARAMESH_PHYSICAL_BOUNDARY)
           isBnd = isBnd.or.(neigh(1,3,i) .LE.PARAMESH_PHYSICAL_BOUNDARY)
           isBnd = isBnd.or.(neigh(1,4,i) .LE.PARAMESH_PHYSICAL_BOUNDARY)
           isBnd = isBnd.or.(neigh(1,5,i) .LE.PARAMESH_PHYSICAL_BOUNDARY)
           isBnd = isBnd.or.(neigh(1,6,i) .LE.PARAMESH_PHYSICAL_BOUNDARY)
           if(isBnd) then
              count = count+1
              listOfBlocks(count)=i
           end if
        end if
     end do
  case(ACTIVE_BLKS)
     do i = 1,lnblocks
        if((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK)) then
           count=count+1
           listOfBlocks(count)=i
        end if
     end do
  case(ALL_BLKS)
     do i = 1,lnblocks
        count=count+1
        listOfBlocks(count)=i
    end do

 case(TRAVERSED)
    do i = 1,lnblocks
       if(gr_oneBlock(i)%blockType==TRAVERSED) then
          count=count+1
          listOfBlocks(count)=i
       end if
     end do

 case(TRAVERSED_AND_ACTIVE)
    do i = 1,lnblocks
       if(gr_oneBlock(i)%blockType==TRAVERSED_AND_ACTIVE) then
          count=count+1
          listOfBlocks(count)=i
       end if
     end do

  case(REFINEMENT)
     if (present(refinementLevel)) then
        ! Deal with the special case of refinement
        do i = 1, lnblocks
           if (lrefine(i) == refinementLevel) then
              count = count + 1
              listOfBlocks(count) = i
           end if
        end do
     else
        call Driver_abortFlash("[Grid_getListofBlocks] with blockType REFINEMENT optional argument refinementlevel must be present")
     end if

  case(INREGION)
     if(present(region_bndBox)) then
        if (present(includePartialBlocks)) then
           partialBlocks = includePartialBlocks
        else
           partialBlocks = .false.
        end if
        do i = 1,lnblocks
           if(nodetype(i)==LEAF) then !! only operate on leaf blocks
              !! if there is need to check refinement level, do so
              if (present(refinementLevel)) then
                 doOverlapTest = (lrefine(i) == refinementLevel)
              else
                 doOverlapTest = .true.
              end if

              if (doOverlapTest) then
                 !! now find if the block has any overlap
                 full_overlap=.true.
                 no_overlap=.false.
                 do j = 1,NDIM
                    full_overlap=full_overlap.and. (bnd_box(LOW,j,i).ge.region_bndBox(LOW,j).and.&
                         (bnd_box(HIGH,j,i).le.region_bndBox(HIGH,j)))
                    no_overlap = no_overlap.or. (bnd_box(HIGH,j,i).le.region_bndBox(LOW,j)).or.&
                         (bnd_box(LOW,j,i).ge.region_bndBox(HIGH,j))
                 end do
                 if((full_overlap).or.(.not.no_overlap.and.partialBlocks))then
                    count=count+1
                    listOfBlocks(count)=i
                 end if
              end if
           end if
        end do
     else
        call Driver_abortFlash("[Grid_getListofBlocks] with blockType INREGION optional argument region_bndBox must be present")
     end if
  case default
     do i = 1,lnblocks
        if(nodetype(i) == blockType) then
           count = count+1
           listOfBlocks(count) = i
        end if
     end do
  end select


  return
end subroutine Grid_getListOfBlocks
