!!****if* source/Grid/GridMain/UG/Grid_getListOfBlocks
!!
!! NAME
!!  Grid_getListOfBlocks
!!
!! SYNOPSIS
!!
!!  Grid_getListOfBlocks(integer(IN)  :: blockType,
!!                       integer(OUT) :: listOfBlocks(MAXBLOCKS), 
!!                       integer(OUT) :: count,
!!                       integer(IN,optional) :: refinementLevel,
!!                       real(IN,optional)    :: region_bndBox(LOW:HIGH,MDIM)
!!                       logical(IN,optional) :: includePartialBlocks)
!!  
!! DESCRIPTION 
!! 
!!  This routine is used to find list of blocks that satisfy certain criterion specified
!!  by blockType. Since UG has only one block per processor, the list consists of
!!  either "0" or "1" block. Some of the valid input arguments don't apply to UG,
!!  for those both output arguments are returned with value zero.
!!
!! ARGUMENTS
!!
!!  blockType - The specification for selecting blocks for the list 
!!              Valid values for Uniform grid are :
!!              ALL_BLKS    all blocks on a local proc.  In this case always 1, same as LEAF  
!! 
!!              IBDRY_BLKS  blocks that are on physical boundary along IAXIS
!!              JBDRY_BLKS  blocks that are on physical boundary along JAXIS
!!              KBDRY_BLKS  blocks that are on physical boundary along KAXIS
!!              ANY_BDRY_BLKS  blocks that have any of their faces on the physical
!!                              boundaries.
!!
!!              LEAF, PARENT_BLK and ANCESTOR are valid values when using Paramesh
!!              they are ignored in UG
!!              INREGION  All blocks within the region defined by
!!                          the accompanying optional argument
!!                          region_bndBox. If the optional argument
!!                          refinementLevel is also present then the blocks
!!                          are also checked to see if they are at the specified
!!                          refinement level, and only those that are get included
!!                          in the list
!!            
!!
!!  listOfBlocks - return 1 if the block is the right type else return 0
!!
!!  count - return 1 if the block is the right type else return 0
!!
!!  refinementLevel - requested refinement level, only valid with blockType = REFINEMENT
!!                    of INREGION
!!
!!  region_bndBox - when blocktype is specified as INREGION this argument defines the
!!                  bounding box of the region
!!  includePartialBlocks - this argument is valid only when region_bndBox is present
!!                  when present and true, the blocks that are partially in the specified region
!!                  are included in the returned list, otherwise they are ignored.
!!
!! EXAMPLE
!!
!!   Consider a 2d domain on 16 processors numbered as in the figure below.
!!
!!    --- --- --- ---
!!   | 0 | 1 | 2 | 3 | 
!!    --- --- --- ---
!!   | 4 | 5 | 6 | 7 | 
!!    --- --- --- ---
!!   | 8 | 9 |10 |11 | 
!!    --- --- --- ---
!!   |12 |13 |14 |15 | 
!!    --- --- --- ---
!!
!!   The table below lists the processors that will return non-zero values in count
!!   and listOfBlocks for each meaningful UG BlockType value
!!
!!    IBDRY_BLKS - 0,3,4,7,8,11,12,15
!!    JBDRY_BLKS - 0,1,2,3,12,13,14,15
!!    KBDRY_BLKS - all blocks since the third dimension doesn't really exist
!!    ANY_BDRY_BLKS - 0,1,2,3,4,7,8,11,12,13,14,15
!!
!!***


subroutine Grid_getListOfBlocks(blockType, listOfBlocks,count,refinementLevel,&
     region_bndBox, includePartialBlocks)


  use Grid_data, ONLY : gr_axisMe, gr_axisNumProcs, gr_blockType
  use Logfile_interface, ONLY : Logfile_stampMessage

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockType
  integer,dimension(MAXBLOCKS),intent(out) :: listOfBlocks
  integer,intent(out) :: count
  integer,intent(IN),optional :: refinementLevel
  real, dimension(LOW:HIGH,MDIM), intent(IN), optional :: region_bndBox
  logical, intent(IN), optional :: includePartialBlocks

  integer :: i

!=============================================================================
  count = 0          
  listOfBlocks=0
  select case (blockType)
  case (IBDRY_BLKS)
     if((gr_axisMe(IAXIS)==0).or.(gr_axisMe(IAXIS)==gr_axisNumProcs(IAXIS)-1)) then
        count=1
        listOfBlocks(count)=1
     end if
  case (JBDRY_BLKS)
     if((gr_axisMe(JAXIS)==0).or.(gr_axisMe(JAXIS)==gr_axisNumProcs(JAXIS)-1)) then
        count=1
        listOfBlocks(count)=1
     end if
  case(KBDRY_BLKS)
     if((gr_axisMe(KAXIS)==0).or.(gr_axisMe(KAXIS)==gr_axisNumProcs(KAXIS)-1)) then
        count=1
        listOfBlocks(count)=1
     end if
  case(ANY_BDRY_BLKS)
     do i = 1,NDIM
        if((gr_axisMe(i)==0).or.(gr_axisMe(i)==gr_axisNumProcs(i)-1)) then
           count=1
           listOfBlocks(count)=1
        end if
     end do
  case(LEAF)
     count = 1
     listOfBlocks(count)=1
  case(ALL_BLKS)
     count = 1
     listOfBlocks(count)=1
  case(PARENT_BLK)
     count = 0
     listOfBlocks(1) = 0
     write(*,*)'WARNING! Grid__GetListOfBlocks(PARENT_BLK,...) is meaningless in Uniform Grid, returning zero'
     call Logfile_stampMessage(&
          'WARNING! Grid__GetListOfBlocks(PARENT_BLK,...) is meaningless in Uniform Grid, returning zero')
  case(ANCESTOR)
     count = 0
     listOfBlocks(1) = 0
     write(*,*)'WARNING! Grid__GetListOfBlocks(ANCESTOR,...) is meaningless in Uniform Grid, returning zero'
     call Logfile_stampMessage("WARNING! UG Grid_getListOfBlocks with ANCESTOR invalid")

  case(TRAVERSED)
     if(gr_blockType==TRAVERSED) then
        count=1
     else
        count=0
     end if
     listOfBlocks(1)=count

  case(TRAVERSED_AND_ACTIVE)
     if(gr_blockType==TRAVERSED_AND_ACTIVE) then
        count=1
     else
        count=0
     end if
     listOfBlocks(1)=count

  case DEFAULT
     count = 1
     listOfBlocks(count)=1
     write(*,*)'WARNING! Grid__GetListOfBlocks(,...) is using a meaningless first argument'
     call Logfile_stampMessage('WARNING! Grid__GetListOfBlocks was called with a meaningless first argument')

  end select
  return
end subroutine Grid_getListOfBlocks
