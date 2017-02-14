!!****if* source/Grid/GridMain/Chombo/Grid_getListOfBlocks
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
!!  region_bndBox - when blocktype is specified as INREGION this argument defines the
!!                  bounding box of the region
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

  use Logfile_interface, ONLY : Logfile_stampMessage
  use Driver_interface, ONLY : Driver_abortFlash
  use chombo_f_c_interface, ONLY : ch_get_block_ids
  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockType
  integer,dimension(MAXBLOCKS),intent(out) :: listOfBlocks
  integer,intent(out) :: count
  integer,intent(IN),optional :: refinementLevel
  real, dimension(LOW:HIGH,MDIM), intent(IN), optional :: region_bndBox
  logical, intent(IN), optional :: includePartialBlocks
  integer :: level

  if (present(refinementLevel)) then
     if (blockType == REFINEMENT) then
        level = refinementLevel
     else
        call Driver_abortFlash &
             ("[Grid_getListOfBlocks]: blockType must be 'REFINEMENT'")
     end if
  else
     level = -1
  end if

  call ch_get_block_ids(blockType, level, listOfBlocks, count)

end subroutine Grid_getListOfBlocks
