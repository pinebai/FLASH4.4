!!****f* source/Grid/Grid_GCTransferOneBlk
!!
!! NAME
!!  Grid_GCTransferOneBlk
!!
!! SYNOPSIS
!!
!!  Grid_GCTransferOneBlk(logical(IN)   :: mode,
!!                        integer(IN)   :: gridDataStruct
!!                        integer(IN)   :: blkIndex
!!                        real, pointer :: blkData(:,:,:,:))
!!  
!!  
!! DESCRIPTION
!!  
!!  This routine can transfers data from a flat array (used for storing saved
!!  guard cell values) to a block array (which will have a full block structure
!!  but invalid interior values; valid values only in the guard cells specified
!!  at the time of their storage), or it can extract data from the mesh arrays one
!!  at a time and store it in the flat array. Normally users are advised to use the
!!  Grid_GCputScratch interface for transferring the data to the flat array, and use 
!!  the current interface for fetching the data from the flat array.
!!
!! ARGUMENTS
!!            
!!  mode           : indicates the direction of data tranfer. 
!!                   For mode=.true. the data is transferred from blkData to
!!                   a flat array allocated by a call to gr_GCAllocScratch
!!                   For mode=.false. the data is transferred to blkData from 
!!                   the saved flat array 
!!  gridDataStruct : integer value specifying data structure. 
!!                   The options are defined in constants.h, the ones
!!                   relevant to this routine are :
!!                   CENTER cell centered variables (default)
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along IAXIS
!!                   made
!!  blkIndex       : index into the list of block, from where blkID can
!!                   be retrieved for the target block
!!  blkData        : pointer to a real four dimensional array. If 
!!
!!  NOTES
!!   variables that start with "gr_" are variables of Grid unit scope
!!   and are stored in the fortran module Grid_data. Variables are not
!!   starting with gr_ are local variables or arguments passed to the 
!!   routine.
!!
!!  SEE ALSO
!!    Grid_GCPutScratch 
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_GCTransferOneBlk(mode,gridDataStruct,blkIndex,blkData)

  implicit none

  logical, intent(IN) :: mode
  integer, intent(IN) :: gridDataStruct, blkIndex

  real, pointer, dimension(:,:,:,:) :: blkData

end subroutine Grid_GCTransferOneBlk
