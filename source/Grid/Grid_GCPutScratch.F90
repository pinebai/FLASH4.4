!!****f* source/Grid/Grid_GCPutScratch
!!
!! NAME
!!  Grid_GCPutScratch
!!
!! SYNOPSIS
!!
!!  Grid_GCPutScratch(              integer(IN) :: gridDataStruct,
!!                                  logical(IN) :: needSetup,
!!                                  logical(IN) :: releaseSetup,
!!                         optional,integer(IN) :: blkCount,
!!                         optional,integer(IN) :: blkList(:),
!!                         optional,integer(IN) :: indCount,
!!                         optional,integer(IN) :: indList(:),
!!                         optional,integer(IN) :: gcCnt(:),
!!                         optional,logical(IN) :: putData)
!!  
!!  
!! DESCRIPTION
!!  
!!  This routine stores guardcells of  the specified indices of 
!!  grid data structures from a list of blocks. If putData is .false.
!!  then the routine is said to be called in initialize only mode.
!!  In this mode "needSetUp" must be true. If putData is missing, it is 
!!  assumed to be true. 
!!  If the argument needSetup is true, then all the optional 
!!  arguments must be included in the call, and the routine will
!!  make a call to gr_allocGCScratch function. The first call to 
!!  Grid_GCPutScratch in a simulation must be made with needSetup=.true.
!!  for correct operation. 
!!  The subsequent calls to Grid_GCPutScratch can be made without the optional
!!  arguments by setting needSetup to false, as long the values of the
!!  optional arguments provided in the call with needSetup=true
!!  are unchanged in the simulation.
!!
!!  If there is any change in the grid, list of blocks, or the number
!!  or order of indices, then both needSetup and releaseSetup must be set
!!  to true, and appropriate values must be provided in the optional 
!!  argument, which must all be present.
!!
!! ARGUMENTS
!!            
!!  gridDataStruct : integer value specifying data structure. 
!!                   The options are defined in constants.h, the ones
!!                   relevant to this routine are :
!!                   CENTER cell centered variables (default)
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along IAXIS
!!  needSetup      : if true a call to gr_allocateGuardScratchSpace will be
!!                   made
!!  releaseSetup   : if true a call to gr_releaseGuardScratchSpace will be
!!                   made before the call to gr_allocateGuardScratchSpace
!!  blkCount       : count of blocks whose guard cells are to be saved
!!  blkList        : list of blocks whose guard cells are to be saved
!!  indCount     : count of indices in the data structure to be saved
!!  indList      : list of indices in the data structure to be saved
!!  gcCnt         : the count of guardcells along each dimension to be saved
!!  putData      : if .false. then the routine is said to be called in 
!!                 initialize only mode.
!!
!!
!!  NOTES
!!   variables that start with "gr_" are variables of Grid unit scope
!!   and are stored in the fortran module Grid_data. Variables are not
!!   starting with gr_ are local variables or arguments passed to the 
!!   routine.
!!
!!  SEE ALSO
!!   
!!   gr_GCAllocScratch and gr_GCTransferOneBlk
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_GCPutScratch(gridDataStruct,needSetup,releaseSetup,&
     &blkCount,blkList,indCount,indList, gcCnt, putData)

  implicit none

  integer, intent(IN) :: gridDataStruct 
  logical, intent(IN) :: needSetup, releaseSetup
  integer, optional, intent(IN) :: blkCount, indCount
  integer, optional, dimension(:), intent(IN) :: blkList
  integer, optional, dimension(:), intent(IN) :: indList
  integer, optional, dimension(NDIM),intent(IN) :: gcCnt  
  logical, optional, intent(IN) :: putData
  
end subroutine Grid_GCPutScratch
