!!****if* source/Grid/GridMain/gr_checkDataType
!!
!! NAME
!!  gr_checkDataType
!!
!! SYNOPSIS
!!
!!  gr_checkDataType (integer(IN) :: gridDataStruct,
!!                    integer(OUT) :: imax,
!!                    integer(OUT) :: jmax,
!!                    integer(OUT) :: kmax,
!!                    logical(IN)  :: isget)
!!  
!! DESCRIPTION 
!!  
!!  This is a helper routine for all the Grid_get/put routines. It tests for
!!  validity of the gridDataStruct arguments (which indicates the data structure 
!!  of the grid package from which to fetch data), 
!!  All the of grid data structures are common to all get and put
!!  routines, however, the get routines also permit cell volume and 
!!  cell face areas as valid quantities to be fetched.
!!  It also calculates the outer bounds of the blocksize for 
!!  diagnostic purposes 
!!
!! ARGUMENTS 
!!
!!
!!  gridDataStruct : integer value specifying the type of data desired.
!!             Valid options are either one of the Grid data structures,
!!             or one of the derived quantities such as the cell volume
!!             or the cell area.
!!
!!             The options are defined in constants.h and they are :
!!                   CENTER cell centered variables
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along IAXIS
!!                   WORK   single, cell centered variable, valid only
!!                          for paramesh
!!                   SCRATCH space for saving previous timestep data if needed
!!              The next two options are valid only if isget=.true.
!!                   CELL_VOLUME volumes of specified cells 
!!                   CELL_FACEAREA face area for the specified cells on 
!!                                 specified face 
!!
!!  imax -     upper bound along IAXIS
!!
!!  jmax -     upper bound along JAXIS
!!
!!  kmax -     upper bound along KAXIS
!!
!!  isget -   if true indicates that this routine is called by a 
!!            get function therefore cell volume and face area are
!!            valid values for gridDataStruct. If isget is false then it is
!!            called by a put routine and those values are invalid.
!!
!!***

subroutine gr_checkDataType(blockID,gridDataStruct,imax,jmax,kmax,isget)
  use Grid_interface, ONLY : Grid_getBlkIndexLimits
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
#include "constants.h"
  integer,intent(IN) :: blockID,gridDataStruct
  integer,intent(OUT) :: imax,jmax,kmax
  logical,intent(IN) :: isget
  logical :: validDataType
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC


  !do some error checking

  !if we are using FACE variables, the cell range is one higher
  !than the number of zones in the block.  set a max value so
  !we can do error checking


  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,gridDataStruct)
  imax = blkLimitsGC(HIGH,IAXIS)
  jmax = blkLimitsGC(HIGH,JAXIS)
  kmax = blkLimitsGC(HIGH,KAXIS)
  
  validDataType = .false.
  if((gridDataStruct == FACEX).or.(gridDataStruct == SCRATCH_FACEX)) then
     imax = imax + 1
     validDataType=.true.
     
  else if((gridDataStruct == FACEY).or.(gridDataStruct == SCRATCH_FACEY)) then
     jmax = jmax + 1
     validDataType=.true.
     
  else if((gridDataStruct == FACEZ).or.(gridDataStruct == SCRATCH_FACEZ)) then
     kmax = kmax + 1
     validDataType=.true.
     
  else if((gridDataStruct == CENTER).or.(gridDataStruct == SCRATCH).or.(gridDataStruct == SCRATCH_CTR)) then
     validDataType=.true.
     
  elseif(isget) then
     if((gridDataStruct==CELL_VOLUME).or.(gridDataStruct == CELL_FACEAREA)) then
        validDataType=.true.
     end if
     
#ifdef FLASH_GRID_PARAMESH
  else if(gridDataStruct == WORK) then
     validDataType=.true.
#endif
     
  end if
  
  if(.not.validDataType) then
     print *, "Grid_checkDataType: gridDataStruct set to improper value"
     if(isget) then
        print *, "gridDataStruct has an undefine value (valid values defined in constants.h)"
        print*,'it actually is ',gridDataStruct
     else
        print *, "gridDataStruct has an undefine value (valid values defined in constants.h)"
     end if
     call Driver_abortFlash(" invalid gridDataStruct in get/put functions ")
  end if
  return
end subroutine gr_checkDataType
