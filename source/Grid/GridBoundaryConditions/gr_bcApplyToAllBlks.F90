!!****if* source/Grid/GridBoundaryConditions/gr_bcApplyToAllBlks
!!
!! NAME
!!  gr_bcApplyToAllBlks
!!
!! SYNOPSIS
!!
!!  gr_bcApplyToAllBlks(integer(IN) :: axis,
!!                      logical(IN) :: isWork)
!!  
!! DESCRIPTION 
!!
!!  This routine is a wrapper around gr_bcApplyToOneFace, and is used by UG and PM2.
!!  It calls gr_bcApplyToOneFace one each for lowerface and upperface, and repeats
!!  the process for all blocks in the grid.
!!  
!! 
!! ARGUMENTS
!!  
!!    axis           - the direction for applying BC, one of IAXIS, JAXIS, or KAXIS
!!    isWork         - is always false for UG. In PM2, if true, indicates that
!!                     the boundary conditions should be applied to work array
!!
!! NOTES
!!  A specific direction is required in axis - no ALLDIR at this point.
!!
!!***

subroutine gr_bcApplyToAllBlks(axis,isWork)
  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_getBlkBC, &
       Grid_getBlkIndexLimits
  use gr_bcInterface, ONLY : gr_bcApplyToOneFace
  use Grid_data,ONLY : gr_numDataStruct,gr_gridDataStruct,gr_gridDataStructSize
  implicit none
#include "constants.h"
#include "Flash.h"
  
  integer, intent(in) :: axis
  logical, intent(in) :: isWork

  logical, dimension(LOW:HIGH) :: loop
  integer, dimension(LOW:HIGH) :: bcType,face
  integer,dimension(LOW:HIGH,MDIM) :: blkLimitsGC,blkLimits,blkBC
  integer :: varCount,gridDataStruct
  integer :: edge,struct
  integer :: blkNum,blockID
  integer,dimension(MDIM) :: regionType
  integer :: idest=0
  integer,dimension(gr_numDataStruct) :: localDataStruct,localStructSize
  integer :: localNum

  if(isWork) then
     localNum=1
     localDataStruct(localNum)=WORK
     localStructSize(localNum)=1
  else
     localNum=gr_numDataStruct
     localDataStruct=gr_gridDataStruct(1:localNum)
     localStructSize=gr_gridDataStructSize(1:localNum)
  end if

  if (localNum > 0) call Grid_getLocalNumBlks(blkNum)   !! Find the local number of blocks

  do struct=1,localNum
     gridDataStruct=localDataStruct(struct)
     varCount=localStructSize(struct)
     
     face(LOW)=LOW
     face(HIGH)=HIGH
     
     do blockID=1,blkNum                 !! and loop over them
        
        call Grid_getBlkBC(blockID,blkBC)  !! For the block find the faces on the boundary
        loop(LOW)=blkBC(LOW,axis)/=NOT_BOUNDARY      !! along the specified dimension
        loop(HIGH)=blkBC(HIGH,axis)/=NOT_BOUNDARY
        
        if(loop(LOW).or.loop(HIGH)) then           !! if at least one face is on the boundary
           bcType=blkBC(:,axis)                    !! then proceed with calculation
           
           !! The next few statements are to get block information to 
           !! prepare for the calculation
           call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,&
                gridDataStruct)
           !! Now that all preparatory work is done, call the routine
           !! that will repackage relevant parts of the block data
           !! to pass on the boundary condition routines that do actual
           !! calculation.
           !! The boundary condition routines will see data arrays
           !! that may be equivalent to some transpose of the block
           !! data, depending on the direction given by axis.
           do edge = LOW,HIGH
              regionType(:)=WHOLE_VECTOR
              regionType(axis)=face(edge)
              if(loop(edge))&
                   call gr_bcApplyToOneFace(axis,bcType(edge),&
                   gridDataStruct,varCount,regionType,&
                   blkLimits,blkLimitsGC,blockID,idest)
           end do
        end if
     end do
  end do
end subroutine gr_bcApplyToAllBlks
