!!****if* source/Grid/GridMain/gr_setDataStructInfo
!!
!! NAME
!!  gr_setDataStructInfo
!!
!! SYNOPSIS
!!
!!  gr_setDataStructInfo()
!!  
!! DESCRIPTION 
!!  
!!  This routine consolidates the information about various grid data structures
!!  into an array format for use by the routines that handle the boundary
!!  conditions. With the array formats, the boundary conditions can be 
!!  handled in a loop, instead of laying out code explicitly for each
!!  data structure. Also makes it expandable
!!  
!!  
!!
!! ARGUMENTS 
!!
!!
!!***


subroutine gr_setDataStructInfo()
#include "constants.h"
#include "Flash.h"
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_data, ONLY : gr_numDataStruct,gr_gridDataStruct,gr_gridDataStructSize, &
                        gr_bcEnableApplyMixedGds

  implicit none 

  gr_numDataStruct=1
  gr_gridDataStruct(gr_numDataStruct)=CENTER
  gr_gridDataStructSize(gr_numDataStruct)=NUNK_VARS
  if(NFACE_VARS>0) then
     gr_numDataStruct=gr_numDataStruct+1
     gr_gridDataStruct(gr_numDataStruct)=FACEX
     gr_gridDataStructSize(gr_numDataStruct)=NFACE_VARS
     if(NDIM>1) then
        gr_numDataStruct=gr_numDataStruct+1
        gr_gridDataStruct(gr_numDataStruct)=FACEY
        gr_gridDataStructSize(gr_numDataStruct)=NFACE_VARS
     end if
     if(NDIM>2) then
        gr_numDataStruct=gr_numDataStruct+1
        gr_gridDataStruct(gr_numDataStruct)=FACEZ
        gr_gridDataStructSize(gr_numDataStruct)=NFACE_VARS
     end if
     if (gr_numDataStruct < size(gr_gridDataStruct)) then
        call RuntimeParameters_get("gr_bcEnableApplyMixedGds", gr_bcEnableApplyMixedGds)
        if (gr_bcEnableApplyMixedGds) then
           gr_numDataStruct=gr_numDataStruct+1
           gr_gridDataStruct(gr_numDataStruct)=CENTER_FACES
           gr_gridDataStructSize(gr_numDataStruct)=sum(gr_gridDataStructSize(1:gr_numDataStruct-1))
        end if
     end if
  end if
  return
end subroutine gr_setDataStructInfo
