!!****if* source/Grid/GridMain/Grid_genReleaseBlkPtr
!!
!! NAME
!!
!!  Grid_genReleaseBlkPtr
!!
!! SYNOPSIS
!!
!!  call Grid_genReleaseBlkPtr(integer(in) :: blockid,
!!                             real, dimension(:,:,:,:), pointer  :: dataptr,
!!                             integer(in),OPTIONAL  :: vardesc,
!!                             real, dimension(:,:,:,:), pointer,OPTIONAL  :: dataptr2,
!!                             real, dimension(:,:,:,:), pointer,OPTIONAL  :: dataptr3)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   blockid : ID of block in current processor
!!
!!   dataptr : data pointer one 
!!
!!   vardesc : varDesc array
!!
!!   dataptr2 : data pointer two
!!
!!   dataptr3 : data pointer three
!!
!!
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_genReleaseBlkPtr(blockID,dataPtr, varDesc,dataPtr2,dataPtr3)

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : Grid_releaseBlkPtr, Grid_ascReleaseBlkPtr

  implicit none
  integer, intent(in) :: blockID
  real, dimension(:,:,:,:), pointer :: dataPtr
  real, dimension(:,:,:,:), pointer,OPTIONAL :: dataPtr2, dataPtr3
  integer,intent(in),OPTIONAL :: varDesc(:)

  integer :: myVarDesc(VARDESC_SIZE)
  integer :: gds, dur

  if (present(varDesc)) then
     if (size(varDesc) .LT. VARDESC_VAR) then
        call Driver_abortFlash('Grid_genReleaseBlkPtr: varDesc argument has invalid length!')
     else
        myVarDesc(1:min(VARDESC_SIZE,size(varDesc))) = varDesc(1:min(VARDESC_SIZE,size(varDesc)))
     end if
     if (size(varDesc) .LT. VARDESC_NUM)      myVarDesc(VARDESC_NUM) = 1
     if (size(varDesc) .LT. VARDESC_GDS)      myVarDesc(VARDESC_GDS) = CENTER
     if (size(varDesc) .LT. VARDESC_DURATION) myVarDesc(VARDESC_DURATION) = VD_DUR_PERM

     gds = myVarDesc(VARDESC_GDS)
     dur = myVarDesc(VARDESC_DURATION)
     select case (gds)
     case(FACES)
        if (dur == VD_DUR_PERM) then
           call Grid_releaseBlkPtr(blockID, dataPtr,  FACEX)
           if (NDIM>1) call Grid_releaseBlkPtr(blockID, dataPtr2, FACEY)
           if (NDIM>2) call Grid_releaseBlkPtr(blockID, dataPtr3, FACEZ)
        else if (dur == VD_DUR_GASC) then
           call Grid_ascReleaseBlkPtr(blockID, dataPtr,  FACEX)
           if (NDIM>1) call Grid_ascReleaseBlkPtr(blockID, dataPtr2, FACEY)
           if (NDIM>2) call Grid_ascReleaseBlkPtr(blockID, dataPtr3, FACEZ)
        else
           call Driver_abortFlash('Grid_genReleaseBlkPtr: called with a "duration" that is not implemented!')
        end if
     case default
        if (dur == VD_DUR_PERM) then
           call Grid_releaseBlkPtr(blockID, dataPtr, gds)
        else if (dur == VD_DUR_GASC) then
           call Grid_ascReleaseBlkPtr(blockID, dataPtr,  gds)
        else
           call Driver_abortFlash('Grid_genReleaseBlkPtr: called with a "duration" that is not implemented!')
        end if
     end select
  end if


end subroutine Grid_genReleaseBlkPtr
