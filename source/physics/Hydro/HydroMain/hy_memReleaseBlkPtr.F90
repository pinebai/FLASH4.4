!!****if* source/physics/Hydro/HydroMain/hy_memReleaseBlkPtr
!!
!! NAME
!!  hy_memReleaseBlkPtr
!!
!! SYNOPSIS
!!
!!  call hy_memReleaseBlkPtr(integer(IN)   :: blockId,
!!                     real(pointer) :: dataPtr(:,:,:,:),
!!                     integer(IN),optional :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  Releases a pointer to a block.
!!  
!! ARGUMENTS 
!!
!!  blockId - ID of the block, should be the same ID was used in the
!!            corresponding Grid_getBlkPtr call.
!!  dataPtr - Pointer to be released.
!!  gridDataStruct - an optional argument that designates the type of grid data
!!                   structure to handle (i.e. facevar, unknown, scratch...)
!!
!!
!! NOTES
!!
!!  This implementation actually does more than just releasing the pointer.
!!
!! SEE ALSO
!!  Grid_getBlkPtr
!!***

subroutine hy_memReleaseBlkPtr(blockId, dataPtr, gridDataStruct)

#include "Flash.h"
#ifdef hy_memReleaseBlkPtr
! disabling drift macro expansion, see: drift
#undef hy_memReleaseBlkPtr
#endif


#include "constants.h"

  use Driver_interface, only: Driver_driftBlock
  use Driver_data, only: dr_driftSrcFile, dr_driftSrcLine

  implicit none

  integer,intent(in) :: blockId
  real, pointer :: dataPtr(:,:,:,:)
  integer,optional, intent(in) :: gridDataStruct

  integer :: gds

  if(present(gridDataStruct)) then
     gds = gridDataStruct
  else
     gds = CENTER
  end if

#if DRIFT_ENABLE
  if(dr_driftSrcLine >= 0) then
    call Driver_driftBlock(dr_driftSrcFile, dr_driftSrcLine, blockId, &
      dataPtr(:,&
        lbound(dataPtr,2)+K1D*NGUARD : ubound(dataPtr,2)-K1D*NGUARD,&
        lbound(dataPtr,3)+K2D*NGUARD : ubound(dataPtr,3)-K2D*NGUARD,&
        lbound(dataPtr,4)+K3D*NGUARD : ubound(dataPtr,4)-K3D*NGUARD),&
      gds)
  end if
  dr_driftSrcLine = -1
#endif


  ! always destroy the pointer, because the other users will have their 
  ! own
  nullify(dataPtr)

end subroutine hy_memReleaseBlkPtr

subroutine hy_memReleaseBlk5Ptr(blockId, data5Ptr, gridDataStruct)

#include "FortranLangFeatures.fh"

  implicit none

  integer,intent(in) :: blockId
  real, POINTER_INTENT_OUT :: data5Ptr(:,:,:,:,:)
  integer,optional, intent(in) :: gridDataStruct

  integer :: gds

  if(present(gridDataStruct)) then
     gds = gridDataStruct
  else
     gds = CENTER
  end if

  nullify(data5Ptr)

end subroutine hy_memReleaseBlk5Ptr
