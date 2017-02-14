!!****if* source/Grid/GridMain/Grid_releaseBlkPtr
!!
!! NAME
!!  Grid_releaseBlkPtr
!!
!! SYNOPSIS
!!
!!  Grid_releaseBlkPtr(integer(IN)   :: blockId,
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

subroutine Grid_releaseBlkPtr(blockId, dataPtr, gridDataStruct)

#include "Flash.h"
#ifdef Grid_releaseBlkPtr
! disabling drift macro expansion, see: drift
#undef Grid_releaseBlkPtr
#endif


#include "constants.h"

  use Driver_interface, only: Driver_driftBlock
  use Driver_data, only: dr_driftSrcFile, dr_driftSrcLine
#ifdef FL_NON_PERMANENT_GUARDCELLS
  use Grid_data, ONLY: gr_blkPtrRefCount, gr_blkPtrRefCount_fc,&
                       gr_ccMask, gr_fcMask
  use physicaldata, ONLY : gcell_on_cc, gcell_on_fc
#endif

  implicit none

  integer,intent(in) :: blockId
  real, pointer :: dataPtr(:,:,:,:)
  integer,optional, intent(in) :: gridDataStruct

  integer :: gds
#ifdef FL_NON_PERMANENT_GUARDCELLS
  integer :: blkPtrRefCount
  integer :: idest, iopt, nlayers, icoord
  logical :: lcc, lfc, lec, lnc, l_srl_only, ldiag
  logical, dimension(NUNK_VARS) :: save_ccMask
#if NFACE_VARS>0
  logical, dimension(3,NFACE_VARS) :: save_fcMask
#endif
#endif

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

#ifdef FL_NON_PERMANENT_GUARDCELLS
  if (gds .eq. CENTER .or. gds .eq. FACEX .or. gds .eq. FACEY .or. gds .eq. FACEZ) then
     idest = 1
     iopt = 1
     if (gds .eq. FACEX .or. gds .eq. FACEY .or. gds .eq. FACEZ) then
        blkPtrRefCount = gr_blkPtrRefCount_fc
#if NFACE_VARS>0
        save_fcMask=gcell_on_fc
        gcell_on_fc=gr_fcMask
#endif
        lcc = .false.
        lfc = .true.
     else
        blkPtrRefCount = gr_blkPtrRefCount
        save_ccMask=gcell_on_cc
        gcell_on_cc=gr_ccMask
        lcc = .true.
        lfc = .false.
     end if
     lec = .false.
     lnc = .false.

     blkPtrRefCount = blkPtrRefCount - 1

     !only put this pointer back in unk when everyone is done with it
     if (blkPtrRefCount .lt. 1) then 
        call amr_1blk_to_perm( lcc,lfc,lec,lnc,blockId,iopt,idest)
     end if

     if (gds .eq. FACEX .or. gds .eq. FACEY .or. gds .eq. FACEZ) then
        gr_blkPtrRefCount_fc = blkPtrRefCount
#if NFACE_VARS>0
        gcell_on_fc=save_fcMask
#endif
     else if (gds .eq. CENTER) then
        gr_blkPtrRefCount = blkPtrRefCount
        gcell_on_cc=save_ccMask
     end if

  end if
#endif !  FL_NON_PERMANENT_GUARDCELLS

  ! always destroy the pointer, because the other users will have their 
  ! own
  nullify(dataPtr)

end subroutine Grid_releaseBlkPtr








