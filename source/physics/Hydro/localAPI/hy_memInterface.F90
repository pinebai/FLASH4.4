!!****ih* source/physics/Hydro/localAPI/hy_memInterface
!!
!! NAME
!!  hy_memInterface
!!
!! SYNOPSIS
!!  use hy_memInterface
!!
!! DESCRIPTION
!!  This is a local interface for memory management routines.
!!
!!***


Module hy_memInterface

  use hy_memScratchData, ONLY: hy_memScratchInitialize, hy_memAllocScratch, hy_memDeallocScratch

  implicit none


  interface hy_memGetBlkPtr
     subroutine hy_memGetBlkPtr(blockID,dataPtr, gridDataStruct)
       implicit none
       integer,intent(in) :: blockId
       real, pointer :: dataPtr(:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memGetBlkPtr
     subroutine hy_memGetBlk5Ptr(blockID,data5Ptr, gridDataStruct)
       implicit none
       integer,intent(in) :: blockId
       real, pointer :: data5Ptr(:,:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memGetBlk5Ptr
  end interface

#include "FortranLangFeatures.fh"

  interface hy_memReleaseBlkPtr
     subroutine hy_memReleaseBlkPtr(blockId, dataPtr, gridDataStruct)
       implicit none
       integer,intent(in) :: blockId
       real, pointer :: dataPtr(:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memReleaseBlkPtr
     subroutine hy_memReleaseBlk5Ptr(blockId, data5Ptr, gridDataStruct)
       implicit none
       integer,intent(in) :: blockId
       real, POINTER_INTENT_OUT :: data5Ptr(:,:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memReleaseBlk5Ptr
  end interface

End Module hy_memInterface
