!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_bhGetNodeStruct
!!
!! NAME
!!
!!  TreeRay_bhGetNodeStruct
!!
!!
!! SYNOPSIS
!!
!!   call TreeRay_bhGetNodeStruct(botnodesize, nodesize)
!!
!! DESCRIPTION
!!
!!   Calculates structure (size, indeces of variables) of the 
!!   treeray section of the tree node
!!
!! ARGUMENTS
!!
!!
!! RESULT
!!
!!
!!***

subroutine TreeRay_bhGetNodeStruct(im, ix, iy, iz, nsize, bnsize)
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use TreeRay_data, ONLY : tr_bhIM, tr_bhIBM, tr_bhIX, tr_bhIY, tr_bhIZ, &
    tr_meshMe, tr_bhUseTreeRay
  use tr_odInterface, ONLY : tr_odGetNodeStruct
  !use tr_osInterface, ONLY : tr_osGetNodeStruct
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(IN) :: im, ix, iy, iz
  integer, intent(INOUT) :: bnsize, nsize
  
  ! this subroutine is called before TreeRay_init and therefore
  ! it has to determine useTreeRay and tr_meshMe by itself
  call RuntimeParameters_get ("useTreeRay", tr_bhUseTreeRay)
  call Driver_getMype(MESH_COMM, tr_meshMe)

  if (.not. tr_bhUseTreeRay) return

  ! set indeces of grid section quantities
  tr_bhIM  = im
  tr_bhIBM = im
  tr_bhIX  = ix
  tr_bhIY  = iy
  tr_bhIZ  = iz


  ! update node sizes and set indices in subroutine 
  if (tr_meshMe == MASTER_PE) &
    print*,  "[TreeRay_GetNodeStruct:] bnsize before=", bnsize 
  call tr_odGetNodeStruct(nsize, bnsize)
  !call tr_osGetNodeStruct(nsize, bnsize)
  if (tr_meshMe == MASTER_PE) &
    print*,  "[TreeRay_GetNodeStruct:] bnsize after=", bnsize 

  return
end subroutine TreeRay_bhGetNodeStruct
