!!****if* source/physics/TreeRay/TreeRayMain/OpticalDepth/tr_odGetNodeStruct
!!
!! NAME
!!
!!  tr_odGetNodeStruct
!!
!!
!! SYNOPSIS
!!
!!   call tr_odGetNodeStruct(botnodesize, nodesize)
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

subroutine tr_odGetNodeStruct(nsize, bnsize)
  use tr_odData, ONLY : tr_odIBH2, tr_odIH2, tr_odIBCO, tr_odICO
  implicit none
#include "Flash.h"
  integer, intent(INOUT) :: bnsize, nsize

#if CHEMISTRYNETWORK == 4
  tr_odIH2  = nsize  + 1
  tr_odIBH2 = bnsize + 1
  nsize  = nsize  + 1
  bnsize = bnsize + 1
#endif
#if CHEMISTRYNETWORK == 5
  tr_odIH2  = nsize  + 1
  tr_odIBH2 = bnsize + 1
  tr_odICO  = nsize  + 2
  tr_odIBCO = bnsize + 2
  nsize  = nsize  + 2
  bnsize = bnsize + 2
#endif

  return
end subroutine tr_odGetNodeStruct

