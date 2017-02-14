!!****f* source/physics/TreeRay/TreeRay_bhGetNodeStruct
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
  implicit none
  integer, intent(IN) :: im, ix, iy, iz
  integer, intent(INOUT) :: bnsize, nsize
  
  return
end subroutine TreeRay_bhGetNodeStruct
