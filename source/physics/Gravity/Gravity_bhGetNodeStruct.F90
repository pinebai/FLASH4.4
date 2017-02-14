!!****f* source/physics/Gravity/Gravity_bhGetNodeStruct
!!
!! NAME
!!
!!  Gravity_bhGetNodeStruct
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhGetNodeStruct(
!!                  integer(in)    :: im,
!!                  integer(in)    :: ix,
!!                  integer(in)    :: iy,
!!                  integer(in)    :: iz,
!!                  integer(inout) :: nsize,
!!                  integer(inout) :: bnsize
!!        )
!!
!! DESCRIPTION
!!
!!   Calculates structure (size, indeces of variables) of the 
!!   gravity section of the tree node.
!!
!! ARGUMENTS
!!
!!  im          : index of the mass in the node array 
!!  ix          : index x coord of the mass centre in the node array 
!!  iy          : index y coord of the mass centre in the node array 
!!  iz          : index z coord of the mass centre in the node array 
!!  nsize       : size of the node array (initially, before the gravity part is
!!                added; at return, after the gravity part is added)
!!  bnsize      : size of the bottom-most node array (initially, before 
!!                the gravity part is added; at return, after the 
!!                gravity part is added)
!!
!!
!!***

subroutine Gravity_bhGetNodeStruct(im, ix, iy, iz, nsize, bnsize)
  implicit none
  integer, intent(IN) :: im, ix, iy, iz
  integer, intent(INOUT) :: nsize, bnsize
  
  return
end subroutine Gravity_bhGetNodeStruct
