!!****f* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistry
!!
!! NAME
!!
!!  PrimordialChemistry
!!
!! SYNOPSIS
!!
!!  call PrimordialChemistry(integer(IN) :: blockCount
!!                           integer(IN) :: blockList(blockCount),
!!                              real(IN) :: dt)
!!
!!
!!
!! DESCRIPTION
!!  Apply chemistry operation on the list of blocks provided as input
!!
!! ARGUMENTS
!!
!!  blockCount : The number of blocks in the list
!!  blockList(:) : The list of blocks on which to apply the operation
!!  dt : the current timestep
!!
!!
!!
!!***


subroutine PrimordialChemistry (blockCount, blockList, dt)

  implicit none
  
  integer, intent(IN) :: blockCount
  integer, dimension(blockCount), intent(IN) :: blockList
  real, intent(IN) :: dt

return
end subroutine PrimordialChemistry
