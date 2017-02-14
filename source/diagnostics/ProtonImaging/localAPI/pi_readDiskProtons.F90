!!****if* source/diagnostics/ProtonImaging/localAPI/pi_readDiskProtons
!!
!! NAME
!!
!!  pi_readDiskProtons
!!
!! SYNOPSIS
!!
!!  call pi_readDiskProtons (integer, intent (in)    :: blockCount,
!!                           integer, intent (in)    :: blockList (:),
!!                           logical, intent (inout) :: moreOnDisk)
!!
!! DESCRIPTION
!!
!!  Reads in a batch of disk protons from the old disk proton file. This routine
!!  calls the appropriate subroutines according to the domain grid geometry specified.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  moreOnDisk     : if true, there are more disk protons on the old disk proton file
!!
!!***

subroutine pi_readDiskProtons (blockCount, blockList, moreOnDisk)

  implicit none

  integer, intent (in)    :: blockCount
  integer, intent (in)    :: blockList (1:blockCount)
  logical, intent (inout) :: moreOnDisk

  moreOnDisk = .false.

  return
end subroutine pi_readDiskProtons
