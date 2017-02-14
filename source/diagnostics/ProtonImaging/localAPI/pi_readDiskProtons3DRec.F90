!!****if* source/diagnostics/ProtonImaging/localAPI/pi_readDiskProtons3DRec
!!
!! NAME
!!
!!  pi_readDiskProtons3DRec
!!
!! SYNOPSIS
!!
!!  call pi_readDiskProtons3DRec (integer, intent (in)    :: blockCount,
!!                                integer, intent (in)    :: blockList (:),
!!                                logical, intent (inout) :: moreOnDisk)
!!
!! DESCRIPTION
!!
!!  Reads in a batch of disk protons from the old disk proton file for 3D rectangular
!!  (cartesian) geometries and determines the list of protons that are present in one
!!  of the blocks on the current processor. Their block ID's are not ordered.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  moreOnDisk     : if true, there are more disk protons on the old disk proton file
!!
!! NOTES
!!
!!***

subroutine pi_readDiskProtons3DRec (blockCount, blockList, moreOnDisk)

  implicit none

  integer, intent (in)    :: blockCount
  integer, intent (in)    :: blockList (1:blockCount)
  logical, intent (inout) :: moreOnDisk

  moreOnDisk = .false.

  return
end subroutine pi_readDiskProtons3DRec
