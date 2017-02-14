!!****if* source/diagnostics/ProtonImaging/localAPI/pi_protonsBlockIDInfo
!!
!! NAME
!!
!!  pi_protonsBlockIDInfo
!!
!! SYNOPSIS
!!
!!  call pi_protonsBlockIDInfo ()
!!
!! DESCRIPTION
!!
!!  Extracts from the protons array the following information about the block ID's:
!!
!!                  1) the # of different block ID's
!!                  2) all unique block ID's
!!                  3) the # of protons for each block ID
!!
!!  This information is extracted from a given linear array of block ID's without
!!  reordering. The algorithm simply loops over the entire array of block ID's and
!!  records any changes in these numbers. This means that even if a pair of block ID's
!!  has identical value but each occurs in separate places of the array, the code
!!  registers this as two separate block ID's. Hence the minimum number of different
!!  block ID's is found only when the block ID array is completely ordered.
!!  Unordered block ID arrays do not lead to a fault, but only diminish performance.
!!
!! NOTES
!!
!!***

subroutine pi_protonsBlockIDInfo ()

  implicit none

  return
end subroutine pi_protonsBlockIDInfo
