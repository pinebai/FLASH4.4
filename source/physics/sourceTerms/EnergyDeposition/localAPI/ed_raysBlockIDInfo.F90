!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_raysBlockIDInfo
!!
!! NAME
!!
!!  ed_raysBlockIDInfo
!!
!! SYNOPSIS
!!
!!  call ed_raysBlockIDInfo ()
!!
!! DESCRIPTION
!!
!!  Extracts from the rays array the following information about the block ID's:
!!
!!                  1) the # of different block ID's
!!                  2) all unique block ID's
!!                  3) the # of rays for each block ID
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

subroutine ed_raysBlockIDInfo ()

  implicit none

  return
end subroutine ed_raysBlockIDInfo
