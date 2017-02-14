!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonUtilities/pi_protonsBlockIDInfo
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

  use ProtonImaging_data,    ONLY : pi_protonBlockID,       &
                                    pi_protonBlockIDCount,  &
                                    pi_protonCount,         &
                                    pi_protonNumberBlockID, &
                                    pi_protons

  implicit none

#include "ProtonImaging.h"

  integer :: blockID
  integer :: blockIDprev
  integer :: number
  integer :: proton
!
!
!     ...Analyze the block ID info of the protons array and store for future reference:
!        1) the # of different block ID's, 2) all unique block ID's and 3) the # of protons
!        for each block ID. This is done in a two-stage pass through the protons array in
!        order to avoid fixed size declarations of the new block ID info arrays involved.
!
!        If no protons are present, do nothing.
!
!           
  if (pi_protonCount == 0) then
      pi_protonBlockIDCount = 0
  else

      if (allocated  (pi_protonNumberBlockID)) then
          deallocate (pi_protonNumberBlockID)
      end if

      if (allocated  (pi_protonBlockID)) then
          deallocate (pi_protonBlockID)
      end if
!
!
!     ...1st stage: determine # of different block ID's.
!
!           
      pi_protonBlockIDCount = 1
      blockIDprev = int (pi_protons (PROTON_BLCK,1))

      do proton = 1,pi_protonCount

         blockID = int (pi_protons (PROTON_BLCK,proton))

         if (blockID /= blockIDprev) then
             pi_protonBlockIDCount = pi_protonBlockIDCount + 1
             blockIDprev = blockID
         end if
      end do

      allocate (pi_protonNumberBlockID (1:pi_protonBlockIDCount))
      allocate (pi_protonBlockID       (1:pi_protonBlockIDCount))
!
!
!     ...2nd stage: store all unique block ID's and the # of protons for each block ID.
!
!           
      number = 0
      pi_protonBlockIDCount = 1
      blockIDprev = int (pi_protons (PROTON_BLCK,1))

      do proton = 1,pi_protonCount

         blockID = int (pi_protons (PROTON_BLCK,proton))

         if (blockID /= blockIDprev) then
             pi_protonNumberBlockID (pi_protonBlockIDCount) = number
             pi_protonBlockID       (pi_protonBlockIDCount) = blockIDprev
             number = 1
             pi_protonBlockIDCount = pi_protonBlockIDCount + 1
             blockIDprev = blockID
         else
             number = number + 1
         end if
     end do

     pi_protonNumberBlockID (pi_protonBlockIDCount) = number
     pi_protonBlockID       (pi_protonBlockIDCount) = blockIDprev

  end if
!
!
!     ...Ready!
!
!
  return
end subroutine pi_protonsBlockIDInfo
