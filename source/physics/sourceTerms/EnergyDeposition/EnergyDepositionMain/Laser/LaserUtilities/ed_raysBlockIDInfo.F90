!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_raysBlockIDInfo
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

  use EnergyDeposition_data,    ONLY : ed_rayBlockID,       &
                                       ed_rayBlockIDCount,  &
                                       ed_rayCount,         &
                                       ed_rayNumberBlockID, &
                                       ed_rays

  implicit none

#include "EnergyDeposition.h"

  integer :: blockID
  integer :: blockIDprev
  integer :: number
  integer :: ray
!
!
!     ...Analyze the block ID info of the rays array and store for future reference:
!        1) the # of different block ID's, 2) all unique block ID's and 3) the # of rays
!        for each block ID. This is done in a two-stage pass through the rays array in
!        order to avoid fixed size declarations of the new block ID info arrays involved.
!
!        If no rays are present, do nothing.
!
!           
  if (ed_rayCount == 0) then
      ed_rayBlockIDCount = 0
  else

      if (allocated  (ed_rayNumberBlockID)) then
          deallocate (ed_rayNumberBlockID)
      end if

      if (allocated  (ed_rayBlockID)) then
          deallocate (ed_rayBlockID)
      end if
!
!
!     ...1st stage: determine # of different block ID's.
!
!           
      ed_rayBlockIDCount = 1
      blockIDprev = int (ed_rays (RAY_BLCK,1))

      do ray = 1,ed_rayCount

         blockID = int (ed_rays (RAY_BLCK,ray))

         if (blockID /= blockIDprev) then
             ed_rayBlockIDCount = ed_rayBlockIDCount + 1
             blockIDprev = blockID
         end if
      end do

      allocate (ed_rayNumberBlockID (1:ed_rayBlockIDCount))
      allocate (ed_rayBlockID       (1:ed_rayBlockIDCount))
!
!
!     ...2nd stage: store all unique block ID's and the # of rays for each block ID.
!
!           
      number = 0
      ed_rayBlockIDCount = 1
      blockIDprev = int (ed_rays (RAY_BLCK,1))

      do ray = 1,ed_rayCount

         blockID = int (ed_rays (RAY_BLCK,ray))

         if (blockID /= blockIDprev) then
             ed_rayNumberBlockID (ed_rayBlockIDCount) = number
             ed_rayBlockID       (ed_rayBlockIDCount) = blockIDprev
             number = 1
             ed_rayBlockIDCount = ed_rayBlockIDCount + 1
             blockIDprev = blockID
         else
             number = number + 1
         end if
     end do

     ed_rayNumberBlockID (ed_rayBlockIDCount) = number
     ed_rayBlockID       (ed_rayBlockIDCount) = blockIDprev

  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ed_raysBlockIDInfo
