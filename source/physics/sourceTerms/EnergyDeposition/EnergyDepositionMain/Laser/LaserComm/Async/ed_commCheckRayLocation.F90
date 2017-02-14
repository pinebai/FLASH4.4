!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Async/ed_commCheckRayLocation
!!
!!  NAME     
!!   ed_commCheckRayLocation
!!
!!  SYNOPSIS
!!   ed_commCheckRayLocation(integer, intent(IN) :: fromRay,
!!                           integer, intent(IN) :: toRay)
!!
!!  DESCRIPTION 
!!    This subroutine checks that the ray position is within the bounding
!!    box of the correct block.
!!    
!!  ARGUMENTS
!!    fromRay : first ray in ed_rays to check
!!    toRay   :  last ray in ed_rays to check
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "EnergyDeposition.h"

subroutine ed_commCheckRayLocation(fromRay, toRay)
  use ed_commData, ONLY : ed_commRayPosIndex
  use EnergyDeposition_data, ONLY : ed_meshMe, ed_rays
  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_getBlkBoundBox
  implicit none
  integer, intent(IN) :: fromRay, toRay
  real, dimension(LOW:HIGH, MDIM) :: bndBox
  integer :: numBlocks, i, procID, blockID, rayPosIndex, d

  if ( fromRay < lbound(ed_rays,2) .or. fromRay > ubound(ed_rays,2) .or. &
       toRay < lbound(ed_rays,2) .or. toRay > ubound(ed_rays,2)) then
     print *, 'ed_meshMe:', ed_meshMe, 'fromRay:', fromRay, 'toRay:', toRay, &
          'lbound:', lbound(ed_rays,2), 'ubound:', ubound(ed_rays,2)
     call Driver_abortFlash('ed_commCheckRayLocation: out of ed_rays bounds')
  end if

  call Grid_getLocalNumBlks(numBlocks)

  do i = fromRay, toRay
     procID = int(ed_rays(RAY_PROC,i))
     if (procID /= ed_meshMe) then
        print *, 'ed_meshMe:', ed_meshMe, 'procID:', procID
        call Driver_abortFlash('ed_commCheckRayLocation: ray on wrong proc')
     end if

     blockID = int(ed_rays(RAY_BLCK,i))
     if (blockID < 1 .or. blockID > numBlocks) then
        print *, 'ed_meshMe:', ed_meshMe, 'blockID:', blockID, &
             'numBlocks:', numBlocks
        call Driver_abortFlash('ed_commCheckRayLocation: invalid block')
     end if

     call Grid_getBlkBoundBox(blockID,bndBox)

     do d = 1, NDIM
        rayPosIndex = ed_commRayPosIndex(d)
        if (ed_rays(rayPosIndex,i) < bndBox(LOW,d) .or. &
             ed_rays(rayPosIndex,i) >= bndBox(HIGH,d)) then
           print *, 'ed_meshMe:', ed_meshMe, 'ray at index:', i, &
                'with tag:', int(ed_rays(RAY_TAGS, i)), 'dimension:', d, &
                'at position:', ed_rays(rayPosIndex,i), &
                'is outside block limits:', bndBox(LOW:HIGH,d), &
                'for block:', blockID
           call Driver_abortFlash('ed_commCheckRayLocation: invalid location')
        end if
     end do
  end do
end subroutine ed_commCheckRayLocation
