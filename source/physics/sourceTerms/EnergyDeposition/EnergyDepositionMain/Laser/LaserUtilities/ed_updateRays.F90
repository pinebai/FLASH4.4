!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_updateRays
!!
!! NAME
!!
!!  ed_updateRays
!!
!! SYNOPSIS
!!
!!  call ed_updateRays (logical (in) :: doMove)
!!
!! DESCRIPTION
!!
!!  This routine updates the block ID info of the current set of rays.
!!  It consists of four stages:
!!
!!                   i) sort ray array on old block ID
!!                  ii) move ray array to new block ID
!!                 iii) sort ray array on new block ID
!!                  iv) extract new block ID info
!!
!!  When it is certain, that the old blockID's are still valid, then
!!  there is the option of calling this routine and having it execute
!!  only steps i) and iv). All nonexistent rays will be removed (i.e.
!!  placed at the very end of the rays array).
!!
!! ARGUMENTS
!!
!!  doMove : logical keyword to invoke movement of the rays (if false -> no move)
!!
!! NOTES
!!
!!***

subroutine ed_updateRays (doMove)

  use ed_commInterface,      ONLY : ed_commGetNewRays,        &
                                    ed_commProgressTransport, &
                                    ed_commSortParticles

  use ed_interface,          ONLY : ed_raysBlockIDInfo

  use EnergyDeposition_data, ONLY : ed_maxRayCount,        &
                                    ed_rayCount,           &
                                    ed_rays

  implicit none

#include "constants.h"
#include "Flash.h"
#include "EnergyDeposition.h"

  logical, intent (in) :: doMove

  logical, parameter :: forceProgress = .true.
!
!
!     ...1st step (easy): sort the rays array. All nonexistent rays get moved to
!                         the end of the array. The number of active rays (ed_rayCount)
!                         gets updated. The block ID of the rays is meaningless at
!                         this stage, since the rays have been moved through that
!                         block and are now in a different block.
!
!
  call ed_commSortParticles ()
!
!
!     ...2nd step (hard): determines the new block ID of the rays using their
!                         current coordinates. Since each processor is associated with
!                         a fixed set of blocks (and thus their ID's), the number of
!                         active rays (ed_rayCount) gets updated again.
!
!
  if (doMove) then

      call ed_commProgressTransport (optionalForceProgress = forceProgress)

      call ed_commGetNewRays ()

!
!
!     ...3rd step (easy): sort the rays array again (just in case, maybe this is not
!                         needed and the 2nd step takes care of that automatically.
!
!
      if (ed_rayCount > 0) then
         call ed_commSortParticles ()
      end if

  end if
!
!
!     ...4th step (easy): extract the block info of the rays array to get ready
!                         for the block ray tracing step.
!
!
  if (ed_rayCount > 0) then
      call ed_raysBlockIDInfo ()
  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ed_updateRays
