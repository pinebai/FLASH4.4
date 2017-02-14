!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonUtilities/pi_updateProtons
!!
!! NAME
!!
!!  pi_updateProtons
!!
!! SYNOPSIS
!!
!!  call pi_updateProtons (logical (in) :: doMove)
!!
!! DESCRIPTION
!!
!!  This routine updates the block ID info of the current set of protons.
!!  It consists of four stages:
!!
!!                   i) sort proton array on old block ID
!!                  ii) move proton array to new block ID
!!                 iii) sort proton array on new block ID
!!                  iv) extract new block ID info
!!
!!  When it is certain, that the old blockID's are still valid, then
!!  there is the option of calling this routine and having it execute
!!  only steps i) and iv). All nonexistent protons will be removed (i.e.
!!  placed at the very end of the protons array).
!!
!! ARGUMENTS
!!
!!  doMove : logical keyword to invoke movement of the protons (if false -> no move)
!!
!! NOTES
!!
!!***

subroutine pi_updateProtons (doMove)

  use Driver_interface,   ONLY : Driver_abortFlash

  use Grid_interface,     ONLY : Grid_moveParticles, &
                                 Grid_sortParticles

  use pi_interface,       ONLY : pi_printProtonsData,  &
                                 pi_protonsBlockIDInfo

  use ProtonImaging_data, ONLY : pi_globalComm,         &
                                 pi_globalMe,           &
                                 pi_particleIndexCount, &
                                 pi_particleIndexList,  &
                                 pi_maxProtonCount,     &
                                 pi_protonCount,        &
                                 pi_protons,            &
                                 pi_protonDeterminism

  implicit none

#include "constants.h"
#include "Flash.h"
#include "ProtonImaging.h"

  include "Flash_mpi.h"

  logical, intent (in) :: doMove

  integer :: protonsPerBlk  (1:MAXBLOCKS,1:1)
!
!
!     ...1st step (easy): sort the protons array. All nonexistent protons get moved to
!                         the end of the array. The number of active protons (pi_protonCount)
!                         gets updated. The block ID of the protons is meaningless at
!                         this stage, since the protons have been moved through that
!                         block and are now in a different block.
!
!
  call Grid_sortParticles (pi_protons,        &
                           PROTON_ATTRCOUNT,  &
                           pi_protonCount,    &
                           1,                 &
                           pi_maxProtonCount, &
                           protonsPerBlk,     &
                           PROTON_BLCK        )
!
!
!     ...2nd step (hard): determines the new block ID of the protons using their
!                         current coordinates. Since each processor is associated with
!                         a fixed set of blocks (and thus their ID's), the number of
!                         active protons (pi_protonCount) gets updated again.
!
!
  if (doMove) then

!      call pi_printProtonsData ('beforeGridMove',pi_globalMe)

      call Grid_moveParticles (pi_protons,            &
                               PROTON_ATTRCOUNT,      &
                               pi_maxProtonCount,     &
                               pi_protonCount,        &
                               pi_particleIndexList,  &
                               pi_particleIndexCount, &
                               pi_protonDeterminism   )

!      call pi_printProtonsData ('afterGridMove',pi_globalMe)

!     call mpi_barrier (pi_globalComm, error)
!     call Driver_abortFlash ("Stopped FLASH smoothly!")
!
!
!     ...3rd step (easy): sort the protons array again (just in case, maybe this is not
!                         needed and the 2nd step takes care of that automatically.
!
!

      if (pi_protonCount > 0) then

          call Grid_sortParticles (pi_protons,        &
                                   PROTON_ATTRCOUNT,  &
                                   pi_protonCount,    &
                                   1,                 &
                                   pi_maxProtonCount, &
                                   protonsPerBlk,     &
                                   PROTON_BLCK        )
      end if

  end if
!
!
!     ...4th step (easy): extract the block info of the protons array to get ready
!                         for the block proton tracing step.
!
!
  if (pi_protonCount > 0) then
      call pi_protonsBlockIDInfo ()
  end if
!
!
!     ...Ready!
!
!
  return
end subroutine pi_updateProtons
