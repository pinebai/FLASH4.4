!!****if* source/Grid/GridMain/paramesh/gr_xyzToBlock
!!
!! NAME
!!  gr_xyzToBlock
!!
!! SYNOPSIS
!!
!!  call gr_xyzToBlock( real (IN)  :: xyz(MDIM),
!!                    integer(OUT) :: procID,
!!                    integer(OUT) :: blkID)
!!  
!! DESCRIPTION 
!!  
!!  This routine returns the identity of the block, on
!!  which the specified physical location falls. It also
!!  returns the ID of the processor on which the block is
!!  residing.
!!
!!  If the physical location lies outside the domain
!!  boundaries, then procID=NONEXISTENT will be returned;
!!  boundary periodicity is not considered.
!!
!! ARGUMENTS 
!!
!! NOTES
!!
!!  This implementation requires BITTREE.
!!***

#include "constants.h"
#include "Flash.h"

Subroutine gr_xyzToBlock(xyz, procID, blkID)
  use gr_interface, ONLY : gr_xyzToBlockLevel
  use Grid_data, ONLY : gr_meshNumProcs
  use Driver_interface, ONLY : Driver_abortFlash
  use tree, ONLY : lrefine_max
#ifdef BITTREE
  use bittree, only : amr_identify_block
#endif

  implicit none

  real, dimension(MDIM),intent(IN) :: xyz
  integer, intent(OUT) :: procID
  integer, intent(OUT) :: blkID

#ifdef BITTREE
  integer, dimension(MDIM) :: ijk
#endif
  integer :: lev, proc, blk

  lev = lrefine_max

#ifdef BITTREE  
  call gr_xyzToBlockLevel(lev, xyz, ijk)
  call amr_identify_block(gr_meshNumProcs, lev, ijk, proc, blk)
  if (lev < 1) then       !if the point was outside of the domain...
     proc = NONEXISTENT         !make sure we communicate nonexistence to the caller.
  end if
#else
  proc = -size(xyz); blk= -1  !avoid compiler warnings about uninitialized variables below
  call Driver_abortFlash("gr_xyzToBlock works only when bittree is enabled")
#endif

  procID=proc
  blkID=blk
End Subroutine gr_xyzToBlock
