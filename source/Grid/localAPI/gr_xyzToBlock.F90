!!****if* source/Grid/localAPI/gr_xyzToBlock
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
!!  returns the processor ID on which the block is residing.
!!  
!!
!!
!! ARGUMENTS 
!!
!! NOTES
!!
!!  This is a stub implementation that does not do anything useful.
!!
!!  A more useful implementation is provided by PARAMESH Grid
!!  implementations. It requires PARAMESH4DEV with BITTREE.
!!***

#include "constants.h"

Subroutine gr_xyzToBlock(xyz, procID, blkID)

  implicit none

  real, dimension(MDIM),intent(IN) :: xyz
  integer, intent(OUT) :: procID
  integer, intent(OUT) :: blkID

  procID=0
  blkID=0  
End Subroutine gr_xyzToBlock
