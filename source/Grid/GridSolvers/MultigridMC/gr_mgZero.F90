!!****if* source/Grid/GridSolvers/MultigridMC/gr_mgZero
!!
!! NAME
!!
!!  gr_mgZero
!!
!! SYNOPSIS
!!
!!  call gr_mgZero(integer(in) :: level,
!!                 integer(in) :: ivar,
!!                 integer(in) :: leaf_only)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   level : 
!!
!!   ivar : 
!!
!!   leaf_only : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     mg_zero()

!  Description: Initializes (zeroes) a multigrid variable on a particular
!               mesh level.

!  Parameters:  level       Level to zero on.
!               ivar        Variable to zero.
!               leaf_only   If nonzero, only zero leaf-node blocks on this
!                           level.


subroutine gr_mgZero (level, ivar, leaf_only)

!===============================================================================

use gr_mgData, ONLY: nodetype_save, ili, jli, kli, iui, jui, kui,&
     ile, jle, kle, iue, jue, kue

use tree, only : lrefine
use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                              Grid_getBlkPtr,       &
                              Grid_releaseBlkPtr

implicit none
#include "constants.h"

integer, intent(in) :: level, ivar, leaf_only
integer :: lnblocks
integer :: lb, i, j, k
real, pointer, dimension(:,:,:,:) :: unk

!===============================================================================

call Grid_getLocalNumBlks(lnblocks)

if (leaf_only == 0) then

  do lb = 1, lnblocks
    if (lrefine(lb) == level) then

      ! Point to blocks center vars:
      call Grid_getBlkPtr(lb,unk,CENTER)

      do k = kli, kui
        do j = jli, jui
          do i = ili, iui
            unk(ivar,i,j,k) = 0.
          enddo
        enddo
      enddo

      ! Release pointers:
      call Grid_releaseBlkPtr(lb,unk,CENTER)

    endif
  enddo

elseif (leaf_only == 1) then

  do lb = 1, lnblocks
    if ((lrefine(lb) == level) .and. (nodetype_save(lb) == 1)) then

      ! Point to blocks center vars:
      call Grid_getBlkPtr(lb,unk,CENTER)
      do k = kli, kui
        do j = jli, jui
          do i = ili, iui
            unk(ivar,i,j,k) = 0.
          enddo
        enddo
      enddo
      ! Release pointers:
      call Grid_releaseBlkPtr(lb,unk,CENTER)


    endif
  enddo

else

  do lb = 1, lnblocks
    if ((lrefine(lb) == level) .and. (nodetype_save(lb) /= 1)) then 

      ! Point to blocks center vars:
      call Grid_getBlkPtr(lb,unk,CENTER)

      unk(ivar,ili:iui,jli:jui,kli:kui) = 0.

      ! Release pointers:
      call Grid_releaseBlkPtr(lb,unk,CENTER)
 
    endif

  enddo

endif

!===============================================================================

return
end
