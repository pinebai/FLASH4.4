!!****if* source/Grid/GridSolvers/MultigridMC/gr_mgRestrict
!!
!! NAME
!!
!!  gr_mgRestrict
!!
!! SYNOPSIS
!!
!!  call gr_mgRestrict(integer(in) :: level,
!!                     integer(in) :: ifrom,
!!                     integer(in) :: ito)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   level : 
!!
!!   ifrom : 
!!
!!   ito : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     mg_restrict()

!  Description: Restrict a multigrid variable from one mesh level to the
!               next coarser level.  Assume that child blocks on the "from"
!               level and leaf blocks on the "to" level contain valid data.
!               On output, all blocks on the "to" level should then contain
!               valid data.

!  Parameters:  level       Level to restrict from.
!               ifrom       Index of the variable to restrict from.
!               ito         Index of the variable to restrict into.


subroutine gr_mgRestrict (level, ifrom, ito)

!===============================================================================

  use gr_mgData, ONLY: nodetype_save, ili, jli, kli, iui, jui, kui

  use Grid_data, ONLY : gr_meshMe,gr_meshNumProcs

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr

  use tree, only : nodetype,lrefine

  use workspace , ONLY: work

  implicit none
#include "constants.h"
  integer, intent(in) :: level, ifrom, ito

  integer :: indexPlotNumber, PTNumber, ierr
  integer :: lb, lnblocks 
  real    :: normf, normt
  real, pointer, dimension(:,:,:,:),save :: unk

!===============================================================================

!               The PARAMESH library's restriction and boundary-update routines
!               only operate on all arrays at once or on the "work" array.
!               So we copy the input variable to work(), restrict work(),
!               then copy to the output variable.  However, we have to trick
!               the mesh package into thinking all of the "from"-level blocks
!               are leaf nodes, so that the supplied restriction routine can
!               do its work.

!               Temporarily alter node type arrays.

  call Grid_getLocalNumBlks(lnblocks)

  do lb = 1, lnblocks
     if ( ((lrefine(lb) == level-1)  .and. &
           (nodetype_save(lb) == 1) ) .or. &
          (lrefine(lb) == level) ) then

          ! Point to blocks center vars:
          call Grid_getBlkPtr(lb,unk,CENTER)

          work(ili:iui,jli:jui,kli:kui,lb,1) = &
                   unk(ifrom,ili:iui,jli:jui,kli:kui)

          ! Point to blocks center vars:
          call Grid_releaseBlkPtr(lb,unk,CENTER)

     endif
  end do

! alter nodetypes

  call amr_mg_restrict(gr_meshNumProcs, gr_meshMe, level)

!               Copy the result into the requested output level and variable.
  do lb = 1, lnblocks
     if (lrefine(lb) == level-1) then

          ! Point to blocks center vars:
          call Grid_getBlkPtr(lb,unk,CENTER)
 
          unk(ito,ili:iui,jli:jui,kli:kui) = & 
          &      work(ili:iui,jli:jui,kli:kui,lb,1)

          ! Point to blocks center vars:
          call Grid_releaseBlkPtr(lb,unk,CENTER)

     endif
  enddo

!===============================================================================

  return
end subroutine gr_mgRestrict
