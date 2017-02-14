!!****if* source/Grid/GridSolvers/Multigrid/gr_hgRestoreNodeTypes
!!
!! NAME
!!  gr_hgRestoreNodeTypes
!!
!! SYNOPSIS
!!  
!!  call gr_hgRestoreNodeTypes()
!!
!! DESCRIPTION
!!
!!  Restore the original nodetypes using the gr_hgSaveNodetype and
!!  gr_hgSaveNewchild arrays
!!
!! NOTE
!!
!!  We can define CAUTIOUS_PARAMESH_CALL to use the more expensive
!!  amr_morton_process().  We have run more tests using this call with 
!!  Multigrid.
!!
!!***

#include "Flash.h"

subroutine gr_hgRestoreNodeTypes ()

  use Grid_data, ONLY : gr_meshNumProcs, gr_meshMe
  use gr_hgdata, ONLY: gr_hgSaveNodetype, gr_hgSaveNewchild, gr_hgMeshRefineMax
  use tree, ONLY : lnblocks, newchild,nodetype, lrefine_max

  use Timers_interface, ONLY: Timers_start, Timers_stop

  implicit none

  integer :: level

  nodetype(1:lnblocks) = gr_hgSaveNodetype(1:lnblocks)
  newchild(1:lnblocks) = gr_hgSaveNewchild(1:lnblocks)

  call Timers_start("amr_morton_process")

#ifdef FLASH_GRID_PARAMESH2
  call get_tree_nodetypes (gr_meshNumProcs, gr_meshMe)
#else

#ifdef CAUTIOUS_PARAMESH_CALL
  call amr_morton_process()
#else
  call amr_get_new_nodetypes (gr_meshNumProcs, gr_meshMe, gr_hgMeshRefineMax)
#endif

#endif

  call Timers_stop("amr_morton_process")

end subroutine gr_hgRestoreNodeTypes

