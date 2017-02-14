!!****if* source/Grid/GridSolvers/MultigridMC/mg_restore_nodetypes
!!
!! NAME
!!
!!  mg_restore_nodetypes
!!
!! SYNOPSIS
!!
!!  call mg_restore_nodetypes(integer(in) :: level)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   level : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     mg_restore_nodetypes

!  Description:  restores the nodetypes and newchild arrays to pre-multigrid states    
!                calls the expensive get_tree_nodetypes
!               



subroutine mg_restore_nodetypes (level)

!===============================================================================


  use gr_mgData, ONLY: nodetype_save, newchild_save
 
  use Grid_data, ONLY : gr_meshMe,gr_meshNumProcs

  use tree, only : nodetype,newchild

  use Grid_interface,    ONLY : Grid_getLocalNumBlks
 
  implicit none

  integer, intent(in) :: level

  integer :: lnblocks, lb

  
  call Grid_getLocalNumBlks(lnblocks)

  do lb = 1, lnblocks
     nodetype(lb) = nodetype_save(lb)
     newchild(lb) = newchild_save(lb)
  enddo

  call amr_get_new_nodetypes (gr_meshNumProcs, gr_meshMe, level)

end subroutine mg_restore_nodetypes

