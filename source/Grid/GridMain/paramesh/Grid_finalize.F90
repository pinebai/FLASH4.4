!!****if* source/Grid/GridMain/paramesh/Grid_finalize
!!
!! NAME
!!
!!  Grid_finalize
!!
!!
!! SYNOPSIS
!!
!!  Grid_finalize()
!!
!!
!! DESCRIPTION
!!
!!  Deallocates memory allocated in the Grid Unit to prepare for shutdowns
!!
!!***

subroutine Grid_finalize()

  use gr_bcInterface, ONLY : gr_bcFinalize
  use gr_ptInterface, ONLY : gr_ptFinalize
  use Grid_data, ONLY : gr_nToLeft, gr_gid
#ifdef FLASH_GRID_PARAMESH3OR4
  use Grid_data, ONLY : gr_gsurr_blks
#endif
  use gr_sbInterface, ONLY: gr_sbFinalize

  implicit none

  if(allocated(gr_gid))deallocate(gr_gid)
  if(allocated(gr_nToLeft))deallocate(gr_nToLeft)
#ifdef FLASH_GRID_PARAMESH3OR4
  if(allocated(gr_gsurr_blks))deallocate(gr_gsurr_blks)
#endif

  call gr_ptFinalize()
  call gr_solversFinalize()
  call gr_bcFinalize()
!  call gr_sbFinalize()
  call Paramesh_finalize()

end subroutine Grid_finalize
