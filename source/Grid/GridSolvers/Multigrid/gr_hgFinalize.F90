!!****if* source/Grid/GridSolvers/Multigrid/gr_hgFinalize
!!
!! NAME
!!  gr_hgFinalize
!!
!! SYNOPSIS
!!
!!  gr_hgFinalize()
!!
!! DESCRIPTION
!!
!!  Deallocate any memory that might have been allocated in Multigrid unit
!!  and other housekeeping to prepare for the end of the unit.
!!
!!***

subroutine gr_hgFinalize

  !==================================================================
  use gr_hgData, ONLY: gr_hgSaveNodeType, gr_hgSaveNewChild, Px, Py, Pz, &
       send_prolong_data, send_prolong_req, send_restrict_data, &
        send_restrict_req, recv_prolong_data, recv_restrict_data
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer :: ierr1, ierr2, ierr3

  deallocate(gr_hgSaveNodetype, stat=ierr1)
  deallocate(gr_hgSaveNewChild, stat=ierr2)
  
  if ((ierr1 /= 0) .or. (ierr2 /= 0)) &
       call Driver_abortFlash("Deallocation error for gr_hgSaveNodetype or gr_hgSaveNewChild")
  
  
  deallocate(Px, stat=ierr1)
  deallocate(Py, stat=ierr2)
  deallocate(Pz, stat=ierr3)
  if ((ierr1 /= 0) .or. (ierr2 /= 0) .or. (ierr3 /= 0)) &
       call Driver_abortFlash("Deallocation error for the P<dim> arrays")


  deallocate(send_prolong_data, stat=ierr1)
  deallocate(recv_prolong_data, stat=ierr2)
  deallocate(send_prolong_req, stat=ierr3)

  if ((ierr1 /= 0) .or. (ierr2 /= 0) .or. (ierr3 /= 0)) &
       call Driver_abortFlash("Allocation error for prolongation")
  

  deallocate(send_restrict_data, stat=ierr1)
  deallocate(recv_restrict_data, stat=ierr2)
  deallocate(send_restrict_req, stat=ierr3)

  if ((ierr1 /= 0) .or. (ierr2 /= 0) .or. (ierr3 /= 0)) &
       call Driver_abortFlash("Deallocation error for restriction")
  

  return
end subroutine gr_hgFinalize
