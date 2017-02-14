!
! Aaron Jackson 2008
!
! deallocates the arrays containing the flame speed table
!
subroutine fl_fsLaminarFinalize

  use fl_fsConeData, ONLY : ldenb, ne22ib, c12ib, stab, dstab
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer :: istat

  deallocate(ldenb,stat=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot deallocate ldenb in fl_fsConeEndTable")
  deallocate(ne22ib,stat=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot deallocate ne22ib in fl_fsConeEndTable")
  deallocate(c12ib,stat=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot deallocate c12ib in fl_fsConeEndTable")

  deallocate(stab,stat=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot deallocate stab in fl_fsConeEndTable")
  deallocate(dstab,stat=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot deallocate dstab in fl_fsConeEndTable")

  return
end subroutine
