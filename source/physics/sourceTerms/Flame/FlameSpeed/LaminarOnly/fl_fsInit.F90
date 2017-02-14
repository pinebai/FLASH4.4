
! Dean Townsley, Aaron Jackson, Alan Calder 2008

subroutine fl_fsInit

  use fl_fsData
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

  call RuntimeParameters_get("fl_fsUseConstFlameSpeed", fl_fsUseConstFlameSpeed)
  call RuntimeParameters_get("fl_fsConstFlameSpeed", fl_fsConstFlameSpeed)
  call RuntimeParameters_get("fl_fsConstFlameWidth", fl_fsConstFlameWidth)
  call RuntimeParameters_get("fl_fsUseTFI", fl_fsUseTFI)

  call fl_fsLaminarInit

  call fl_fsTFIInit

end subroutine fl_fsInit

