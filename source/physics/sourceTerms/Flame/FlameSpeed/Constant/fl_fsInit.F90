
! Dean Townsley 2008

subroutine fl_fsInit

  use fl_fsData
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

  call RuntimeParameters_get("fl_fsConstFlameSpeed", fl_fsConstFlameSpeed)
  call RuntimeParameters_get("fl_fsConstFlameWidth", fl_fsConstFlameWidth)

end subroutine fl_fsInit
