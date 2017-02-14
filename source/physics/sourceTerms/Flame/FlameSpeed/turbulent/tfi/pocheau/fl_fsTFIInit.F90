! Aaron Jackson 2010

#include "constants.h"

subroutine fl_fsTFIInit

  use fl_fsTFIData
  use fl_fsData, ONLY : fl_fsUseTFI
  use Turb_data, ONLY : turb_useTurb
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

  if (fl_fsUseTFI .and. (.not. turb_useTurb) ) &
     call Driver_abortFlash("fl_fsUseTFI requires turb_useTurb")

  call RuntimeParameters_get("fl_fsTFICt", fl_fsTFICt)

end subroutine fl_fsTFIInit
