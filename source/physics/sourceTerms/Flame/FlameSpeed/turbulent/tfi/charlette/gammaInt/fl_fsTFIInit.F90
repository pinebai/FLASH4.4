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

  call RuntimeParameters_get("fl_fsTFIPrandtl", fl_fsTFIPrandtl)
  call RuntimeParameters_get("fl_fsTFIetol", fl_fsTFIetol)
  call RuntimeParameters_get("fl_fsTFIBeta", fl_fsTFIBeta)
  if ((fl_fsTFIBeta .lt. 0.0) .or. (fl_fsTFIBeta .gt. 1.0)) &
     call Driver_abortFlash("fl_fsTFIBeta must be between 0 and 1 for Charlette model")

end subroutine fl_fsTFIInit
