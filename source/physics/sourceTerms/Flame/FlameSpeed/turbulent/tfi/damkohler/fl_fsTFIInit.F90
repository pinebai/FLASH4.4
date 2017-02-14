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

  call RuntimeParameters_get("fl_fsTFIBeta", fl_fsTFIBeta)
  ! go ahead and combine Beta with other coefficients for alpha
  ! Cms = 0.28 from Yeung et al. quoted by Colin et al.
  fl_fsTFIBeta = fl_fsTFIBeta * 2.0e0*log(2.0e0)/3.0e0/0.28e0

end subroutine fl_fsTFIInit
