
! Dean Townsley 2008

#include "constants.h"

subroutine fl_fsInit

  use fl_fsData
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_interface, only : Grid_getMinCellSize, Grid_getGeometry
  implicit none
  real delta

  call Grid_getGeometry(fl_fsGeom)

  call RuntimeParameters_get("fl_fsUseConstFlameSpeed", fl_fsUseConstFlameSpeed)
  call RuntimeParameters_get("fl_fsConstFlameSpeed", fl_fsConstFlameSpeed)
  call RuntimeParameters_get("fl_fsConstFlameWidth", fl_fsConstFlameWidth)
  call RuntimeParameters_get("fl_fsUseTFI", fl_fsUseTFI)
 
  call RuntimeParameters_get("fl_fsM", fl_fsMDelta)
  call Grid_getMinCellSize(delta)
  fl_fsMDelta = fl_fsMDelta*delta
 
  call RuntimeParameters_get("fl_fsQuench", fl_fsQuench)
  call RuntimeParameters_get("fl_fsQuenchDens0", fl_fsQuenchDens0)
  call RuntimeParameters_get("fl_fsQuenchDens1", fl_fsQuenchDens1)
  fl_fsQuenchInvDDens = 1.0/(fl_fsQuenchDens1-fl_fsQuenchDens0)
 
  call RuntimeParameters_get("fl_fsGcdFlameSuppress", fl_fsGcdFlameSuppress)
  call RuntimeParameters_get("fl_fsGcdFlameSuppressTime", fl_fsGcdFlameSuppressTime)
  call RuntimeParameters_get("fl_fsGcdFlameSuppressTheta", fl_fsGcdFlameSuppressCosTheta)
  fl_fsGcdFlameSuppressCosTheta = -cos(fl_fsGcdFlameSuppressCosTheta*PI/180.0)
 
  call RuntimeParameters_get("fl_fsBuoyCompSuppress", fl_fsBuoyCompSuppress)
  call RuntimeParameters_get("fl_fsBuoyCompSuppressTime", fl_fsBuoyCompSuppressTime)
  call RuntimeParameters_get("fl_fsBuoyCompSuppressTheta", fl_fsBuoyCompSuppressCosTheta)
  fl_fsBuoyCompSuppressCosTheta = -cos(fl_fsBuoyCompSuppressCosTheta*PI/180.0)

  call fl_fsAtwoodInitTable

  call fl_fsLaminarInit

  call fl_fsTFIInit
 
end subroutine fl_fsInit
