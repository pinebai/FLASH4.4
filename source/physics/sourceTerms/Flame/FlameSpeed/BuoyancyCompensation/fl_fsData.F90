! Dean Townsley 2008

module fl_fsData

  implicit none

  integer fl_fsGeom

  real, save :: fl_fsConstFlameSpeed, fl_fsConstFlameWidth
  logical, save :: fl_fsUseConstFlameSpeed, fl_fsUseTFI

  real, save    :: fl_fsMDelta

  logical, save :: fl_fsQuench
  real, save    :: fl_fsQuenchDens0=0.0, fl_fsQuenchDens1, fl_fsQuenchInvDDens

  logical, save :: fl_fsGcdFlameSuppress
  real, save    :: fl_fsGcdFlameSuppressTime, fl_fsGcdFlameSuppressCosTheta

  logical, save :: fl_fsBuoyCompSuppress
  real, save    :: fl_fsBuoyCompSuppressTime, fl_fsBuoyCompSuppressCosTheta


end module
