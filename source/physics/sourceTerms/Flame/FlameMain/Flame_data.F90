!  Dean Townsley 2008

#include "Flash.h"
Module Flame_data

  implicit none

  ! runtime parameters
  logical, save :: fl_useFlame
  real, save :: fl_epsilon_0, fl_epsilon_1, fl_kpp_fact, fl_b
  real, save :: fl_initProfileAdjustWidth

  ! some constants for diffusion-reaction equation
  real, save :: fl_R_over_s, fl_kappa_over_s, fl_width

  integer, save :: fl_gcMaskSize = NUNK_VARS
  logical, save,dimension(NUNK_VARS) :: fl_gcMask
  logical, save :: fl_gcDoEos, fl_gcDoLogMask

end module
