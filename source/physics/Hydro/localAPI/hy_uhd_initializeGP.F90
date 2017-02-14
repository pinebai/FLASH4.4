Subroutine hy_uhd_counterGP_init(radiusGP, counterGP, blkLimitsGC)

  implicit none

#include "constants.h"

  !!-----Arguments------------------------------------------------
  real,    intent(IN)  :: radiusGP
  integer, intent(OUT) :: counterGP
  integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimitsGC
  !!--------------------------------------------------------------
  counterGP = 0
end subroutine hy_uhd_counterGP_init



! *****************************************************************
Subroutine hy_uhd_initGP(RinvGP, WpGP, WmGP, blkLimitsGC)

! * Main init of GP interpolation. Here we set up one-time 
! * initialized arrays and matrices. This is the only place
! * where we invert R_mn(counter,counter). 
! *****************************************************************

  implicit none


  !!-----Arguments------------------------------------------------
  real, intent(INOUT), dimension(:,:) :: RinvGP
  real, intent(INOUT), dimension(:,:) :: WpGP
  real, intent(INOUT), dimension(:,:) :: WmGP
  integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimitsGC
  !!--------------------------------------------------------------

end subroutine hy_uhd_initGP
