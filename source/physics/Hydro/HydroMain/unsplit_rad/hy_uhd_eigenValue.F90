!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_eigenValue
!!
!! NAME
!!
!!  hy_uhd_eigenValue
!!
!! SYNOPSIS
!!
!!  hy_uhd_eigenValue( real (OUT)          :: EigValue(HY_WAVENUM),
!!                     real (IN)           :: U_normal,
!!                     real (IN)           :: C_fast,
!!                     real (IN), optional :: C_alfn,
!!                     real (IN), optional :: C_slow,
!!                     real (OUT),optional :: C_hyp)
!!                     
!!
!! DESCRIPTION
!!
!!  This routine calculates MHD/Hydro eigenvalues.
!!
!! ARGUMENTS
!!
!!  EigValue   - Eigenvalues (wave speeds)
!!  U_normal   - Fluid velocity in normal direction
!!  C_fast     - Fast magnetoacoustic speed for MHD/Sound speed for Hydro
!!  C_alfn     - Alfven speed (needed for MHD only)
!!  C_slow     - Slow magnetoacoustic speed (needed for MHD only)
!!  C_hyp    - advection wave speed for GLM-MHD
!!
!!***


Subroutine hy_uhd_eigenValue(EigValue,U_normal,C_fast,C_alfn,C_slow,C_hyp)
    
  implicit none

#include "Flash.h"
#include "UHD.h"

  !! Arguments type declaration --------------------------
  real,dimension(HY_WAVENUM), intent(OUT) :: EigValue
  real,intent(IN) :: U_normal,C_fast
  real,intent(IN), optional :: C_alfn,C_slow,C_hyp
  !! -----------------------------------------------------

  EigValue(HY_FASTLEFT) = U_normal-C_fast
  EigValue(HY_SLOWLEFT) = U_normal
  EigValue(HY_ENTROPY)  = U_normal
  EigValue(HY_SLOWRGHT) = U_normal
  EigValue(HY_FASTRGHT) = U_normal+C_fast

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  EigValue(HY_ALFNLEFT) = U_normal-C_alfn
  EigValue(HY_ALFNRGHT) = U_normal+C_alfn
  EigValue(HY_SLOWLEFT) = EigValue(HY_SLOWLEFT)-C_slow
  EigValue(HY_SLOWRGHT) = EigValue(HY_SLOWRGHT)+C_slow
#ifdef FLASH_UGLM_MHD
  EigValue(HY_GLMPLEFT) =-C_hyp
  EigValue(HY_GLMPRGHT) = C_hyp
#endif
#endif


End Subroutine hy_uhd_eigenValue
