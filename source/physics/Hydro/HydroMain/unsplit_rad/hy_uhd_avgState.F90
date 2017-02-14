!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_avgState
!!
!! NAME
!!
!!  hy_uhd_avgState
!!
!! SYNOPSIS
!!
!!  hy_uhd_avgState( integer(IN) :: sweepDir,
!!                   real(IN)  :: VL(HY_VARINUM3),
!!                   real(IN)  :: VR(HY_VARINUM3),
!!                   real(OUT) :: Vavg(HY_VARINUM2) )
!!
!! DESCRIPTION
!!
!!  This routine computes proper average state values at each interface
!!  using a simple arithmatic average.
!!  The calculated average state values are used in the Roe and Lax-Friedrichs
!!  Riemann solvers.
!!
!! ARGUMENTS
!!
!!  sweepDir - sweep direction
!!  VL    -  a vector array for the left state  
!!            (DENS,VELX,VELY,VELZ,MAGX,MAGY,MAGZ,PRES + GAMC,GAME,EINT)
!!  VR    -  a vector array for the right state 
!!            (DENS,VELX,VELY,VELZ,MAGX,MAGY,MAGZ,PRES + GAMC,GAME,EINT)
!!  Vavg  -  a vector array for the computed average state
!!            (DENS,VELX,VELY,VELZ,MAGX,MAGY,MAGZ,PRES + GAMC,GAME)
!!
!!*** 
#include "Flash.h"
#include "UHD.h"

Subroutine hy_uhd_avgState(sweepDir,VL,VR,Vavg)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  use Hydro_data,           ONLY : hy_forceHydroLimit
#endif
  use hy_uhd_slopeLimiters, ONLY : signum

  implicit none

#include "Flash.h"
#include "UHD.h"
  !! Arguments type declaration -----------------
  integer, intent(IN) :: sweepDir
  real, dimension(HY_VARINUM3), intent(IN)  :: VL,VR
  real, dimension(HY_VARINUM2), intent(OUT) :: Vavg
  !! --------------------------------------------
  real :: sig

  !! Arithmetic averages
  Vavg(HY_DENS:HY_GAME) = .5*(VL(HY_DENS:HY_GAME)+VR(HY_DENS:HY_GAME))

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  if (hy_forceHydroLimit) Vavg(HY_MAGX:HY_MAGZ) = 0.
#endif

  !! Use upwinding for game and gamc that are averaged along the
  !! entropy wave.
  sig = signum(Vavg(sweepDir+1))
  Vavg(HY_GAMC:HY_GAME) = 0.5*( (1.+sig)*VL(HY_GAMC:HY_GAME) &
                               +(1.-sig)*VR(HY_GAMC:HY_GAME) )

End Subroutine hy_uhd_avgState
