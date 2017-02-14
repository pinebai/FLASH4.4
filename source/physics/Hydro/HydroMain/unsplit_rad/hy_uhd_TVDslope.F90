!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_TVDslope
!!
!! NAME
!!
!!  hy_uhd_TVDslope
!!
!! SYNOPSIS
!!
!!  hy_uhd_TVDslope(integer(IN) :: dir,
!!                  real(IN)    :: VL (HY_VARINUMMAX),
!!                  real(IN)    :: V0 (HY_VARINUMMAX),
!!                  real(IN)    :: VR (HY_VARINUMMAX),
!!                  real(IN)    :: lambda (HY_WAVENUM),
!!                  real(IN)    :: leig(HY_VARINUM,HY_WAVENUM),
!!                  real(OUT)   :: delbar(HY_VARINUMMAX))
!!
!! DESCRIPTION
!!
!!   This routine calculates limited slopes depending on the user's choice of the slope limiter. 
!!   Possible choices are 
!!
!!      "minmod", 
!!      "vanLeer", 
!!      "mc", 
!!      "hybrid" (see Balsara, ApJ. Suppl, 2004), 
!!      "limited" (see Toro, pp 481, 2nd Ed, for detail).
!!
!!   The hybrid limiter is a combination of any limiters, preferably 
!!   minmod and vanLeer (or mc) limiters applied to different
!!   characteristic/primitive variables. Users can choose to limit
!!   either primitive or characteristic variables. In general,
!!   the characteristic limiting provides more accurate solutions
!!   than the primitive limiting near discontinuities.
!!   
!!
!! ARGUMENTS
!!
!!   dir         - sweep direction
!!   VL          - primitive variables at the first left  neighboring cell
!!   V0          - cell centered variables in primitive form
!!   VR          - primitive variables at the first right neighboring cell
!!   lambda      - left eigenvalue evaluated at the current cell
!!   leig        - left eigenvector evaluated at the current cell
!!   delbar      - limited slopes for primitive variables + gamc,game,gravity
!!
!!
!! REFERENCES
!!
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!  * Balsara, ApJ. Suppl, 151:149--184, 2004
!!
!!***

Subroutine hy_uhd_TVDslope(dir,VL,V0,VR,lambda,leig,delbar)

  use hy_uhd_slopeLimiters, ONLY : minmod,mc,vanLeer
  use Hydro_data,           ONLY : hy_limiter,hy_charLimiting,hy_LimitedSlopeBeta

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Arguments type declaration ------------------------------------------
  integer, intent(IN) :: dir
  real, dimension(HY_VARINUMMAX),intent(IN)  :: VL,V0,VR
  real, dimension(HY_WAVENUM),   intent(IN)  :: lambda
  real, dimension(HY_VARINUM,HY_WAVENUM), intent(IN) :: leig
  real, dimension(HY_VARINUMMAX),intent(OUT) :: delbar
  !! ---------------------------------------------------------------------
  integer :: k,extraVarinum
  real    :: del1,del2
  real, dimension(HY_VARINUMMAX) :: delV_L, delV_R

  delbar(:) = 0.0
  extraVarinum = HY_END_VARS-HY_GAMC+1

  ! Define undivided deltas
  delV_L(HY_DENS:HY_VARINUMMAX) = V0(HY_DENS:HY_VARINUMMAX)-VL(HY_DENS:HY_VARINUMMAX)
  delV_R(HY_DENS:HY_VARINUMMAX) = VR(HY_DENS:HY_VARINUMMAX)-V0(HY_DENS:HY_VARINUMMAX)



  ! Get slope limiter
  if (hy_limiter==MINMOD) then
     if (hy_charLimiting) then

        do k=1,HY_WAVENUM
           delbar(k)=minmod(dot_product(leig(1:HY_VARINUM,k),delV_L(HY_DENS:HY_VARINUM)),&
                            dot_product(leig(1:HY_VARINUM,k),delV_R(HY_DENS:HY_VARINUM)))
        enddo

!#ifndef FLASH_EOS_GAMMA
        ! for gamc, game, eint, and gravity
        do k=HY_GAMC,HY_END_VARS
           !*** MC (monotonized  central) limiter ***
           delbar(k)=minmod(delV_L(k),delV_R(k))
        enddo
!#endif
     else
        do k=1,HY_VARINUM+extraVarinum
           !*** MinMod limiter ***
           delbar(k)=minmod(delV_L(k),delV_R(k))

        enddo
     endif

  elseif (hy_limiter==MC) then
     if (hy_charLimiting) then
        do k=1,HY_WAVENUM
           delbar(k)=mc(dot_product(leig(1:HY_VARINUM,k),delV_L(HY_DENS:HY_VARINUM)),&
                        dot_product(leig(1:HY_VARINUM,k),delV_R(HY_DENS:HY_VARINUM)))
        enddo
!#ifndef FLASH_EOS_GAMMA
        ! for gamc, game, eint, and gravity
        do k=HY_GAMC,HY_END_VARS
           !*** MC (monotonized  central) limiter ***
           delbar(k)=mc(delV_L(k),delV_R(k))          
        enddo
!#endif
     else   
        do k=1,HY_VARINUM+extraVarinum
           !*** MC (monotonized  central) limiter ***
           delbar(k)=mc(delV_L(k),delV_R(k))
        enddo
     endif

  elseif (hy_limiter==VANLEER) then
     if (hy_charLimiting) then
        do k=1,HY_WAVENUM
           delbar(k)=vanLeer(dot_product(leig(1:HY_VARINUM,k),delV_L(HY_DENS:HY_VARINUM)),&
                             dot_product(leig(1:HY_VARINUM,k),delV_R(HY_DENS:HY_VARINUM)))
        enddo
!#ifndef FLASH_EOS_GAMMA
        ! for gamc, game, eint, and gravity
        do k=HY_GAMC,HY_END_VARS
           !*** MC (monotonized  central) limiter ***
           delbar(k)=vanLeer(delV_L(k),delV_R(k))
        enddo
!#endif
     else
        do k=1,HY_VARINUM+extraVarinum
           !*** vanLeer limiter ***
           delbar(k)=vanLeer(delV_L(k),delV_R(k))
        enddo
     endif

  elseif (hy_limiter==HYBRID) then
     if (hy_charLimiting) then
        !minmod for genuinely nonlinear waves and vanLeer for lineary degenerate waves
        delbar(HY_FASTLEFT)=&
             minmod(dot_product(leig(1:HY_VARINUM,HY_FASTLEFT),delV_L(HY_DENS:HY_VARINUM)),&
                    dot_product(leig(1:HY_VARINUM,HY_FASTLEFT),delV_R(HY_DENS:HY_VARINUM)))
#ifdef FLASH_USM_MHD
        delbar(HY_ALFNLEFT)=&
                 mc(dot_product(leig(1:HY_VARINUM,HY_ALFNLEFT),delV_L(HY_DENS:HY_VARINUM)),&
                    dot_product(leig(1:HY_VARINUM,HY_ALFNLEFT),delV_R(HY_DENS:HY_VARINUM)))
#endif
        delbar(HY_SLOWLEFT)=&
             minmod(dot_product(leig(1:HY_VARINUM,HY_SLOWLEFT),delV_L(HY_DENS:HY_VARINUM)),&
                    dot_product(leig(1:HY_VARINUM,HY_SLOWLEFT),delV_R(HY_DENS:HY_VARINUM)))

        delbar(HY_ENTROPY)=&
                 mc(dot_product(leig(1:HY_VARINUM,HY_ENTROPY), delV_L(HY_DENS:HY_VARINUM)),&
                    dot_product(leig(1:HY_VARINUM,HY_ENTROPY), delV_R(HY_DENS:HY_VARINUM)))

        delbar(HY_SLOWRGHT)=&
             minmod(dot_product(leig(1:HY_VARINUM,HY_SLOWRGHT),delV_L(HY_DENS:HY_VARINUM)),&
                    dot_product(leig(1:HY_VARINUM,HY_SLOWRGHT),delV_R(HY_DENS:HY_VARINUM)))
#ifdef FLASH_USM_MHD
        delbar(HY_ALFNRGHT)=&
                 mc(dot_product(leig(1:HY_VARINUM,HY_ALFNRGHT),delV_L(HY_DENS:HY_VARINUM)),&
                    dot_product(leig(1:HY_VARINUM,HY_ALFNRGHT),delV_R(HY_DENS:HY_VARINUM)))
#endif
        delbar(HY_FASTRGHT)=&
             minmod(dot_product(leig(1:HY_VARINUM,HY_FASTRGHT),delV_L(HY_DENS:HY_VARINUM)),&
                    dot_product(leig(1:HY_VARINUM,HY_FASTRGHT),delV_R(HY_DENS:HY_VARINUM)))
!#ifndef FLASH_EOS_GAMMA
        ! for gamc, game, eint, and gravity
        do k=HY_GAMC,HY_END_VARS
           !*** MC (monotonized  central) limiter ***
           delbar(k)=mc(delV_L(k),delV_R(k))
        enddo
!#endif
     else   
        !minmod for velocity fields and pressure
        !(i.e., four magnetosonic genuinely nonlinear waves)
        !compressible limiter for density and magnetic fields 
        !(i.e., three linearly degenerate waves: two Alfven waves, one entropy wave)

        do k=1,HY_VARINUM+extraVarinum
           select case (k)
#ifdef FLASH_USM_MHD
           case (HY_DENS,HY_MAGX,HY_MAGY,HY_MAGZ,HY_GAMC,HY_GAME,HY_GRAV)
#else
           case (HY_DENS,HY_GAMC,HY_GAME,HY_GRAV)
#endif
              delbar(k)=    mc(delV_L(k),delV_R(k))
           case (HY_VELX,HY_VELY,HY_VELZ,HY_PRES)
              delbar(k)=minmod(delV_L(k),delV_R(k))
           end select
        enddo
     endif

  elseif (hy_limiter==LIMITED) then
     if (hy_charLimiting) then
        do k=1,HY_WAVENUM
           del1=dot_product(leig(1:HY_VARINUM,k),delV_L(HY_DENS:HY_VARINUM))
           del2=dot_product(leig(1:HY_VARINUM,k),delV_R(HY_DENS:HY_VARINUM))
           if (del2 > 0.) then
              delbar(k) = max(0.,min(hy_LimitedSlopeBeta*del1,del2),&
                                 min(del1,hy_LimitedSlopeBeta*del2))
           else
              delbar(k) = min(0.,max(hy_LimitedSlopeBeta*del1,del2),&
                                 max(del1,hy_LimitedSlopeBeta*del2))
           endif
        enddo
!#ifndef FLASH_EOS_GAMMA
        ! for gamc, game, eint, and gravity
        do k=HY_GAMC,HY_END_VARS
           del1=delV_L(k)
           del2=delV_R(k)
           if (del2 > 0.) then
              delbar(k) = max(0.,min(hy_LimitedSlopeBeta*del1,del2),&
                                 min(del1,hy_LimitedSlopeBeta*del2))
           else
              delbar(k) = min(0.,max(hy_LimitedSlopeBeta*del1,del2),&
                                 max(del1,hy_LimitedSlopeBeta*del2))
           endif
        enddo
!#endif
     else
        do k=1,HY_VARINUM+extraVarinum
           del1=delV_L(k)
           del2=delV_R(k)
           if (del2 > 0.) then
              delbar(k) = max(0.,min(hy_LimitedSlopeBeta*del1,del2),&
                                 min(del1,hy_LimitedSlopeBeta*del2))
           else
              delbar(k) = min(0.,max(hy_LimitedSlopeBeta*del1,del2),&
                                 max(del1,hy_LimitedSlopeBeta*del2))
           endif
        enddo
     endif
  endif


!!$#ifdef FLASH_USM_MHD
!!$  if (.not. hy_charLimiting) then
!!$     ! No slope limiter for normal field in order to keep
!!$     ! continuous profile across cell boundaries.
!!$     select case(dir)
!!$     case(DIR_X)
!!$        delbar(HY_MAGX)=0.
!!$     case(DIR_Y)
!!$        delbar(HY_MAGY)=0.
!!$     case(DIR_Z)
!!$        delbar(HY_MAGZ)=0.
!!$     end select
!!$  endif
!!$#endif


End Subroutine hy_uhd_TVDslope
