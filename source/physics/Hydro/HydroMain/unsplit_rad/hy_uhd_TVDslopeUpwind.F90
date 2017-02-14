!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_TVDslopeUpwind
!!
!! NAME
!!
!!  hy_uhd_TVDslopeUpwind
!!
!! SYNOPSIS
!!
!!  hy_uhd_TVDslopeUpwind(integer(IN) :: dir,
!!                        real(IN)    :: VLL(HY_VARINUMMAX),
!!                        real(IN)    :: VL (HY_VARINUMMAX),
!!                        real(IN)    :: V0 (HY_VARINUMMAX),
!!                        real(IN)    :: VR (HY_VARINUMMAX),
!!                        real(IN)    :: VRR(HY_VARINUMMAX),
!!                        real(IN)    :: lambdaL(HY_WAVENUM),
!!                        real(IN)    :: lambda (HY_WAVENUM),
!!                        real(IN)    :: lambdaR(HY_WAVENUM),
!!                        real(IN)    :: leig(HY_VARINUM,HY_WAVENUM),
!!                        real(OUT)   :: delbar(HY_VARINUMMAX))
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
!!   VLL         - primitive variables at the second left  neighboring cell
!!   VL          - primitive variables at the first left  neighboring cell
!!   V0          - cell centered variables in primitive form
!!   VR          - primitive variables at the first right neighboring cell
!!   VRR         - primitive variables at the second right neighboring cell
!!   lambdaL     - left eigenvalue evaluated at the first left neighboring cell
!!   lambda      - left eigenvalue evaluated at the current cell
!!   lambdaR     - left eigenvalue evaluated at the first right neighboring cell
!!   leig        - left eigenvector evaluated at the current cell
!!   delbar      - limited slopes for primitive variables + gamc,game,gravity
!!
!!
!!
!! REFERENCES
!!
!!  * Lee, D., An Upwind Slope Limiter for PPM that Preserves Monotonicity in Magnetohydrodynamics,
!!     Astronum, 2010, a full paper in preparation.
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!  * Balsara, ApJ. Suppl, 151:149--184, 2004
!!
!!***


Subroutine hy_uhd_TVDslopeUpwind(dir,VLL,VL,V0,VR,VRR,lambdaL,lambda,lambdaR,leig,delbar)

  use Hydro_data,           ONLY : hy_limiter,hy_ContactSteepening
  use hy_uhd_slopeLimiters, ONLY : minmod,mc,vanLeer
  use Hydro_data,           ONLY : hy_charLimiting,hy_LimitedSlopeBeta

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Arguments type declaration ------------------------------------------
  integer, intent(IN) :: dir
  real, dimension(HY_VARINUMMAX),intent(IN)  :: VLL,VL,V0,VR,VRR
  real, dimension(HY_WAVENUM),   intent(IN)  :: lambdaL,lambda,lambdaR
  real, dimension(HY_VARINUM,HY_WAVENUM),intent(IN) :: leig
  real, dimension(HY_VARINUMMAX),intent(OUT) :: delbar
  !! ---------------------------------------------------------------------
  integer :: k,extraVarinum
  real, PARAMETER :: epsilon = 1.e-12
  real, dimension(HY_VARINUMMAX) :: delL1,delR1,del


  ! initialize with zeros
  delbar = 0.

!!$#ifndef GRAVITY
!!$  HY_END_VARS = HY_EINT
!!$#else
!!$  HY_END_VARS = HY_GRAV
!!$#endif
!!$
!!$#ifdef FLASH_UHD_3T  
!!$  HY_END_VARS = HY_ERAD
!!$#endif


  extraVarinum = HY_END_VARS-HY_GAMC+1


  delL1(HY_DENS:HY_VARINUM) = V0(HY_DENS:HY_VARINUM)-VL(HY_DENS:HY_VARINUM) !lambda >0
  delR1(HY_DENS:HY_VARINUM) = VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM) !lambda <0
  del=0.

  if (hy_limiter==MC) then

     if (hy_charLimiting) then !mc char limiting
        do k=1,HY_WAVENUM

           ! upwinding begins
           ! (1) normal velocity > 0.
           if (lambda(HY_ENTROPY) > epsilon) then

              delbar(k)=mc(dot_product(leig(1:HY_VARINUM,k), V0(HY_DENS:HY_VARINUM)- VL(HY_DENS:HY_VARINUM)),&
                           dot_product(leig(1:HY_VARINUM,k), VL(HY_DENS:HY_VARINUM)-VLL(HY_DENS:HY_VARINUM)))

              ! (1)-a: linearly degenerate wave steepening 
              select case(k)
#ifdef FLASH_USM_MHD
              case(HY_ALFNLEFT,HY_ENTROPY,HY_ALFNRGHT)
#else
              case(HY_ENTROPY)
#endif                 
                 del(HY_DENS:HY_VARINUM) = delR1(HY_DENS:HY_VARINUM)
                 delbar(k)=mc(delbar(k),dot_product(leig(1:HY_VARINUM,k), del(HY_DENS:HY_VARINUM) ))
              end select

           ! (2) normal velocity < 0.
           elseif (lambda(HY_ENTROPY) < -epsilon) then
              delbar(k)=mc(dot_product(leig(1:HY_VARINUM,k), VRR(HY_DENS:HY_VARINUM)-VR(HY_DENS:HY_VARINUM)),&
                           dot_product(leig(1:HY_VARINUM,k),  VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))

              ! (2)-a: linearly degenerate wave steepening 
              select case(k)
#ifdef FLASH_USM_MHD
              case(HY_ALFNLEFT,HY_ENTROPY,HY_ALFNRGHT)
#else
              case(HY_ENTROPY)
#endif
                 del(HY_DENS:HY_VARINUM) = delL1(HY_DENS:HY_VARINUM)
                 delbar(k)=mc(delbar(k),dot_product(leig(1:HY_VARINUM,k), del(HY_DENS:HY_VARINUM)))
              end select

           ! (3) normal velocity = 0.
           else ! lambda(HY_ENTROPY) = 0
              delbar(k)=mc(dot_product(leig(1:HY_VARINUM,k), delR1(HY_DENS:HY_VARINUM)),&
                           dot_product(leig(1:HY_VARINUM,k), delL1(HY_DENS:HY_VARINUM)))
           endif

        enddo


        ! for gamc, game, and gravity
        do k=HY_GAMC,HY_END_VARS
           !*** MC (monotonized  central) limiter ***
           if (lambda(HY_ENTROPY) > epsilon) then
              delbar(k)=mc(V0(k)- VL(k), VL(k)-VLL(k))

           elseif (lambda(HY_ENTROPY) < -epsilon) then
              delbar(k)=mc(VRR(k)-VR(k), VR(k)-V0(k))

           else
              delbar(k)=mc(VR(k)-V0(k), V0(k)-VL(k))
           endif
        enddo

        
     else ! mc primitive limiting
        do k=1,HY_VARINUM+extraVarinum
           if (lambda(HY_ENTROPY) > epsilon) then
              delbar(k)=mc(V0(k)- VL(k), VL(k)-VLL(k))

           elseif (lambda(HY_ENTROPY) < -epsilon) then
              delbar(k)=mc(VRR(k)-VR(k), VR(k)-V0(k))

           else
              delbar(k)=mc(VR(k)-V0(k), V0(k)-VL(k))
           endif
        enddo

     endif



  elseif (hy_limiter==MINMOD) then

     if (hy_charLimiting) then !minmod char limiting
        do k=1,HY_WAVENUM

           ! upwinding begins
           ! (1) normal velocity > 0.
           if (lambda(HY_ENTROPY) > epsilon) then

              delbar(k)=minmod(dot_product(leig(1:HY_VARINUM,k), V0(HY_DENS:HY_VARINUM)- VL(HY_DENS:HY_VARINUM)),&
                               dot_product(leig(1:HY_VARINUM,k), VL(HY_DENS:HY_VARINUM)-VLL(HY_DENS:HY_VARINUM)))

              ! (1)-a: linearly degenerate wave steepening 
              select case(k)
#ifdef FLASH_USM_MHD
              case(HY_ALFNLEFT,HY_ENTROPY,HY_ALFNRGHT)
#else
              case(HY_ENTROPY)
#endif                 
                 del(HY_DENS:HY_VARINUM) = delR1(HY_DENS:HY_VARINUM)
                 delbar(k)=minmod(delbar(k),dot_product(leig(1:HY_VARINUM,k), del(HY_DENS:HY_VARINUM) ))
              end select

           ! (2) normal velocity < 0.
           elseif (lambda(HY_ENTROPY) < -epsilon) then
              delbar(k)=minmod(dot_product(leig(1:HY_VARINUM,k), VRR(HY_DENS:HY_VARINUM)-VR(HY_DENS:HY_VARINUM)),&
                               dot_product(leig(1:HY_VARINUM,k),  VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))

              ! (2)-a: linearly degenerate wave steepening 
              select case(k)
#ifdef FLASH_USM_MHD
              case(HY_ALFNLEFT,HY_ENTROPY,HY_ALFNRGHT)
#else
              case(HY_ENTROPY)
#endif
                 del(HY_DENS:HY_VARINUM) = delL1(HY_DENS:HY_VARINUM)
                 delbar(k)=minmod(delbar(k),dot_product(leig(1:HY_VARINUM,k), del(HY_DENS:HY_VARINUM)))
              end select

           ! (3) normal velocity = 0.
           else ! lambda(HY_ENTROPY) = 0
              delbar(k)=minmod(dot_product(leig(1:HY_VARINUM,k), delR1(HY_DENS:HY_VARINUM)),&
                               dot_product(leig(1:HY_VARINUM,k), delL1(HY_DENS:HY_VARINUM)))
           endif

        enddo


        ! for gamc, game, and gravity
        do k=HY_GAMC,HY_END_VARS
           !*** MC (monotonized  central) limiter ***
           if (lambda(HY_ENTROPY) > epsilon) then
              delbar(k)=minmod(V0(k)- VL(k), VL(k)-VLL(k))

           elseif (lambda(HY_ENTROPY) < -epsilon) then
              delbar(k)=minmod(VRR(k)-VR(k), VR(k)-V0(k))

           else
              delbar(k)=minmod(VR(k)-V0(k), V0(k)-VL(k))

           endif
        enddo

        
     else ! minmod primitive limiting
        do k=1,HY_VARINUM+extraVarinum
           if (lambda(HY_ENTROPY) > epsilon) then
              delbar(k)=minmod(V0(k)- VL(k), VL(k)-VLL(k))

           elseif (lambda(HY_ENTROPY) < -epsilon) then
              delbar(k)=minmod(VRR(k)-VR(k), VR(k)-V0(k))

           else
              delbar(k)=minmod(VR(k)-V0(k), V0(k)-VL(k))
           endif
        enddo
     endif


  elseif (hy_limiter==VANLEER) then

     if (hy_charLimiting) then !vanLeer char limiting
        do k=1,HY_WAVENUM

           ! upwinding begins
           ! (1) normal velocity > 0.
           if (lambda(HY_ENTROPY) > epsilon) then

              delbar(k)=vanLeer(dot_product(leig(1:HY_VARINUM,k), V0(HY_DENS:HY_VARINUM)- VL(HY_DENS:HY_VARINUM)),&
                                dot_product(leig(1:HY_VARINUM,k), VL(HY_DENS:HY_VARINUM)-VLL(HY_DENS:HY_VARINUM)))

              ! (1)-a: linearly degenerate wave steepening 
              select case(k)
#ifdef FLASH_USM_MHD
              case(HY_ALFNLEFT,HY_ENTROPY,HY_ALFNRGHT)
#else
              case(HY_ENTROPY)
#endif                 
                 del(HY_DENS:HY_VARINUM) = delR1(HY_DENS:HY_VARINUM)
                 delbar(k)=vanLeer(delbar(k),dot_product(leig(1:HY_VARINUM,k), del(HY_DENS:HY_VARINUM) ))
              end select

           ! (2) normal velocity < 0.
           elseif (lambda(HY_ENTROPY) < -epsilon) then
              delbar(k)=vanLeer(dot_product(leig(1:HY_VARINUM,k), VRR(HY_DENS:HY_VARINUM)-VR(HY_DENS:HY_VARINUM)),&
                                dot_product(leig(1:HY_VARINUM,k),  VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))

              ! (2)-a: linearly degenerate wave steepening 
              select case(k)
#ifdef FLASH_USM_MHD
              case(HY_ALFNLEFT,HY_ENTROPY,HY_ALFNRGHT)
#else
              case(HY_ENTROPY)
#endif
                 del(HY_DENS:HY_VARINUM) = delL1(HY_DENS:HY_VARINUM)
                 delbar(k)=vanLeer(delbar(k),dot_product(leig(1:HY_VARINUM,k), del(HY_DENS:HY_VARINUM)))
              end select

           ! (3) normal velocity = 0.
           else ! lambda(HY_ENTROPY) = 0
              delbar(k)=vanLeer(dot_product(leig(1:HY_VARINUM,k), delR1(HY_DENS:HY_VARINUM)),&
                                dot_product(leig(1:HY_VARINUM,k), delL1(HY_DENS:HY_VARINUM)))
           endif

        enddo


        ! for gamc, game, and gravity
        do k=HY_GAMC,HY_END_VARS
           !*** MC (monotonized  central) limiter ***
           if (lambda(HY_ENTROPY) > epsilon) then
              delbar(k)=vanLeer(V0(k)- VL(k), VL(k)-VLL(k))

           elseif (lambda(HY_ENTROPY) < -epsilon) then
              delbar(k)=vanLeer(VRR(k)-VR(k), VR(k)-V0(k))

           else
              delbar(k)=vanLeer(VR(k)-V0(k), V0(k)-VL(k))
           endif
        enddo

        
     else ! vanLeer primitive limiting
        do k=1,HY_VARINUM+extraVarinum
           if (lambda(HY_ENTROPY) > epsilon) then
              delbar(k)=vanLeer(V0(k)- VL(k), VL(k)-VLL(k))

           elseif (lambda(HY_ENTROPY) < -epsilon) then
              delbar(k)=vanLeer(VRR(k)-VR(k), VR(k)-V0(k))

           else
              delbar(k)=vanLeer(VR(k)-V0(k), V0(k)-VL(k))
           endif
        enddo
     endif
 



  elseif (hy_limiter==HYBRID) then

     if (hy_charLimiting) then

        ! (1) normal velocity > 0.
        if (lambda(HY_ENTROPY) > epsilon) then

           delbar(HY_FASTLEFT)=&
                minmod(dot_product(leig(1:HY_VARINUM,HY_FASTLEFT),V0(HY_DENS:HY_VARINUM)- VL(HY_DENS:HY_VARINUM)),&
                       dot_product(leig(1:HY_VARINUM,HY_FASTLEFT),VL(HY_DENS:HY_VARINUM)-VLL(HY_DENS:HY_VARINUM)))

#ifdef FLASH_USM_MHD
           delbar(HY_ALFNLEFT)=&
                mc(dot_product(leig(1:HY_VARINUM,HY_ALFNLEFT),V0(HY_DENS:HY_VARINUM)-VL(HY_DENS:HY_VARINUM)),&
                   dot_product(leig(1:HY_VARINUM,HY_ALFNLEFT),VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))

#endif
           delbar(HY_SLOWLEFT)=&
                minmod(dot_product(leig(1:HY_VARINUM,HY_SLOWLEFT),V0(HY_DENS:HY_VARINUM)- VL(HY_DENS:HY_VARINUM)),&
                       dot_product(leig(1:HY_VARINUM,HY_SLOWLEFT),VL(HY_DENS:HY_VARINUM)-VLL(HY_DENS:HY_VARINUM)))

           delbar(HY_ENTROPY)=&
                mc(dot_product(leig(1:HY_VARINUM,HY_ENTROPY), V0(HY_DENS:HY_VARINUM)-VL(HY_DENS:HY_VARINUM)),&
                   dot_product(leig(1:HY_VARINUM,HY_ENTROPY), VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))

           delbar(HY_SLOWRGHT)=&
                minmod(dot_product(leig(1:HY_VARINUM,HY_SLOWRGHT),V0(HY_DENS:HY_VARINUM)- VL(HY_DENS:HY_VARINUM)),&
                       dot_product(leig(1:HY_VARINUM,HY_SLOWRGHT),VL(HY_DENS:HY_VARINUM)-VLL(HY_DENS:HY_VARINUM)))

#ifdef FLASH_USM_MHD
           delbar(HY_ALFNRGHT)=&
                mc(dot_product(leig(1:HY_VARINUM,HY_ALFNRGHT),V0(HY_DENS:HY_VARINUM)-VL(HY_DENS:HY_VARINUM)),&
                   dot_product(leig(1:HY_VARINUM,HY_ALFNRGHT),VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))
#endif
           delbar(HY_FASTRGHT)=&
                minmod(dot_product(leig(1:HY_VARINUM,HY_FASTRGHT),V0(HY_DENS:HY_VARINUM)- VL(HY_DENS:HY_VARINUM)),&
                       dot_product(leig(1:HY_VARINUM,HY_FASTRGHT),VL(HY_DENS:HY_VARINUM)-VLL(HY_DENS:HY_VARINUM)))


           ! (2) normal velocity < 0.
        elseif (lambda(HY_ENTROPY) < -epsilon) then
           delbar(HY_FASTLEFT)=&
                minmod(dot_product(leig(1:HY_VARINUM,HY_FASTLEFT),VRR(HY_DENS:HY_VARINUM)-VR(HY_DENS:HY_VARINUM)),&
                       dot_product(leig(1:HY_VARINUM,HY_FASTLEFT),VR(HY_DENS:HY_VARINUM)- V0(HY_DENS:HY_VARINUM)))
#ifdef FLASH_USM_MHD
           delbar(HY_ALFNLEFT)=&
                mc(dot_product(leig(1:HY_VARINUM,HY_ALFNLEFT),V0(HY_DENS:HY_VARINUM)-VL(HY_DENS:HY_VARINUM)),&
                   dot_product(leig(1:HY_VARINUM,HY_ALFNLEFT),VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))
#endif
           delbar(HY_SLOWLEFT)=&
                minmod(dot_product(leig(1:HY_VARINUM,HY_SLOWLEFT),VRR(HY_DENS:HY_VARINUM)-VR(HY_DENS:HY_VARINUM)),&
                       dot_product(leig(1:HY_VARINUM,HY_SLOWLEFT),VR(HY_DENS:HY_VARINUM)- V0(HY_DENS:HY_VARINUM)))

           delbar(HY_ENTROPY)=&
                mc(dot_product(leig(1:HY_VARINUM,HY_ENTROPY), V0(HY_DENS:HY_VARINUM)-VL(HY_DENS:HY_VARINUM)),&
                   dot_product(leig(1:HY_VARINUM,HY_ENTROPY), VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))

           delbar(HY_SLOWRGHT)=&
                minmod(dot_product(leig(1:HY_VARINUM,HY_SLOWRGHT),VRR(HY_DENS:HY_VARINUM)-VR(HY_DENS:HY_VARINUM)),&
                       dot_product(leig(1:HY_VARINUM,HY_SLOWRGHT),VR(HY_DENS:HY_VARINUM)- V0(HY_DENS:HY_VARINUM)))
#ifdef FLASH_USM_MHD
           delbar(HY_ALFNRGHT)=&
                mc(dot_product(leig(1:HY_VARINUM,HY_ALFNRGHT),V0(HY_DENS:HY_VARINUM)-VL(HY_DENS:HY_VARINUM)),&
                   dot_product(leig(1:HY_VARINUM,HY_ALFNRGHT),VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))
#endif
           delbar(HY_FASTRGHT)=&
               minmod(dot_product(leig(1:HY_VARINUM,HY_FASTRGHT),VRR(HY_DENS:HY_VARINUM)-VR(HY_DENS:HY_VARINUM)),&
                      dot_product(leig(1:HY_VARINUM,HY_FASTRGHT),VR(HY_DENS:HY_VARINUM) -V0(HY_DENS:HY_VARINUM)))


        else
           ! (3) normal velocity = 0.
           ! No upwinding for wave strength = 0.
           delbar(HY_FASTLEFT)=&
                minmod(dot_product(leig(1:HY_VARINUM,HY_FASTLEFT),V0(HY_DENS:HY_VARINUM)-VL(HY_DENS:HY_VARINUM)),&
                       dot_product(leig(1:HY_VARINUM,HY_FASTLEFT),VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))
#ifdef FLASH_USM_MHD
           delbar(HY_ALFNLEFT)=&
                mc(dot_product(leig(1:HY_VARINUM,HY_ALFNLEFT),V0(HY_DENS:HY_VARINUM)-VL(HY_DENS:HY_VARINUM)),&
                   dot_product(leig(1:HY_VARINUM,HY_ALFNLEFT),VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))
#endif
           delbar(HY_SLOWLEFT)=&
                minmod(dot_product(leig(1:HY_VARINUM,HY_SLOWLEFT),V0(HY_DENS:HY_VARINUM)-VL(HY_DENS:HY_VARINUM)),&
                       dot_product(leig(1:HY_VARINUM,HY_SLOWLEFT),VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))

           delbar(HY_ENTROPY)=&
                mc(dot_product(leig(1:HY_VARINUM,HY_ENTROPY), V0(HY_DENS:HY_VARINUM)-VL(HY_DENS:HY_VARINUM)),&
                   dot_product(leig(1:HY_VARINUM,HY_ENTROPY), VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))

           delbar(HY_SLOWRGHT)=&
                minmod(dot_product(leig(1:HY_VARINUM,HY_SLOWRGHT),V0(HY_DENS:HY_VARINUM)-VL(HY_DENS:HY_VARINUM)),&
                       dot_product(leig(1:HY_VARINUM,HY_SLOWRGHT),VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))
#ifdef FLASH_USM_MHD
           delbar(HY_ALFNRGHT)=&
                mc(dot_product(leig(1:HY_VARINUM,HY_ALFNRGHT),V0(HY_DENS:HY_VARINUM)-VL(HY_DENS:HY_VARINUM)),&
                   dot_product(leig(1:HY_VARINUM,HY_ALFNRGHT),VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))
#endif
           delbar(HY_FASTRGHT)=&
                minmod(dot_product(leig(1:HY_VARINUM,HY_FASTRGHT),V0(HY_DENS:HY_VARINUM)-VL(HY_DENS:HY_VARINUM)),&
                       dot_product(leig(1:HY_VARINUM,HY_FASTRGHT),VR(HY_DENS:HY_VARINUM)-V0(HY_DENS:HY_VARINUM)))
        endif



        ! for gamc, game, and gravity
        do k=HY_GAMC,HY_END_VARS
           !*** MC (monotonized  central) limiter ***
           if (lambda(HY_ENTROPY) > epsilon) then
              delbar(k)=mc(V0(k)- VL(k), VL(k)-VLL(k))

           elseif (lambda(HY_ENTROPY) < -epsilon) then
              delbar(k)=mc(VRR(k)-VR(k), VR(k)-V0(k))

           else
              delbar(k)=mc(VR(k)-V0(k), V0(k)-VL(k))

           endif
        enddo


     else ! mc primitive limiting
        do k=1,HY_VARINUM+extraVarinum
           select case (k)
#ifdef FLASH_USM_MHD
           case (HY_DENS,HY_MAGX,HY_MAGY,HY_MAGZ,HY_GAMC,HY_GAME,HY_GRAV)
#else
           case (HY_DENS,HY_GAMC,HY_GAME,HY_GRAV)
#endif
              if (lambda(HY_ENTROPY) > 0.) then
                 delbar(k)=    mc(V0(k)-VL(k),VL(k)-VLL(k))
              elseif (lambda(HY_ENTROPY) < 0.) then
                 delbar(k)=    mc(VRR(k)-VR(k),VR(k)-V0(k))
              elseif (lambda(HY_ENTROPY) == 0.) then
                 delbar(k)=    mc(V0(k)-VL(k),VR(k)-V0(k))
              endif
           case (HY_VELX,HY_VELY,HY_VELZ,HY_PRES)
              if (lambda(HY_ENTROPY) > 0.) then
                 delbar(k)=    minmod(V0(k)-VL(k),VL(k)-VLL(k))
              elseif (lambda(HY_ENTROPY) < 0.) then
                 delbar(k)=    minmod(VRR(k)-VR(k),VR(k)-V0(k))
              elseif (lambda(HY_ENTROPY) == 0.) then
                 delbar(k)=    minmod(V0(k)-VL(k),VR(k)-V0(k))
              endif
           end select
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



End Subroutine hy_uhd_TVDslopeUpwind
