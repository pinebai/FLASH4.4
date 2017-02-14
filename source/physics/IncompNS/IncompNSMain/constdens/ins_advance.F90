SUBROUTINE ins_predictor(uni,vni,wni,unew,vnew,wnew,uold,vold,&
        wold,p,dt,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,gama,rhoa,alfa)

      ! This routine computes the intermediate velocities based on
      ! the explicit second-order Adams-Bashforth scheme (gamma=1.5,
      ! rhoa=-0.5,alfa=1), or a third order Runge-Kutta method.

      use IncompNS_data, ONLY : ins_prescoeff,ins_dpdx,ins_dpdy,ins_dpdz, &
                                ins_gravX,ins_gravY,ins_gravZ

      implicit none
#include "Flash.h"
#include "constants.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dt,dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: unew,vnew,wnew,&
                                            uold,vold,wold,&
                                            p
      REAL, DIMENSION(:,:,:), INTENT(INOUT) :: uni,vni,wni

      REAL :: gama,rhoa,alfa

      uni(ix1:ix2+1,jy1:jy2,kz1:kz2) = uni(ix1:ix2+1,jy1:jy2,kz1:kz2) + &
         dt*(gama*unew(ix1:ix2+1,jy1:jy2,kz1:kz2)                     + &
             rhoa*uold(ix1:ix2+1,jy1:jy2,kz1:kz2)                     - &
             ins_prescoeff*alfa*(( p(ix1:ix2+1,jy1:jy2,kz1:kz2) -       &
                                   p(ix1-1:ix2,jy1:jy2,kz1:kz2) )/dx) - &
                           alfa*ins_dpdx + alfa*ins_gravX)

      vni(ix1:ix2,jy1:jy2+1,kz1:kz2) = vni(ix1:ix2,jy1:jy2+1,kz1:kz2) + &
         dt*(gama*vnew(ix1:ix2,jy1:jy2+1,kz1:kz2)                     + &
             rhoa*vold(ix1:ix2,jy1:jy2+1,kz1:kz2)                     - &
             ins_prescoeff*alfa*(( p(ix1:ix2,jy1:jy2+1,kz1:kz2) -       &
                                   p(ix1:ix2,jy1-1:jy2,kz1:kz2) )/dy) - &
                           alfa*ins_dpdy + alfa*ins_gravY)    

#if NDIM == MDIM
      wni(ix1:ix2,jy1:jy2,kz1:kz2+1) = wni(ix1:ix2,jy1:jy2,kz1:kz2+1) + &
         dt*(gama*wnew(ix1:ix2,jy1:jy2,kz1:kz2+1)                     + &
             rhoa*wold(ix1:ix2,jy1:jy2,kz1:kz2+1)                     - &
             ins_prescoeff*alfa*(( p(ix1:ix2,jy1:jy2,kz1:kz2+1) -       &
                                   p(ix1:ix2,jy1:jy2,kz1-1:kz2) )/dz) - &
                           alfa*ins_dpdz + alfa*ins_gravZ)
#endif

END SUBROUTINE ins_predictor

!########################################################################

SUBROUTINE ins_divergence(uni,vni,wni,ix1,ix2,jy1,jy2,kz1,kz2,&
         dx,dy,dz,divv)

      ! This routine computes the divergence of the velocity field.

      implicit none

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: uni,vni,wni
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: divv


      divv(ix1:ix2,jy1:jy2,kz1:kz2) =           & 
         ( uni(ix1+1:ix2+1,jy1:jy2,kz1:kz2) -   &
           uni(ix1:ix2,jy1:jy2,kz1:kz2) )/dx +  &       
         ( vni(ix1:ix2,jy1+1:jy2+1,kz1:kz2) -   &
           vni(ix1:ix2,jy1:jy2,kz1:kz2) )/dy


#if NDIM == MDIM
      divv(ix1:ix2,jy1:jy2,kz1:kz2) =           & 
           divv(ix1:ix2,jy1:jy2,kz1:kz2)     +  &
         ( wni(ix1:ix2,jy1:jy2,kz1+1:kz2+1) -   &
           wni(ix1:ix2,jy1:jy2,kz1:kz2) )/dz 
#endif

END SUBROUTINE ins_divergence

!########################################################################

SUBROUTINE ins_corrector(uni,vni,wni,p,ix1,ix2,jy1,jy2,kz1,kz2, &
        dt,dx,dy,dz,alfa)

      ! This routine computes the corrected divergence-free velocities.
    
      implicit none

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dt,dx,dy,dz,alfa
      REAL, DIMENSION(:,:,:), INTENT(IN) :: p
      REAL, DIMENSION(:,:,:), INTENT(INOUT) :: uni,vni,wni

      REAL :: coef

      coef = dt*alfa

      uni(ix1+1:ix2,jy1:jy2,kz1:kz2) =             & 
         uni(ix1+1:ix2,jy1:jy2,kz1:kz2) -          &
         coef*( p(ix1+1:ix2,jy1:jy2,kz1:kz2) -     &
                p(ix1:ix2-1,jy1:jy2,kz1:kz2) )/dx

      vni(ix1:ix2,jy1+1:jy2,kz1:kz2) =             &
         vni(ix1:ix2,jy1+1:jy2,kz1:kz2) -          &
         coef*( p(ix1:ix2,jy1+1:jy2,kz1:kz2) -     &
                p(ix1:ix2,jy1:jy2-1,kz1:kz2) )/dy

#if NDIM == MDIM
      wni(ix1:ix2,jy1:jy2,kz1+1:kz2) =             &
         wni(ix1:ix2,jy1:jy2,kz1+1:kz2) -          &
         coef*( p(ix1:ix2,jy1:jy2,kz1+1:kz2) -     &
                p(ix1:ix2,jy1:jy2,kz1:kz2-1) )/dz
#endif

END SUBROUTINE ins_corrector

