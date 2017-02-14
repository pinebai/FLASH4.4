!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_entropyFix
!!
!! NAME
!!
!!  hy_uhd_entropyFix
!!
!!
!! SYNOPSIS
!!
!!  hy_uhd_entropyFix ( real(INOUT) :: lambda (HY_WAVENUM),
!!                      real(IN)    :: lambdaL(HY_WAVENUM),
!!                      real(IN)    :: lambdaR(HY_WAVENUM))
!!
!! ARGUMENTS
!!
!!   lambda  - an array for eigenvalues of the average state
!!   lambdaL - an array for eigenvalues of the left state
!!   lambdaR - an array for eigenvalues of the right state
!!
!! DESCRIPTION
!! 
!!   The entropy fix is required for an approximate Roe-type Riemann solver
!!   to keep waves away from violating the entropy condition.
!!   When the flow is expanding, this will keep the magnetosonic wave 
!!   speeds away from zero.
!!
!!   Two different entropy fix methods by (i) Harten and (ii) Harten & Hymann (default)
!!   are provided. Users can also choose different method by setting, for example,
!!
!!        entropyFixMethod = "HARTEN" 
!!
!!   in flash.par. 
!!
!! REFERENCES
!!
!!   * Harten, JCP, 49:357-393, 1983
!!   * Harten and Hyman, JCP, 50:235-269
!!
!!***

Subroutine hy_uhd_entropyFix(lambda,lambdaL,lambdaR)

  use Hydro_data, ONLY : hy_entropyFixMethod

  implicit none

#include "Flash.h"
#include "UHD.h"

  !! Arguments type declaration -----------------------
  real, dimension(HY_WAVENUM), intent(INOUT) :: lambda
  real, dimension(HY_WAVENUM), intent(IN)    :: lambdaL,lambdaR
  !! --------------------------------------------------

  integer :: i
  real, dimension(HY_WAVENUM) :: dlambda
  real    :: L,R

  do i=1,HY_WAVENUM

     select case(i)
     ! Apply entropy fix only to genuinely nonlinear waves

#ifdef FLASH_USM_MHD
     case(HY_FASTLEFT,HY_SLOWLEFT,HY_SLOWRGHT,HY_FASTRGHT)
#else
     case(HY_FASTLEFT,HY_FASTRGHT)
#endif
        if ((lambdaR(i) > 0.) .and. (lambdaL(i) < 0.)) then
           if (hy_entropyFixMethod == HARTENHYMAN) then
              ! Harten-Hyman's entropy fix
              L = lambdaL(i)*(lambdaR(i)-lambda(i) )/(lambdaR(i)-lambdaL(i))
              R = lambdaR(i)*(lambda(i) -lambdaL(i))/(lambdaR(i)-lambdaL(i))
              lambda(i) = R-L

           elseif (hy_entropyFixMethod == HARTEN) then
              ! Harten's entropy fix
              dlambda(i)=max(4.0*(lambdaR(i)-lambdaL(i)),0.0 )
              if (abs(lambda(i)) < 0.5*dlambda(i)) then
                 if (lambda(i) >= 0.0) then
                    lambda(i) = (lambda(i)**2/dlambda(i)+.25 *dlambda(i))
                 else
                    lambda(i) =-(lambda(i)**2/dlambda(i)+.25 *dlambda(i))
                 endif
              endif
           endif
        endif

     end select

  enddo

End Subroutine hy_uhd_entropyFix
