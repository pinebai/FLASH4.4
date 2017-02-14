!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/grv_elintF
!!
!! NAME
!!
!!  grv_elintF
!!
!!
!! SYNOPSIS
!!
!!   real = grv_elintF(
!!            real(IN) :: theta,
!!            real(IN) :: k
!!   )
!!
!! DESCRIPTION
!!
!!   Calculates elliptic integral of the first kind F(theta,k).
!!
!! ARGUMENTS
!!
!!  theta - Jacobi amplitude
!!  k     - elliptic modulus (eccentricity)
!!
!! RESULT
!!
!!  Returns value of the elliptic integral F(theta,k).
!!
!!***

real function grv_elintF(theta, k)
  implicit none
#include "constants.h"
  real, intent(in) :: theta, k
  real :: fe
  integer :: n
  real :: a, b, c, d, g, a0, b0, d0, fac, ck, r, pi_half

  pi_half = 0.5*PI
  a0 = 1.0
  b0 = sqrt(1.0-k*k)
  d0 = theta
  r  = k*k
  if (k .ge. 1.0) then
     fe=dlog((1.0+sin(theta))/(cos(theta)+1d-99))
  else
     fac=1.0
     do n = 1,40
        a = 0.5*(a0+b0)
        b = sqrt(a0*b0)
        c = 0.5*(a0-b0)
        fac = 2.0*fac
        r = r + fac*c*c
        if (abs(d0 - pi_half) > 1d-99) then
           d = d0+atan((b0/a0)*tan(d0))
           g = g+c*sin(d)
           d0 = d+PI*int(d/PI+.5)
        endif
        a0 = a
        b0 = b
        !print *, "F: ", a0, b0, d0, a, b, c, d, k, r, fac, theta, PI
        if (c < 1.0e-7) go to 15
     enddo
15   ck = PI/(2.0*a)
     if (abs(d0 - pi_half) < 1d-99) then
        fe = ck
     else
        fe = d/(fac*a)
     endif
  endif
  grv_elintF = fe

end function grv_elintF


