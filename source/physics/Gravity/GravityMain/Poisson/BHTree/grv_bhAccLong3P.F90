!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/grv_bhAccLong3P
!!
!! NAME
!!
!!  grv_bhAccLong3P
!!
!!
!! SYNOPSIS
!!
!!   real(0:12) ef = grv_bhAccLong3P(
!!                       real(in) :: Linv,
!!                       real(in) :: ewald_dzeta,
!!                       real(in) :: ratio_pinv1,
!!                       real(in) :: ratio_pinv2,
!!                       integer(in) :: hi,
!!                       integer(in) :: hj,
!!                       integer(in) :: hk,
!!                       real(in) :: x,
!!                       real(in) :: y,
!!                       real(in) :: z
!!                   )
!!
!! DESCRIPTION
!!
!!   Calculates potential, acceleration and partial derivatives of acceleration 
!!   for the Ewald field for long-range contributions in the case of periodic
!!   boundary conditions.
!!
!! ARGUMENTS
!!
!!   Linv         : inverted size of the computational domain in one 
!!                  of three periodic directions
!!   ewald_dzeta  : coefficient dzeta needed to calculate the Ewald field
!!   ratio_pinv1  : ratio between sizes of computational domain in one periodic direction
!!   ratio_pinv2  : ratio between sizes of computational domain in the other periodic direction
!!   hi           : wavenumber in the direction x
!!   hj           : wavenumber in the direction y
!!   hk           : wavenumber in the direction z
!!   x            : x-coordinate of the point where the Ewald field is determined
!!   y            : y-coordinate of the point where the Ewald field is determined
!!   z            : z-coordinate of the point where the Ewald field is determined
!!
!! RESULT
!!
!!   Returns an array with 13 real numbers. 
!!   Index 0 is a value of potential
!!   indices 1 to 3 are values of acceleration
!!   indices 4 to 12 are values of partial derivatives of acceleration
!!   Note that three derivatives (at 7, 10 and 11) are not evaluated because 
!!   of the symmetry.
!!   
!!
!!***

function grv_bhAccLong3P(Linv, ewald_dzeta, ratio_pinv1, ratio_pinv2, hi, hj, hk, x, y, z)

! computes long range interactions in the case 3p

  implicit none

#include "constants.h"

  real, parameter :: pi = PI

  real, intent(IN) :: Linv, ewald_dzeta, ratio_pinv1, ratio_pinv2, x, y, z
  integer, intent(IN) :: hi, hj, hk
  real :: grv_bhAccLong3P(0:12)

  real :: sin_lx, cos_lx, k2_inv, exp_nh2
  real :: Linv2, Linv_acc, Linv_acc_der

  if ((hi**2+hj**2+hk**2).eq.0) then 
    grv_bhAccLong3P(:) = 0.0D0
    return
  endif

  Linv2 = Linv*ratio_pinv1*ratio_pinv2/pi
  Linv_acc = 2.0*(Linv**2)*ratio_pinv1*ratio_pinv2
  Linv_acc_der = 4.0*pi*(Linv**3)*ratio_pinv1*ratio_pinv2
  sin_lx = sin(2*pi*Linv*(hi*x+hj*y*ratio_pinv1+hk*z*ratio_pinv2))
  cos_lx = cos(2*pi*Linv*(hi*x+hj*y*ratio_pinv1+hk*z*ratio_pinv2))
  k2_inv = 1.0D0/(hi**2+(hj*ratio_pinv1)**2+(hk*ratio_pinv2)**2)
  exp_nh2 = exp(-ewald_dzeta*(hi**2+(hj*ratio_pinv1)**2+(hk*ratio_pinv2)**2))


! potential
  grv_bhAccLong3P(0) = Linv2*exp_nh2*cos_lx*k2_inv

! acceleration
  grv_bhAccLong3P(1) = Linv_acc*exp_nh2*hi*sin_lx*k2_inv
  grv_bhAccLong3P(2) = ratio_pinv1*Linv_acc*exp_nh2* hj*sin_lx*k2_inv
  grv_bhAccLong3P(3) = ratio_pinv2*Linv_acc*exp_nh2* hk*sin_lx*k2_inv
! partial derivatives of acceleration
  grv_bhAccLong3P(4) = Linv_acc_der*exp_nh2 &
  &                * hi * hi * cos_lx *k2_inv
  grv_bhAccLong3P(5) = Linv_acc_der*ratio_pinv1*exp_nh2 &
  &                * hi * hj * cos_lx *k2_inv
  grv_bhAccLong3P(6) = Linv_acc_der*ratio_pinv2*exp_nh2 &
  &                * hi * hk * cos_lx *k2_inv
  grv_bhAccLong3P(8) = Linv_acc_der*(ratio_pinv1**2)*exp_nh2 &
  &                * hj * hj * cos_lx *k2_inv
  grv_bhAccLong3P(9) = Linv_acc_der*(ratio_pinv1*ratio_pinv2)*exp_nh2 &
  &                * hj * hk * cos_lx *k2_inv
  grv_bhAccLong3P(12) = Linv_acc_der*(ratio_pinv2**2)*exp_nh2 &
  &                * hk * hk * cos_lx *k2_inv

end

