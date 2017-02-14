!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/grv_bhAccLong1P2I
!!
!! NAME
!!
!!  grv_bhAccLong1P2I
!!
!!
!! SYNOPSIS
!!
!!   real(0:12) ef = grv_bhAccLong1P2I(
!!                       real(in) :: Linv,
!!                       real(in) :: ewald_dzeta,
!!                       real(in) :: ewald_eta,
!!                       integer(in) :: hi,
!!                       real(in) :: x,
!!                       real(in) :: y,
!!                       real(in) :: z
!!                   )
!!
!! DESCRIPTION
!!
!!   Calculates potential, acceleration and partial derivatives of acceleration 
!!   for the Ewald field for long-range contributions in the case of mixed boundary 
!!   conditions - periodic in 1 direction and isolated in the other two.
!!
!! ARGUMENTS
!!
!!   Linv         : inverted size of the computational domain in the periodic
!!                  direction
!!   ewald_dzeta  : coefficient dzeta needed to calculate the Ewald field
!!   ewald_eta    : coefficient eta needed to calculate the Ewald field
!!   hi           : wavenumber in periodic direction
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

function grv_bhAccLong1P2I(Linv, ewald_dzeta, ewald_eta, hi, x, y, z)
  implicit none
#include "constants.h"

  real, intent(IN) :: Linv, ewald_dzeta, ewald_eta, x, y, z
  integer, intent(IN) :: hi
  real :: grv_bhAccLong1P2I(0:12)

  real, parameter :: pi = PI

! parameters controlling numerical integration
  real, parameter :: int_xmin = 1.0e-8
  real, parameter :: int_xmax = 5.0
  integer, parameter :: int_n = 300

  real :: Linv2, Linv_acc, Linv_acc_der
  real :: sin_lx, cos_lx, exp_nh2, K_integral, M_integral, N_integral
! called function
  real :: grv_IntSimpson
! functions parsed as arguments
  real, external :: grv_Coef1P2I, grv_Coef1P2I_der, grv_Coef1P2I_der2

  Linv2 = 2.0*Linv
  Linv_acc = pi*(Linv2**2)
  Linv_acc_der = (pi**2)*Linv2**3

  sin_lx = sin(Linv2*pi*hi*x)
  cos_lx = cos(Linv2*pi*hi*x)
  exp_nh2 = exp(-ewald_dzeta*(hi**2))
  K_integral = grv_IntSimpson(grv_Coef1P2I,hi, &
  &            ewald_eta,ewald_dzeta,int_xmin,int_xmax,int_n)
  M_integral = grv_IntSimpson(grv_Coef1P2I_der,hi,ewald_eta, &
  &            ewald_dzeta,int_xmin,int_xmax,int_n)
  N_integral = grv_IntSimpson(grv_Coef1P2I_der2,hi,ewald_eta, &
  &            ewald_dzeta,int_xmin,int_xmax,int_n)


  ! potential
  grv_bhAccLong1P2I(0) = Linv2*K_integral*exp_nh2*cos_lx
  ! acceleration
  grv_bhAccLong1P2I(1) = Linv_acc*hi*exp_nh2*sin_lx*K_integral
  grv_bhAccLong1P2I(2) = Linv_acc*y*exp_nh2*cos_lx*M_integral/sqrt(y**2+z**2+1.0d-99)
  grv_bhAccLong1P2I(3) = Linv_acc*z*exp_nh2*cos_lx*M_integral/sqrt(y**2+z**2+1.0d-99)
  ! partial derivatives of acceleration
  grv_bhAccLong1P2I(4) = Linv_acc_der*(hi**2)*exp_nh2*cos_lx*K_integral
  grv_bhAccLong1P2I(5) = -Linv_acc_der*hi*exp_nh2*sin_lx*M_integral*y/sqrt(y**2+z**2+1.0d-99)
  grv_bhAccLong1P2I(6) = -Linv_acc_der*hi*exp_nh2*sin_lx*M_integral*z/sqrt(y**2+z**2+1.0d-99)
  grv_bhAccLong1P2I(8) = Linv_acc*exp_nh2*cos_lx*M_integral*z**2/(sqrt(y**2+z**2+1.0d-99)**3) + &
  &     Linv_acc_der*exp_nh2*cos_lx*N_integral*y**2/(y**2+z**2+1.0d-99)
  grv_bhAccLong1P2I(9) = -Linv_acc*exp_nh2*cos_lx*M_integral*y*z/(sqrt(y**2+z**2+1.0d-99)**3) + &
  &     Linv_acc_der*exp_nh2*cos_lx*N_integral*y*z/(y**2+z**2+1.0d-99)
  grv_bhAccLong1P2I(12) = Linv_acc*exp_nh2*cos_lx*M_integral*y**2/(sqrt(y**2+z**2+1.0d-99)**3) + &
  &     Linv_acc_der*exp_nh2*cos_lx*N_integral*z**2/(y**2+z**2+1.0d-99)

  ! in the case of ewald_eta==0 in N_integral is division by zero - but final result is zero
  if (ewald_eta.lt.tiny(ewald_eta)) then
    grv_bhAccLong1P2I(8) = 0.0D0
    grv_bhAccLong1P2I(9) = 0.0D0
    grv_bhAccLong1P2I(12) = 0.0D0
  endif
end

