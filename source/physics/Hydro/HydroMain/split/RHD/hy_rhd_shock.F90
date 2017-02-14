!!****if* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_shock
!!
!! NAME
!!
!!   hy_rhd_shock
!!
!! SYNOPSIS
!!   hy_rhd_shock(real    (IN) :: tau0 
!!                real    (IN) :: u0,
!!                real    (IN) :: p0,
!!                real    (IN) :: g0,
!!                real    (IN) :: V0,
!!                real    (IN) :: h0,
!!                real    (IN) :: p1,
!!                real    (OUT):: u1,
!!                real    (OUT):: dudp,
!!                real    (OUT):: zeta,
!!                integer  (IN):: istate)
!!
!!
!! DESCRIPTION
!!
!!   Compute the shock adiabat for the Riemann problem
!!
!!  
!! ARGUMENTS
!!
!!  tau0   -  proper specific volume in the pre-shock state (1/rho0)
!!  u0     -  pre-shock normal velocity 
!!  p0     -  pre-shock pressure
!!  g0     -  pre-shock lorentz factor
!!  V0     -  pre-shock specific volume  ( tau0/g0 )
!!  h0     -  pre-shock specific enthalpy
!!  p1     -  input post shock pressure
!!  u1     -  post shock velocity
!!  dudp   -  derivative of u with respect to p, for p = p1
!!  zeta   -  a useful function 
!!  istate - specifies whether the state is left (-1) or right (+1)
!!  
!!***************************************************************************
  
subroutine hy_rhd_shock(tau0, u0, p0, g0, V0, h0, p1, u1,  dudp, zeta, istate)

  use Hydro_data, ONLY : hy_gamma
  implicit none
  
  integer, INTENT(in)  :: istate
  real, INTENT(in)     :: tau0, u0,p0, g0, V0, h0, p1
  real, INTENT(out)    :: u1, dudp, zeta
  
  real  ::  a,b,c, da,db,dc, dp
  real  ::  tau1,h1,d_htau1,g1, j2, dx, gamma_r
  
  dp = p1 - p0
  gamma_r = hy_gamma/(hy_gamma - 1.e0)
  
  !! ***********************************************************
  !!        Use Taub Adiabat to find post-shock Enthalpy
  !!           and mass flux; here j2 --> 1/(j)^2
  !! *********************************************************** 
  
  !!    call rhd_taub_adiabat(dp, p1, tau0, h0, j2, d_htau1)
  
  a = 1.e0 - dp/(gamma_r*p1)
  b = 1.e0 - a;
  c = -h0*(h0 + tau0*dp)
  
  h1   = 0.5d0/a*(-b + sqrt(b*b - 4.e0*a*c))
  tau1 = (h1 - 1.e0)/(gamma_r*p1)
  g1   = 2.e0*h1*gamma_r/(2.e0*h1 - 1.e0)
  
  j2 = h0*gamma_r*tau0 + (h1*tau1 + h0*tau0)*(1.e0/(h0 + h1) - 1.e0)
  j2 = j2/(gamma_r*p1)
  
  d_htau1 = (h1*tau1 + h0*tau0 - g1*h1*tau1)
  d_htau1 = d_htau1/(g1*p1 - dp)
  
  g1   = istate*sqrt(V0*V0 + (1.e0 - u0*u0)*j2)
  zeta = (V0*u0 + g1) / (1.e0 - u0*u0)
  
  !! ***********************************************************
  !!              Get post-shock velocity
  !! *********************************************************** 
  
  b   = 1.e0/(h0*g0 + (p1 - p0)*(zeta*u0 + V0))
  u1  = (h0*g0*u0 + zeta*(p1 - p0))*b
  
  !! ***********************************************************
  !!              Get du/dp for next iteration
  !! *********************************************************** 
  
  a    = -0.5d0*(d_htau1 + j2) / g1
  dudp = (zeta + a - u1*(zeta*u0 + V0 + a*u0))*b
  
end subroutine hy_rhd_shock
 




