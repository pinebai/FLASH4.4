!!****if* source/Simulation/SimulationMain/SedovSelfGravity/sim_setAnalyticSedov
!!
!! NAME
!!
!!  sim_setAnalyticSedov
!!
!!
!! SYNOPSIS
!!   call sim_setAnalyticSedov(integer(IN) :: N,
!!                        real(N)(inout) :: r, 
!!                        real(N)(out) :: rho, 
!!                        real(N)(out) :: p, 
!!                        real(N)(out) :: v, 
!!                        real(in) :: t, 
!!                        real(in) :: gamma, 
!!                        real(in) :: E,
!!                        real(in) :: p_ambient, 
!!                        real(in) :: rho_ambient)
!!
!! DESCRIPTION
!!     Given a set of arrays to store radius, density, pressure, and
!!     velocity, together with a time and a ratio of specific heats,
!!     generate the analytical solution to the Sedov problem.
!!
!!     Currently this routine is limited to gamma=1.4 (there is hardwired
!!     dimensionless constant beta) and will complain if it gets
!!     anything else.
!!
!! ARGUMENTS
!!   N - number of data points in the profile
!!   r - the radius for the profile
!!   rho - density profile
!!   p - pressure profile
!!   v - velocity profile
!!   t - initial time since explosion
!!   gamma - the ideal gas gamma
!!   E - explosion energy
!!   p_ambient - ambient pressure
!!   rho_ambient - ambient density
!!
!!***

subroutine sim_setAnalyticSedov (N, r, rho, p, v, t, gamma, E, &
                                 p_ambient, rho_ambient)

!======================================================================
  implicit none
  integer, intent(IN) :: N
  real,dimension(N),intent(inout) :: r
  real,dimension(N),intent(out) ::rho,p,v
  real, intent(IN) ::  t, gamma, E, p_ambient, rho_ambient

  real    beta, R0, dr, xi, nu1, nu2, nu3, nu4, nu5, VV, G, Z, &
       kappa, zeta, epsilon, c2sqr, k, gamp1, gam7, gamm1
  integer i

!==========================================================================

  if (gamma /= 1.4) then
     write (*,*) 'Warning!  init_block() found gamma<>1.4 and t>0.'
     write (*,*) '          Analytical initial conditions will be'
     write (*,*) '          wrong.  Assuming beta=1.033...'
  endif
  
!               Compute dimensionless scaling constant and explosion radius.
  
  beta = 1.033
  R0   = beta * (E*t*t/rho_ambient)**0.2
  
  !               Compute exponents for self-similar solution.
  
  nu1 = - (13.*gamma*gamma - 7.*gamma + 12.) / &
       ((3.*gamma - 1.) * (2.*gamma+1.))
  nu2 = 5. * (gamma - 1.) / (2.*gamma + 1.)
  nu3 = 3. / (2.*gamma + 1.)
  nu4 = - nu1 / (2. - gamma)
  nu5 = - 2. / (2. - gamma)
  
  !               Other useful combinations of gamma.
  
  gamp1 = gamma + 1.E0
  gamm1 = gamma - 1.E0
  gam7  = 7.E0 - gamma
  k     = gamp1 / gamm1
  
  !===============================================================

  !               Generate the solution.
  
  do i = 1, N
     
     xi = r(i) / R0                ! Fraction of explosion radius.
     
     if (xi .le. 1.) then          ! Inside the explosion.
        
        ! Compute V(xi) using bisection.
        ! See Landau & Lifshitz for notation.
        call sim_computeSedovV (xi, gamma, nu1, nu2, VV)
        
        G = k * (k*(gamma*VV-1.))**nu3 * & 
                (gamp1/gam7*(5.-(3.*gamma-1.)*VV))**nu4 * & 
                (k*(1.-VV))**nu5
        rho(i) = rho_ambient * G
        v(i)   = 2.*r(i)*VV / (5.*t)
        
        if (xi .le. 1.E-6) then     ! Use asymptotic r->0 solution.
           kappa = ( (0.5*gamp1/gamma)**2. * & 
                     (gamp1/gam7* & 
                     (5.-(3.*gamma-1.)/gamma))**nu1 )**(1./nu2) * & 
                      gamm1/gamp1 / gamma
           epsilon = k**(nu5+1.) * (k*gamma*kappa)**nu3 * & 
                     (gamp1/gam7*(3.*gamma-1.))**nu4 * & 
                     ((2.*gamma+1)/gamma/(3.*gamma-1.)) * & 
                     (gamm1/gamma)
           zeta = gamm1*gamm1/(2.*gamma*gamma*kappa)
           p(i) = rho_ambient/gamma * 0.16*(R0/t)**2 * epsilon*zeta
        else
           Z = gamma * gamm1 * (1.-VV) * VV**2 / (2.*(gamma*VV-1.))
           c2sqr = 0.16*(r(i)/t)**2 * Z
           p(i) = rho(i) * c2sqr / gamma
        endif
        
     else                          ! Outside the explosion.
        
        rho(i) = rho_ambient
        p(i)   = p_ambient
        v(i)   = 0.
        
     endif
     
  enddo
  
!=======================================================================
  
  return
end subroutine sim_setAnalyticSedov
