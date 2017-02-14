!!****if* source/physics/Cosmology/CosmologyMain/csm_integrateFriedmann
!!
!! NAME
!!
!!  csm_integrateFriedmann
!!
!! SYNOPSIS
!!
!!  csm_integrateFriedmann(real, intent(IN)  :: a,
!!                         real, intent(IN)  :: t,
!!                         real, intent(IN)  :: dt,
!!                         real, intent(OUT)  :: anew)
!!
!! DESCRIPTION
!!
!!  This uses a fourth-order Runge-Kutta method to integrate the Friedmann
!!  equation and compute the new cosmological redshift scaling factor.  This 
!!  was handled by a third-order method provided by SVODE in Flash 2. This will
!!  advance the Friedmann equation one step to calculate the scale factor.
!!
!! ARGUMENTS
!!
!!   a : The initial scaling factor
!!
!!   t : The current time to solve for
!!
!!   dt : The current change in timestep
!!
!!   anew : The new scaling factor
!!
!!
!!
!!***


subroutine csm_integrateFriedmann (a,t,dt,anew)


  implicit none
  
  real, intent(IN)  :: a
  real, intent(IN)  :: t
  real, intent(IN)  :: dt
  real, intent(OUT) :: anew

  real :: dt2, dt6, dtn
  real :: a1
  real :: dadt,dadt1,dadt2

  ! Use a Runge-Kutta method

  dt2 = 0.5*dt
  dt6 = dt/6.
  dtn = dt + dt2
  
  call csm_friedmannDeriv(a,dadt)

  a1 = a + dt2*dadt

  call csm_friedmannDeriv(a1,dadt1)
  
  a1 = a + dt2*dadt1

  call csm_friedmannDeriv(a1,dadt2)

  a1 = a + dt*dadt2
  dadt2 = dadt1 + dadt2

  call csm_friedmannDeriv(a1,dadt1)

  anew = a + dt6*(dadt + dadt1 + 2.*dadt2)
  
  
return
end subroutine csm_integrateFriedmann
