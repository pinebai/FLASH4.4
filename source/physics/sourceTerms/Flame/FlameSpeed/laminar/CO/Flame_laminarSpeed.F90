!!****if* source/physics/sourceTerms/Flame/FlameSpeed/laminar/CO/Flame_laminarSpeed
!!
!! NAME
!!
!!  Flame_laminarSpeed
!!
!! SYNOPSIS
!!
!!  call Flame_laminarSpeed(real(in) :: dens,
!!                          real(out) :: s,
!!                          real, optional(out) :: ds,
!!                          real, dimension(:), optional(in) :: info)
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!!
!!
!!  In this implementation the info array is ignored, since the fitting form
!!  is only valid for X_C12 = X_O16 = 0.5
!!
!! This is Shimon Asida's fit of Laminar flame speed as function of fuel density
!! From data tabulated in Timmes & Woosley 1992, ApJ, 396, 649
!! Fit only for carbon fraction Xc = 0.5 
!!
!! ARGUMENTS
!!
!!   dens : 
!!
!!   s : 
!!
!!   ds : 
!!
!!   info : 
!!
!!
!!
!!***


subroutine Flame_laminarSpeed(dens, s, ds, info)

  implicit none

  real, intent(in)   :: dens
  real, intent(out)  :: s
  real, optional, intent(out) :: ds !DEV: no meaningful output provided
  real, dimension(:), optional, intent(in) :: info

  real, parameter    :: c1 = -43.0e0
  real, parameter    :: c2 =  4.534e0
  real, parameter    :: c3 = -0.08333e0

  real :: ldens

  ldens = log(dens)
  s = exp(c1 + c2*ldens + c3*(ldens**2)) 

  if (present(ds)) ds = 0.0     !DEV: no meaningful output provided

end subroutine Flame_laminarSpeed
