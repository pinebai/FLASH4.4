!!****if* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_flux
!! 
!! NAME
!! 
!!  hy_rhd_flux 
!! 
!! SYNOPSIS
!!  hy_rhd_flux( real   (IN) :: u(:,;),
!!               real  (OUT) :: flux(:,:),
!!               real  (OUT) :: ucns(:,:),
!!               integer(IN) :: ibeg,
!!               integer(IN) :: iend,
!!               integer(IN) :: n,
!!               integer(IN) :: dir)
!!
!! DESCRIPTION
!!
!!   Compute the dir-component of the relativistic fluxes fluxes from
!!   a 1-D vector of primitive quantities u. For example, when dir = 1
!!   the flux vector is given as:
!!
!!      Fx = (D v_x  ;  m_x v_x + p  ;  m_y v_x  ;  m_z v_x  ;  m_x)
!!
!!    where  D is the Lab-density, m is the momentum, v is the velocity and
!!    p is the pressure.
!!
!!
!! ARGUMENTS
!! 
!!  u     -     a 1-D vector of primitive quantities
!!  flux  -     a 1-D vector containing the fluxes
!!  ucns  -     a 1-D vector containing the conservative map of u
!!  ibeg  -     initial point
!!  iend  -     final point
!!  n     -     number of points
!!  dir   -     flux component; dir = 1,2,3 for the 3 axis
!! 
!!***

subroutine hy_rhd_flux (u, flux, ucns, ibeg, iend, n, dir)  

  use Hydro_data, ONLY : hy_gamma
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "RHD.h"  

  integer, INTENT(in)                  :: ibeg, iend, n, dir
  real, DIMENSION(NUNK_VARS,n), INTENT(in)  :: u
  real, DIMENSION(NUNK_VARS,n), INTENT(out) :: ucns
  real, DIMENSION(NFLUXES,n), INTENT(out) :: flux
  integer             :: i, ivn
  real                :: lor2, vn, vel2, rhoh
  integer  :: VELN_FLUX,VELN_VAR


  select case(dir)
  case (SWEEP_X)
     VELN_FLUX = XMOM_FLUX
     VELN_VAR  = VELX_VAR
  case (SWEEP_Y)
     VELN_FLUX = YMOM_FLUX
     VELN_VAR  = VELY_VAR
  case (SWEEP_Z)
     VELN_FLUX = ZMOM_FLUX
     VELN_VAR  = VELZ_VAR
  end select


  flux = 0.d0

  do i = ibeg, iend
    vn   = u(VELN_VAR,i)
    vel2 = u(VELX_VAR,i)*u(VELX_VAR,i) + u(VELY_VAR,i)*u(VELY_VAR,i) + u(VELZ_VAR,i)*u(VELZ_VAR,i)

    if (vel2 .ge. 1.d0) then 
      call Driver_abortFlash("v^2 > 1 in hy_rhd_flux")
    end if

    lor2  = 1.d0/(1.d0 - vel2)
    rhoh = u(DENS_VAR,i) + hy_gamma/(hy_gamma - 1.d0)*u(PRES_VAR,i)
       
    ucns(DENS_VAR,i)     = u(DENS_VAR,i)*sqrt(lor2)
    ucns(VELX_VAR:VELZ_VAR,i) = rhoh*lor2*u(VELX_VAR:VELZ_VAR,i)
    ucns(ENER_VAR,i)     = rhoh*lor2 - u(PRES_VAR,i)

    flux(DENS_FLUX,i)     = ucns(DENS_VAR,i)*vn
    flux(XMOM_FLUX:ZMOM_FLUX,i) = ucns(VELX_VAR:VELZ_VAR,i)*vn
    flux(VELN_FLUX,i)     = flux(VELN_FLUX,i) + u(PRES_VAR,i)
    flux(ENER_FLUX,i)     = ucns(VELN_VAR,i)
  end do

end subroutine hy_rhd_flux






