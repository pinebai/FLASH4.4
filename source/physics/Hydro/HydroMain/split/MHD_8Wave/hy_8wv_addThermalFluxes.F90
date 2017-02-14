!!****if* source/physics/Hydro/HydroMain/split/MHD_8Wave/hy_8wv_addThermalFluxes
!!
!! NAME
!!
!!  hy_8wv_addThermalFluxes
!!
!!
!! SYNOPSIS
!!
!!
!!  hy_8wv_addThermalFluxes(integer(IN) :: i1,
!!                         integer(IN) :: i2,
!!                         real(IN)    :: T(3,:,:,:),
!!                         real(OUT)   :: Flux(:,:),
!!                         real(IN)    :: cond(:),
!!                         real(IN)    :: x(:),
!!                         real(IN)    :: y(:),
!!                         real(IN)    :: z(:),
!!                         integer(IN) :: nx,
!!                         integer(IN) :: ny,
!!                         integer(IN) :: nz,
!!                         integer(IN) :: sweepDir)
!!
!! DESCRIPTION
!!
!!  Adds thermal flux contribution to total MHD fluxes
!!
!! ARGUMENTS
!!
!!  i1,i2      - indices of the line along which the sweep is made
!!  T          - array containing temperature
!!  Flux       - array containing MHD fluxes
!!  cond       - array containing thermal conductivity coefficient
!!  x,y,z      - coordinate arrays
!!  nx,ny,nz   - sizes of coordinate arrays
!!  sweepDir   - direction of sweep
!!
!!***

subroutine hy_8wv_addThermalFluxes(i1,i2,T,Flux,cond,x,y,z,nx,ny,nz,sweepDir)

  implicit none

#include "constants.h"
#include "Flash.h"

  !! Argument list ---------------------------------------------
  integer, INTENT(IN) :: i1,i2,nx,ny,nz,sweepDir
  real, DIMENSION(nx,ny,nz), INTENT(IN) :: T
  real, DIMENSION(NFLUXES,max(nx,ny,nz)), INTENT(OUT) :: Flux
  real, DIMENSION(max(nx,ny,nz)), INTENT(IN) :: x,y,z,cond
  !! -----------------------------------------------------------

  integer :: i,j,k
  real    :: kappa_th,idx,idy,idz

  select case(sweepDir)
  case(SWEEP_X)

    j = i1
    k = i2

    do i = 3, nx-1

      idx = 1./(x(i)-x(i-1))

      kappa_th = 0.5*(cond(i)+cond(i-1))

      Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)-idx*kappa_th*(T(i,j,k)-T(i-1,j,k))

    end do

#if NDIM >= 2

  case(SWEEP_Y)

    i = i1
    k = i2

    do j = 3, ny-1

      idy = 1./(y(j)-y(j-1))

      kappa_th = 0.5*(cond(j)+cond(j-1))

      Flux(ENER_FLUX,j) = Flux(ENER_FLUX,j)-idy*kappa_th*(T(i,j,k)-T(i,j-1,k))

    end do

#if NDIM == 3

  case(SWEEP_Z)

    i = i1
    j = i2

    do k = 3, nz-1

      idz = 1./(z(k)-z(k-1))

      kappa_th = 0.5*(cond(k)+cond(k-1))

      Flux(ENER_FLUX,k) = Flux(ENER_FLUX,k)-idz*kappa_th*(T(i,j,k)-T(i,j,k-1))

    end do

#endif
#endif

  end select

end subroutine hy_8wv_addThermalFluxes
