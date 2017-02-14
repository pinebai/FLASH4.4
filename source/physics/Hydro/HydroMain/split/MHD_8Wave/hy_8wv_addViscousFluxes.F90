!!****if* source/physics/Hydro/HydroMain/split/MHD_8Wave/hy_8wv_addViscousFluxes
!!
!! NAME
!!
!!  hy_8wv_addViscousFluxes
!!
!!
!! SYNOPSIS
!!
!!
!!  hy_8wv_addViscousFluxes(integer(IN) :: i1,
!!                            integer(IN) :: i2,
!!                            real(IN)    :: B(3,:,:,:),
!!                            real(OUT)   :: Flux(:,:),
!!                            real(IN)    :: visc(:),
!!                            real(IN)    :: x(:),
!!                            real(IN)    :: y(:),
!!                            real(IN)    :: z(:),
!!                            integer(IN) :: nx,
!!                            integer(IN) :: ny,
!!                            integer(IN) :: nz,
!!                            integer(IN) :: sweepDir)
!!
!! DESCRIPTION
!!
!!  Adds viscous flux contributions to total MHD fluxes
!!
!! ARGUMENTS
!!
!!  i1,i2      - indices of the line along which the sweep is made
!!  B          - array containing magnetic field components
!!  Flux       - array containing MHD fluxes
!!  visc       - array containing viscosity coefficients
!!  x,y,z      - coordinate arrays
!!  nx,ny,nz   - sizes of coordinate arrays
!!  sweepDir   - direction of sweep
!!
!!***

subroutine hy_8wv_addViscousFluxes(i1,i2,V,Flux,visc,x,y,z,nx,ny,nz,sweepDir)

  implicit none

#include "constants.h"
#include "Flash.h"

  !! Argument list ---------------------------------------------
  integer, INTENT(IN) :: i1,i2,nx,ny,nz,sweepDir
  real, DIMENSION(3,nx,ny,nz), INTENT(IN) :: V
  real, DIMENSION(NFLUXES,max(nx,ny,nz)), INTENT(OUT) :: Flux
  real, DIMENSION(max(nx,ny,nz)), INTENT(IN) :: x,y,z,visc
  !! -----------------------------------------------------------

  integer :: i,j,k
  real    :: vsc,idx,idy,idz

  select case(sweepDir)
  case(SWEEP_X)

    j = i1
    k = i2

#if NDIM >= 2
    idy = 1./(y(j+1)-y(j-1))
#if NDIM == 3
    idz = 1./(z(k+1)-z(k-1))
#endif
#endif

    do i = 3, nx-1

      idx = 1./(x(i)-x(i-1))

      vsc = 0.5*(visc(i)+visc(i-1))

      Flux(XMOM_FLUX:ZMOM_FLUX,i) = Flux(XMOM_FLUX:ZMOM_FLUX,i)-idx*vsc*(V(:,i,j,k)-V(:,i-1,j,k))

      Flux(XMOM_FLUX,i) = Flux(XMOM_FLUX,i)-idx*vsc*(V(IAXIS,i,j,k)-V(IAXIS,i-1,j,k))/3.

      Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)-idx*vsc* &
        ((4./3.)*V(IAXIS, i ,j,k)**2+V(JAXIS, i ,j,k)**2+V(KAXIS, i ,j,k)**2- &
         (4./3.)*V(IAXIS,i-1,j,k)**2-V(JAXIS,i-1,j,k)**2-V(KAXIS,i-1,j,k)**2)/2.

#if NDIM >= 2

      Flux(XMOM_FLUX,i) = Flux(XMOM_FLUX,i)+idy*vsc*(V(JAXIS, i ,j+1,k)-V(JAXIS, i ,j-1,k)+ &
                                         V(JAXIS,i-1,j+1,k)-V(JAXIS,i-1,j-1,k))/3.

      Flux(YMOM_FLUX,i) = Flux(YMOM_FLUX,i)-idy*vsc*(V(IAXIS, i ,j+1,k)-V(IAXIS, i ,j-1,k)+ &
                                         V(IAXIS,i-1,j+1,k)-V(IAXIS,i-1,j-1,k))/2.

      Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)+idy*vsc* &
        ((2./3.)*(V(IAXIS, i ,j,k)*(V(JAXIS, i ,j+1,k)-V(JAXIS, i ,j-1,k))+ &
                  V(IAXIS,i-1,j,k)*(V(JAXIS,i-1,j+1,k)-V(JAXIS,i-1,j-1,k))) &
                -(V(JAXIS, i ,j,k)*(V(IAXIS, i ,j+1,k)-V(IAXIS, i ,j-1,k))+ &
                  V(JAXIS,i-1,j,k)*(V(IAXIS,i-1,j+1,k)-V(IAXIS,i-1,j-1,k))))/2.

#if NDIM == 3

      Flux(XMOM_FLUX,i) = Flux(XMOM_FLUX,i)+idz*vsc*(V(KAXIS, i ,j,k+1)-V(KAXIS, i ,j,k-1)+ &
                                         V(KAXIS,i-1,j,k+1)-V(KAXIS,i-1,j,k-1))/3.

      Flux(ZMOM_FLUX,i) = Flux(ZMOM_FLUX,i)-idz*vsc*(V(IAXIS, i ,j,k+1)-V(IAXIS, i ,j,k-1)+ &
                                         V(IAXIS,i-1,j,k+1)-V(IAXIS,i-1,j,k-1))/2.

      Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)+idz*vsc* &
        ((2./3.)*(V(IAXIS, i ,j,k)*(V(KAXIS, i ,j,k+1)-V(KAXIS, i ,j,k-1))+ &
                  V(IAXIS,i-1,j,k)*(V(KAXIS,i-1,j,k+1)-V(KAXIS,i-1,j,k-1))) &
                -(V(KAXIS, i ,j,k)*(V(IAXIS, i ,j,k+1)-V(IAXIS, i ,j,k-1))- &
                  V(KAXIS,i-1,j,k)*(V(IAXIS,i-1,j,k+1)-V(IAXIS,i-1,j,k-1))))/2.

#endif
#endif

    end do

#if NDIM >= 2

  case(SWEEP_Y)

    i = i1
    k = i2

    idx = 1./(x(i+1)-x(i-1))
#if NDIM == 3
    idz = 1./(z(k+1)-z(k-1))
#endif

    do j = 3, ny-1

      idy = 1./(y(j)-y(j-1))
      vsc = 0.5*(visc(j)+visc(j-1))

      Flux(XMOM_FLUX:ZMOM_FLUX,j) = Flux(XMOM_FLUX:ZMOM_FLUX,j)-idy*vsc*(V(:,i,j,k)-V(:,i,j-1,k))

      Flux(XMOM_FLUX,j) = Flux(XMOM_FLUX,j)-idx*vsc*(V(JAXIS,i+1, j ,k)-V(JAXIS,i-1, j ,k)+ &
                                         V(JAXIS,i+1,j-1,k)-V(JAXIS,i-1,j-1,k))/2.

      Flux(YMOM_FLUX,j) = Flux(YMOM_FLUX,j)-idy*vsc*(V(JAXIS,i,j,k)-V(JAXIS,i,j-1,k))/3.+ &
                                idx*vsc*(V(IAXIS,i+1, j ,k)-V(IAXIS,i-1, j ,k)+ &
                                         V(IAXIS,i+1,j-1,k)-V(IAXIS,i-1,j-1,k))/3.

      Flux(ENER_FLUX,j) = Flux(ENER_FLUX,j)+idx*vsc* &
        ((2./3.)*(V(JAXIS,i, j ,k)*(V(IAXIS,i+1, j ,k)-V(IAXIS,i-1, j ,k))+ &
                  V(JAXIS,i,j-1,k)*(V(IAXIS,i+1,j-1,k)-V(IAXIS,i-1,j-1,k))) &
                -(V(IAXIS,i, j ,k)*(V(JAXIS,i+1, j ,k)-V(JAXIS,i-1, j ,k))+ &
                  V(IAXIS,i,j-1,k)*(V(JAXIS,i+1,j-1,k)-V(JAXIS,i-1,j-1,k))))/2. &
                               -idy*vsc* &
        ((4./3.)*V(JAXIS,i, j ,k)**2+V(IAXIS,i, j ,k)**2+V(KAXIS,i, j ,k)**2- &
         (4./3.)*V(JAXIS,i,j-1,k)**2-V(IAXIS,i,j-1,k)**2-V(KAXIS,i,j-1,k)**2)/2.

#if NDIM == 3

      Flux(YMOM_FLUX,j) = Flux(YMOM_FLUX,j)+idz*vsc*(V(KAXIS,i, j ,k+1)-V(KAXIS,i, j ,k-1)+ &
                                         V(KAXIS,i,j-1,k+1)-V(KAXIS,i,j-1,k-1))/3.

      Flux(ZMOM_FLUX,j) = Flux(ZMOM_FLUX,j)-idz*vsc*(V(JAXIS,i, j ,k+1)-V(JAXIS,i, j ,k-1)+ &
                                         V(JAXIS,i,j-1,k+1)-V(JAXIS,i,j-1,k-1))/2.

      Flux(ENER_FLUX,j) = Flux(ENER_FLUX,j)+idz*vsc* &
        ((2./3.)*(V(JAXIS,i, j ,k)*(V(KAXIS,i, j ,k+1)-V(KAXIS,i, j ,k-1))+ &
                  V(JAXIS,i,j-1,k)*(V(KAXIS,i,j-1,k+1)-V(KAXIS,i,j-1,k-1))) &
                -(V(KAXIS,i, j ,k)*(V(JAXIS,i, j ,k+1)-V(JAXIS,i, j ,k-1))- &
                  V(KAXIS,i,j-1,k)*(V(JAXIS,i,j-1,k+1)-V(JAXIS,i,j-1,k-1))))/2.

#endif

    end do

#if NDIM == 3

  case(SWEEP_Z)

    i = i1
    j = i2

    idx = 1./(x(i+1)-x(i-1))
    idy = 1./(y(j+1)-y(j-1))

    do k = 3, nz-1

      idz = 1./(z(k)-z(k-1))

      vsc = 0.5*(visc(k)+visc(k-1))

      Flux(XMOM_FLUX:ZMOM_FLUX,k) = Flux(XMOM_FLUX:ZMOM_FLUX,k)-idz*vsc*(V(:,i,j,k)-V(:,i,j,k-1))

      Flux(XMOM_FLUX,k) = Flux(XMOM_FLUX,k)-idx*vsc*(V(KAXIS,i+1,j, k )-V(KAXIS,i-1,j, k )+ &
                                         V(KAXIS,i+1,j,k-1)-V(KAXIS,i-1,j,k-1))/2.

      Flux(YMOM_FLUX,k) = Flux(YMOM_FLUX,k)-idy*vsc*(V(KAXIS,i,j+1, k )-V(KAXIS,i,j-1, k )+ &
                                         V(KAXIS,i,j+1,k-1)-V(KAXIS,i,j-1,k-1))/2.

      Flux(ZMOM_FLUX,k) = Flux(ZMOM_FLUX,k)-idz*vsc*(V(KAXIS,i,j,k)-V(KAXIS,i,j,k-1))/3.+ &
                                idx*vsc*(V(IAXIS,i+1,j, k )-V(IAXIS,i-1,j, k )+ &
                                         V(IAXIS,i+1,j,k-1)-V(IAXIS,i-1,j,k-1))/3.+ &
                                idy*vsc*(V(JAXIS,i,j+1, k )-V(JAXIS,i,j-1, k )+ &
                                         V(JAXIS,i,j+1,k-1)-V(JAXIS,i,j-1,k-1))/3.

      Flux(ENER_FLUX,k) = Flux(ENER_FLUX,k)+idx*vsc* &
        ((2./3.)*(V(KAXIS,i,j, k )*(V(IAXIS,i+1,j, k )-V(IAXIS,i-1,j, k ))+ &
                  V(KAXIS,i,j,k-1)*(V(IAXIS,i+1,j,k-1)-V(IAXIS,i-1,j,k-1))) &
                -(V(IAXIS,i,j, k )*(V(KAXIS,i+1,j, k )-V(KAXIS,i-1,j, k ))+ &
                  V(IAXIS,i,j,k-1)*(V(KAXIS,i+1,j,k-1)-V(KAXIS,i-1,j,k-1))))/2. &
                               +idy*vsc* &
        ((2./3.)*(V(KAXIS,i,j, k )*(V(JAXIS,i,j+1, k )-V(JAXIS,i,j-1, k ))+ &
                  V(KAXIS,i,j,k-1)*(V(JAXIS,i,j+1,k-1)-V(JAXIS,i,j-1,k-1))) &
                -(V(JAXIS,i,j, k )*(V(KAXIS,i,j+1, k )-V(KAXIS,i,j-1, k ))- &
                  V(JAXIS,i,j,k-1)*(V(KAXIS,i,j+1,k-1)-V(KAXIS,i,j-1,k-1))))/2. &
                               -idz*vsc* &
        ((4./3.)*V(KAXIS,i,j, k )**2+V(IAXIS,i,j, k )**2+V(JAXIS,i,j, k )**2- &
         (4./3.)*V(KAXIS,i,j,k-1)**2-V(IAXIS,i,j,k-1)**2-V(JAXIS,i,j,k-1)**2)/2.

    end do

#endif
#endif

  end select

end subroutine hy_8wv_addViscousFluxes
