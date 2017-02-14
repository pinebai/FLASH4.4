!!****if* source/physics/Hydro/HydroMain/split/MHD_8Wave/hy_8wv_addResistiveFluxes
!!
!! NAME
!!
!!  hy_8wv_addResistiveFluxes
!!
!!
!! SYNOPSIS
!!
!!
!!  hy_8wv_addResistiveFluxes(integer(IN) :: i1,
!!                              integer(IN) :: i2,
!!                              real(IN)    :: B(3,:,:,:),
!!                              real(OUT)   :: Flux(:,:),
!!                              real(IN)    :: visc(:),
!!                              real(IN)    :: x(:),
!!                              real(IN)    :: y(:),
!!                              real(IN)    :: z(:),
!!                              integer(IN) :: nx,
!!                              integer(IN) :: ny,
!!                              integer(IN) :: nz,
!!                              integer(IN) :: sweepDir)
!!
!!
!! DESCRIPTION
!!
!!  Adds resistive fluxe contributions to total MHD fluxes
!!
!!
!! ARGUMENTS
!!
!!  i1,i2      - indices of the line along which the sweep is made
!!  B          - array containing magnetic field components
!!  Flux       - array containing MHD fluxes
!!  visc       - array containing magnetic viscosity coefficients
!!  x,y,z      - coordinate arrays
!!  nx,ny,nz   - sizes of coordinate arrays
!!  sweepDir   - direction of sweep
!!
!!***

subroutine hy_8wv_addResistiveFluxes(i1,i2,B,Flux,visc,x,y,z,nx,ny,nz,sweepDir)

  implicit none

#include "constants.h"
#include "Flash.h"

  !! Argument list ---------------------------------------------
  integer, INTENT(IN) :: i1,i2,nx,ny,nz,sweepDir
  real, DIMENSION(3,nx,ny,nz), INTENT(IN) :: B
  real, DIMENSION(NFLUXES,max(nx,ny,nz)), INTENT(OUT) :: Flux
  real, DIMENSION(max(nx,ny,nz)), INTENT(IN) :: x,y,z,visc
  !! -----------------------------------------------------------

  integer :: i,j,k
  real    :: mvisc,jx,jy,jz,db

  select case(sweepDir)
  case(SWEEP_X)

    j  = i1
    k  = i2

    do i = 3, nx-1

      mvisc = 0.5*(visc(i)+visc(i-1))

      jy = -(B(KAXIS,i,j,k)-B(KAXIS,i-1,j,k))/(x(i)-x(i-1))
      jz =  (B(JAXIS,i,j,k)-B(JAXIS,i-1,j,k))/(x(i)-x(i-1))

      db =  (B(IAXIS,i,j,k)-B(IAXIS,i-1,j,k))/(x(i)-x(i-1))

#if NDIM >= 2
      jz = jz-0.5*(B(IAXIS, i ,j+1, k )-B(IAXIS, i ,j-1, k )+ &
                   B(IAXIS,i-1,j+1, k )-B(IAXIS,i-1,j-1, k ))/(y(j+1)-y(j-1))

      db = db+0.5*(B(JAXIS, i ,j+1, k )-B(JAXIS, i ,j-1, k )+ &
                   B(JAXIS,i-1,j+1, k )-B(JAXIS,i-1,j-1, k ))/(y(j+1)-y(j-1))
#if NDIM == 3
      jy = jy+0.5*(B(IAXIS, i , j ,k+1)-B(IAXIS, i , j ,k-1)+ &
                   B(IAXIS,i-1, j ,k+1)-B(IAXIS,i-1, j ,k-1))/(z(k+1)-z(k-1))
      db = db+0.5*(B(KAXIS, i , j ,k+1)-B(KAXIS, i , j ,k-1)+ &
                   B(KAXIS,i-1, j ,k+1)-B(KAXIS,i-1, j ,k-1))/(z(k+1)-z(k-1))
#endif
#endif

      Flux(MAGX_FLUX,i) = Flux(MAGX_FLUX,i)-mvisc*db
      Flux(MAGY_FLUX,i) = Flux(MAGY_FLUX,i)-mvisc*jz
      Flux(MAGZ_FLUX,i) = Flux(MAGZ_FLUX,i)+mvisc*jy

      Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)-0.5*mvisc* &
       ((B(JAXIS,i,j,k)**2+B(KAXIS,i,j,k)**2-B(JAXIS,i-1,j,k)**2-B(KAXIS,i-1,j,k)**2)/(x(i)-x(i-1)))

#if NDIM >= 2
      Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)+0.5*mvisc* &
        (B(JAXIS, i ,j,k)*(B(IAXIS, i ,j+1, k )-B(IAXIS, i ,j-1, k ))+ &
         B(JAXIS,i-1,j,k)*(B(IAXIS,i-1,j+1, k )-B(IAXIS,i-1,j-1, k )))/(y(j+1)-y(j-1))
#if NDIM == 3
      Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)+0.5*mvisc* &
        (B(KAXIS, i ,j,k)*(B(IAXIS, i , j ,k+1)-B(IAXIS, i , j ,k-1))+ &
         B(KAXIS,i-1,j,k)*(B(IAXIS,i-1, j ,k+1)-B(IAXIS,i-1, j ,k-1)))/(z(k+1)-z(k-1))
#endif
#endif

    end do

#if NDIM >= 2

  case(SWEEP_Y)

    i = i1
    k = i2

    do j = 3, ny-1

      mvisc = 0.5*(visc(j)+visc(j-1))

      jx =  (B(KAXIS,i,j,k)-B(KAXIS,i,j-1,k))/(y(j)-y(j-1))
      jz = -(B(IAXIS,i,j,k)-B(IAXIS,i,j-1,k))/(y(j)-y(j-1))+ &
              0.5*(B(JAXIS,i+1, j , k )-B(JAXIS,i-1, j , k )+ &
                   B(JAXIS,i+1,j-1, k )-B(JAXIS,i-1,j-1, k ))/(x(i+1)-x(i-1))

      db =  (B(JAXIS,i,j,k)-B(JAXIS,i,j-1,k))/(y(j)-y(j-1))+ &
              0.5*(B(IAXIS,i+1, j , k )-B(IAXIS,i-1, j , k )+ &
                   B(IAXIS,i+1,j-1, k )-B(IAXIS,i-1,j-1, k ))/(x(i+1)-x(i-1))

#if NDIM == 3
      jx = jx-0.5*(B(JAXIS, i , j ,k+1)-B(JAXIS, i , j ,k-1)+ &
                   B(JAXIS, i ,j-1,k+1)-B(JAXIS, i ,j-1,k-1))/(z(k+1)-z(k-1))
      db = db+0.5*(B(KAXIS, i , j ,k+1)-B(KAXIS, i , j ,k-1)+ &
                   B(KAXIS, i ,j-1,k+1)-B(KAXIS, i ,j-1,k-1))/(z(k+1)-z(k-1))
#endif

      Flux(MAGY_FLUX,j) = Flux(MAGY_FLUX,j)-mvisc*db
      Flux(MAGX_FLUX,j) = Flux(MAGX_FLUX,j)+mvisc*jz
      Flux(MAGZ_FLUX,j) = Flux(MAGZ_FLUX,j)-mvisc*jx

      Flux(ENER_FLUX,j) = Flux(ENER_FLUX,j)-0.5*mvisc* &
       ((B(IAXIS,i,j,k)**2+B(KAXIS,i,j,k)**2-B(IAXIS,i,j-1,k)**2-B(KAXIS,i,j-1,k)**2)/(y(j)-y(j-1))- &
        (B(IAXIS,i, j ,k)*(B(JAXIS,i+1, j , k )-B(JAXIS,i-1, j , k ))+ &
         B(IAXIS,i,j-1,k)*(B(JAXIS,i+1,j-1, k )-B(JAXIS,i-1,j-1, k )))/(x(i+1)-x(i-1)))

#if NDIM ==3
      Flux(ENER_FLUX,j) = Flux(ENER_FLUX,j)+0.5*mvisc* &
        (B(KAXIS,i, j ,k)*(B(JAXIS, i , j ,k+1)-B(JAXIS, i , j ,k-1))+ &
         B(KAXIS,i,j-1,k)*(B(JAXIS, i ,j-1,k+1)-B(JAXIS, i ,j-1,k-1)))/(z(k+1)-z(k-1))
#endif

    end do

#if NDIM == 3

  case(SWEEP_Z)

    i = i1
    j = i2

    do k = 3, nz-1

      mvisc = 0.5*(visc(k)+visc(k-1))

      jx = -(B(JAXIS,i,j,k)-B(JAXIS,i,j,k-1))/(z(k)-z(k-1))+ &
              0.5*(B(KAXIS, i ,j+1, k )-B(KAXIS, i ,j-1, k )+ &
                   B(KAXIS, i ,j+1,k-1)-B(KAXIS, i ,j-1,k-1))/(y(j+1)-y(j-1))
      jy =  (B(IAXIS,i,j,k)-B(IAXIS,i,j,k-1))/(z(k)-z(k-1))- &
              0.5*(B(KAXIS,i+1, j , k )-B(KAXIS,i-1, j , k )+ &
                   B(KAXIS,i+1, j ,k-1)-B(KAXIS,i-1, j ,k-1))/(x(i+1)-x(i-1))

      db =  (B(KAXIS,i,j,k)-B(KAXIS,i,j,k-1))/(z(k)-z(k-1))+ &
              0.5*(B(JAXIS, i ,j+1, k )-B(JAXIS, i ,j-1, k )+ &
                   B(JAXIS, i ,j+1,k-1)-B(JAXIS, i ,j-1,k-1))/(y(j+1)-y(j-1))+ &
              0.5*(B(IAXIS,i+1, j , k )-B(IAXIS,i-1, j , k )+ &
                   B(IAXIS,i+1, j ,k-1)-B(IAXIS,i-1, j ,k-1))/(x(i+1)-x(i-1))

      Flux(MAGZ_FLUX,k) = Flux(MAGZ_FLUX,k)-mvisc*db
      Flux(MAGX_FLUX,k) = Flux(MAGX_FLUX,k)-mvisc*jy
      Flux(MAGY_FLUX,k) = Flux(MAGY_FLUX,k)+mvisc*jx

      Flux(MAGX_FLUX:MAGZ_FLUX,k) = Flux(MAGX_FLUX:MAGZ_FLUX,k)-mvisc*(B(:,i,j,k)-B(:,i,j,k-1))/(z(k)-z(k-1))

      Flux(ENER_FLUX,k) = Flux(ENER_FLUX,k)-0.5*mvisc* &
       ((B(IAXIS,i,j,k)**2+B(JAXIS,i,j,k)**2-B(IAXIS,i,j,k-1)**2-B(JAXIS,i,j,k-1)**2)/(z(k)-z(k-1))- &
        (B(IAXIS,i,j, k )*(B(KAXIS,i+1, j , k )-B(KAXIS,i-1, j , k ))+ &
         B(IAXIS,i,j,k-1)*(B(KAXIS,i+1, j ,k-1)-B(KAXIS,i-1, j ,k-1)))/(x(i+1)-x(i-1))- &
        (B(JAXIS,i,j, k )*(B(KAXIS, i ,j+1, k )-B(KAXIS, i ,j-1, k ))+ &
         B(JAXIS,i,j,k-1)*(B(KAXIS, i ,j+1,k-1)-B(KAXIS, i ,j-1,k-1)))/(y(j+1)-y(j-1)))

    end do

#endif
#endif

  end select

end subroutine hy_8wv_addResistiveFluxes
