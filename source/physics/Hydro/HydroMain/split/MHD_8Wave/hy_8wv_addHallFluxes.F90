!!****if* source/physics/Hydro/HydroMain/split/MHD_8Wave/hy_8wv_addHallFluxes
!!
!! NAME
!!
!!  hy_8wv_addHallFluxes
!!
!!
!! SYNOPSIS
!!
!!
!!  hy_8wv_addHallFluxes(integer(IN)   :: i1,
!!                         integer(IN)   :: i2,
!!                         real(IN)      :: dens(:),
!!                         real(IN)      :: B(3,:,:,:),
!!                         real(OUT)     :: Flux(:,:),
!!                         real(OUT)     :: speed,
!!                         real(IN)      :: x(:),
!!                         real(IN)      :: y(:),
!!                         real(IN)      :: z(:),
!!                         integer(IN)   :: nx,
!!                         integer(IN)   :: ny,
!!                         integer(IN)   :: nz,
!!                         integer(IN)   :: sweepDir)
!!
!!
!! DESCRIPTION
!!
!!  Adds Hall flux contributions to total MHD fluxes
!!
!!
!! ARGUMENTS
!!
!!  i1,i2       - indices of the line along which the sweep is made
!!  dens        - array containing density
!!  B           - array containing magnetic field components
!!  Flux        - array containing MHD fluxes
!!  speed       - Current "speed" needed to limit timestep
!!  x,y,z       - coordinate arrays
!!  nx,ny,nz    - sizes of coordinate arrays
!!  sweepDir    - direction of sweep
!!
!!***

subroutine hy_8wv_addHallFluxes(i1,i2,dens,B,Flux,speed,x,y,z,nx,ny,nz,sweepDir)

  use Hydro_data, ONLY : hy_hall_parameter, hy_hyperResistivity

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, INTENT(IN) :: i1,i2,nx,ny,nz,sweepDir
  real, DIMENSION(3,nx,ny,nz), INTENT(IN) :: B
  real, DIMENSION(NFLUXES,max(nx,ny,nz)), INTENT(OUT) :: Flux
  real, DIMENSION(max(nx,ny,nz)), INTENT(IN) :: x,y,z,dens
  real, INTENT(OUT) :: speed

  integer :: i,j,k
  real, DIMENSION(3) :: Ba,Ja
  real :: d,hd3

  speed = TINY(speed)

  select case(sweepDir)

  case(SWEEP_X)
    j = i1
    k = i2

    do i = 3, nx-1
      Ja(IAXIS) = 0.0
      Ja(JAXIS) = -(B(KAXIS,i,j,k)-B(KAXIS,i-1,j,k))/(x(i)-x(i-1))
      Ja(KAXIS) =  (B(JAXIS,i,j,k)-B(JAXIS,i-1,j,k))/(x(i)-x(i-1))
#if N_DIM >= 2
      Ja(IAXIS) = Ja(IAXIS)+0.5*(B(KAXIS, i ,j+1, k ) - B(KAXIS, i ,j-1, k )+ &
                                 B(KAXIS,i-1,j+1, k ) - B(KAXIS,i-1,j-1, k ))/(y(j+1)-y(j-1))
      Ja(KAXIS) = Ja(KAXIS)-0.5*(B(IAXIS, i ,j+1, k ) - B(IAXIS, i ,j-1, k )+ &
                                 B(IAXIS,i-1,j+1, k ) - B(IAXIS,i-1,j-1, k ))/(y(j+1)-y(j-1))
#if N_DIM == 3
      Ja(IAXIS) = Ja(IAXIS)-0.5*(B(JAXIS, i , j ,k+1) - B(JAXIS, i , j ,k-1)+ &
                                 B(JAXIS,i-1, j ,k+1) - B(JAXIS,i-1, j ,k-1))/(z(k+1)-z(k-1))
      Ja(JAXIS) = Ja(JAXIS)+0.5*(B(IAXIS, i , j ,k+1) - B(IAXIS, i , j ,k-1)+ &
                                 B(IAXIS,i-1, j ,k+1) - B(IAXIS,i-1, j ,k-1))/(z(k+1)-z(k-1))
#endif
#endif

      d  = 0.5*(dens(i)+dens(i-1))
      Ba = 0.5*(B(:,i,j,k)+B(:,i-1,j,k))
      Ja  = (hy_hall_parameter/d)*Ja

      Flux(MAGX_FLUX:MAGZ_FLUX,i) = Flux(MAGX_FLUX:MAGZ_FLUX,i)+Ba(IAXIS)*Ja-Ja(IAXIS)*Ba
      Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)+Ba(IAXIS)*dot_product(Ja,Ba)-Ja(IAXIS)*dot_product(Ba,Ba)

      ! Add hyper-resistive terms
      hd3 = hy_hyperResistivity/(x(i)-x(i-1))**3
      Flux(MAGX_FLUX:MAGZ_FLUX,i) = &
           Flux(MAGX_FLUX:MAGZ_FLUX,i)+hd3*(B(:,i+1,j,k)-3.*B(:,i,j,k)+3.*B(:,i-1,j,k)-B(:,i-2,j,k))

      speed  = max(speed,abs(Ja(IAXIS)),abs(Ja(JAXIS)),abs(Ja(KAXIS)),8.*hd3)
    end do

#if N_DIM >= 2

  case(SWEEP_Y)
    i = i1
    k = i2

    do j = 3, ny-1
      Ja(IAXIS) =  (B(KAXIS,i,j,k)-B(KAXIS,i,j-1,k))/(y(j)-y(j-1))
      Ja(JAXIS) = -0.5*(B(KAXIS,i+1, j ,k) - B(KAXIS,i-1, j ,k)+ &
                        B(KAXIS,i+1,j-1,k) - B(KAXIS,i-1,j-1,k))/(x(i+1)-x(i-1))
      Ja(KAXIS) =  0.5*(B(JAXIS,i+1, j ,k) - B(JAXIS,i-1, j ,k)+ &
                        B(JAXIS,i+1,j-1,k) - B(JAXIS,i-1,j-1,k))/(x(i+1)-x(i-1))- &
                       (B(IAXIS, i , j ,k) - B(IAXIS, i ,j-1,k))/(y( j )-y(j-1))
#if N_DIM == 3
      Ja(IAXIS) = Ja(IAXIS)-0.5*(B(JAXIS,i, j ,k+1) - B(JAXIS,i, j ,k-1)+ &
                                 B(JAXIS,i,j-1,k+1) - B(JAXIS,i,j-1,k-1))/(z(k+1)-z(k-1))
      Ja(JAXIS) = Ja(JAXIS)+0.5*(B(IAXIS,i, j ,k+1) - B(IAXIS,i, j ,k-1)+ &
                                 B(IAXIS,i,j-1,k+1) - B(IAXIS,i,j-1,k-1))/(z(k+1)-z(k-1))
#endif

      d  = 0.5*(dens(j)+dens(j-1))
      Ba = 0.5*(B(:,i,j,k)+B(:,i,j-1,k))
      Ja  = (hy_hall_parameter/d)*Ja

      Flux(MAGX_FLUX:MAGZ_FLUX,j) = Flux(MAGX_FLUX:MAGZ_FLUX,j)+Ba(JAXIS)*Ja-Ja(JAXIS)*Ba
      Flux(ENER_FLUX,j)     = Flux(ENER_FLUX,j)+Ba(JAXIS)*dot_product(Ja,Ba)-Ja(JAXIS)*dot_product(Ba,Ba)

      ! Add hyper-resistive terms
      hd3 = hy_hyperResistivity/(y(j)-y(j-1))**3
      Flux(MAGX_FLUX:MAGZ_FLUX,j) = &
           Flux(MAGX_FLUX:MAGZ_FLUX,j)+hd3*(B(:,i,j+1,k)-3.*B(:,i,j,k)+3.*B(:,i,j-1,k)-B(:,i,j-2,k))

      speed  = max(speed,abs(Ja(IAXIS)),abs(Ja(JAXIS)),abs(Ja(KAXIS)),8.*hd3)
    end do

#if N_DIM == 3

  case(SWEEP_Z)
    i = i1
    j = i2

    do k = 3, nz-1
      Ja(IAXIS) = 0.5*(B(KAXIS,i,j+1, k ) - B(KAXIS,i,j-1, k )+ &
                       B(KAXIS,i,j+1,k-1) - B(KAXIS,i,j-1,k-1))/(y(j+1)-y(j-1))- &
                      (B(JAXIS,i, j , k ) - B(JAXIS,i, j ,k-1))/(z(k)-z(k-1))
      Ja(JAXIS) =     (B(IAXIS,i, j , k ) - B(IAXIS,i, j ,k-1))/(z(k)-z(k-1))- &
                  0.5*(B(KAXIS,i+1,j, k ) - B(KAXIS,i-1,j, k )- &
                       B(KAXIS,i+1,j,k-1) - B(KAXIS,i-1,j,k-1))/(x(i+1)-x(i-1))
      Ja(KAXIS) = 0.5*(B(JAXIS,i+1,j, k ) - B(JAXIS,i-1,j, k )- &
                       B(JAXIS,i+1,j,k-1) - B(JAXIS,i-1,j,k-1))/(x(i+1)-x(i-1))- &
                  0.5*(B(IAXIS,i,j+1, k ) - B(IAXIS,i,j-1, k )- &
                       B(IAXIS,i,j+1,k-1) - B(IAXIS,i,j-1,k-1))/(y(j+1)-y(j-1))

      d  = 0.5*(dens(k)+dens(k-1))
      Ba = 0.5*(B(:,i,j,k)+B(:,i,j,k-1))
      Ja  = (hy_hall_parameter/d)*Ja

      Flux(MAGX_FLUX:MAGZ_FLUX,k) = Flux(MAGX_FLUX:MAGZ_FLUX,k)+Ba(KAXIS)*Ja-Ja(KAXIS)*Ba
      Flux(ENER_FLUX,k) = Flux(ENER_FLUX,k)+Ba(KAXIS)*dot_product(Ja,Ba)-Ja(KAXIS)*dot_product(Ba,Ba)

      ! Add hyper-resistive terms
      hd3 = hy_hyperResistivity/(z(k)-z(k-1))**3
      Flux(MAGX_FLUX:MAGZ_FLUX,k) = &
           Flux(MAGX_FLUX:MAGZ_FLUX,k)+hd3*(B(:,i,j,k+1)-3.*B(:,i,j,k)+3.*B(:,i,j,k-1)-B(:,i,j,k-2))

      speed  = max(speed,abs(Ja(IAXIS)),abs(Ja(JAXIS)),abs(Ja(KAXIS)),8.*hd3)
    end do

#endif

#endif

  end select

end subroutine hy_8wv_addHallFluxes
