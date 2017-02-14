!!****if* source/physics/Hydro/HydroMain/split/MHD_8Wave/divbDiffuse/hy_8wv_sources
!!
!! NAME
!!
!!  hy_8wv_sources
!!
!! SYNOPSIS
!!
!!  hy_8wv_sources(real(IN)    :: Uc(NUNK_VARS,n),
!!                 real(IN)    :: Um(NUNK_VARS,n),
!!                 real(IN)    :: Up(NUNK_VARS,n),
!!                 real(OUT)   :: S(NUNK_VARS,n),
!!                 real(IN)    :: grav(n),
!!                 real(IN)    :: xc(n),
!!                 real(IN)    :: dx(n),
!!                 integer(IN) :: n,
!!                 integer(IN) :: dir)
!!
!!
!! DESCRIPTION
!!
!!  Computes source terms (gravity, etc...) in the MHD equations
!!
!!
!! ARGUMENTS
!!
!!  Uc,Um,Up   - Arrays of interpolated data
!!  S          - Storage array for source terms
!!  grav       - Array containing gravitational acceleration
!!               component in the direction of the sweep
!!  xc         - Array of cell center coordinates
!!  dx         - Array of cell sizes
!!  n          - Size of arrays in the sweep direction
!!  dir        - Sweep direction
!!
!!***

subroutine hy_8wv_sources(Uc,Um,Up,S,grav,xc,dx,n,dir)

  implicit none

#include "Flash.h"
#include "constants.h"

  !!$ Argument list -------------------------------------
  integer, INTENT(IN) :: n,dir
  real, DIMENSION(NUNK_VARS,n), INTENT(OUT) :: S
  real, DIMENSION(NUNK_VARS,n), INTENT(IN) :: Uc,Um,Up
  real, DIMENSION(n), INTENT(IN) :: grav,xc,dx
  !!$ ---------------------------------------------------

  integer :: i,VELN_VAR,MAGN_VAR
  real :: dBdx

  select case(dir)
  case (SWEEP_X)
     VELN_VAR = VELX_VAR
     MAGN_VAR = MAGX_VAR
  case (SWEEP_Y)
     VELN_VAR = VELY_VAR
     MAGN_VAR = MAGY_VAR
  case (SWEEP_Z)
     VELN_VAR = VELZ_VAR
     MAGN_VAR = MAGZ_VAR
  end select

  S = 0.0
  do i = 3, n-2
    dBdx = 0.5*(Up(MAGN_VAR,i)-Um(MAGN_VAR,i)+Um(MAGN_VAR,i+1)-Up(MAGN_VAR,i-1))/dx(i)

    S(VELX_VAR:VELZ_VAR,i) = -dBdx*Uc(MAGX_VAR:MAGZ_VAR,i)
    S(VELN_VAR,i) =  S(VELN_VAR,i)+Uc(DENS_VAR,i)*grav(i)
    S(ENER_VAR,i) = Uc(VELN_VAR,i)*Uc(DENS_VAR,i)*grav(i)- &
      dBdx*dot_product(Uc(VELX_VAR:VELZ_VAR,i),Uc(MAGX_VAR:MAGZ_VAR,i))
    S(MAGX_VAR:MAGZ_VAR,i) = -dBdx*Uc(VELX_VAR:VELZ_VAR,i)

    ! Source term for internal energy
    S(EINT_VAR,i) =((Um(DENS_VAR,i)*Um(EINT_VAR,i)-Up(DENS_VAR,i)*Up(EINT_VAR,i))*Uc(VELN_VAR,i)+ &
                    (Um(DENS_VAR,i)*Um(EINT_VAR,i)+Up(DENS_VAR,i)*Up(EINT_VAR,i)+ &
                     Um(PRES_VAR,i)+Up(PRES_VAR,i))*0.5*(Um(VELN_VAR,i)-Up(VELN_VAR,i)))/dx(i)

  end do

end subroutine hy_8wv_sources
