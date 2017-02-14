!!****if* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_enthalpy
!! 
!! NAME
!! 
!!  hy_rhd_enthalpy 
!! 
!! SYNOPSIS
!!  hy_rhd_enthalpy( real   (IN) :: u(:,;),
!!                   real  (OUT) :: h(:),
!!                   integer(IN) :: ibeg,
!!                   integer(IN) :: iend,
!!                   integer(IN) :: n)
!! DESCRIPTION
!!
!!  Calculate the enthalpy given pressure and
!!  proper density
!!
!! ARGUMENTS
!! 
!!  u     -  2-D array of primitive values
!!  h     -  1-D array containing enthalpies
!!  ibeg  -  initial point
!!  iend  -  final point
!!  n     -  number of points
!!
!!***

subroutine hy_rhd_enthalpy(u, h, ibeg, iend, n)

  use Hydro_data, ONLY :hy_gamma

  implicit none

#include "Flash.h"
#include "RHD.h"

  !! Argument list ---------------------------------------
  integer, INTENT(in) :: ibeg,iend,n
  real, DIMENSION(NUNK_VARS,n), INTENT(in) :: u
  real, DIMENSION(n), INTENT(out)     :: h
  !! -----------------------------------------------------

  real :: gamma_r
  integer :: i

  gamma_r = hy_gamma/(hy_gamma - 1.d0)
  do i = ibeg, iend
    h(i)  = 1.d0 + gamma_r*u(PRES_VAR,i)/u(DENS_VAR,i)
  end do

end subroutine hy_rhd_enthalpy

