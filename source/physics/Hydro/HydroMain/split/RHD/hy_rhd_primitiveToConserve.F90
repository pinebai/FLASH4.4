!!****if* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_primitiveToConserve
!!
!! NAME
!!
!!   hy_rhd_primitiveToConserve
!!
!! SYNOPSIS
!! 
!!  hy_rhd_primitiveToConserve( real (INOUT) :: U(:,:),
!!                              integer (IN) :: ibeg,
!!                              integer (IN) :: iend,
!!                              integer (IN) :: n)
!!
!! DESCRIPTION
!!
!!   Take an array of primitive quantities, U, 
!!   and convert it into an array of conservative
!!   quantities 
!!
!!
!! ARGUMENTS
!!
!!   U -   in input, an array of primitive quantities
!!         on output, an array of conserved quantities
!!   ibeg -    initial point for the conversion
!!   iend -    final   point for the conversion
!!   n    -    number of point in the specified direction
!!
!!
!!***

subroutine hy_rhd_primitiveToConserve (U, ibeg, iend, n)

  use Driver_interface, ONLY : Driver_abortFlash
  use hy_rhd_interface, ONLY: hy_rhd_enthalpy

  implicit none

#include "Flash.h"
#include "RHD.h"

  integer, INTENT(in) :: ibeg, iend, n
  real, DIMENSION(NUNK_VARS,n), INTENT(inout)  :: U

  integer                :: i
  real, DIMENSION(n)     :: h
  real, DIMENSION(NUNK_VARS)  :: up
  real       ::   g, g2, scrh
 
  call hy_rhd_enthalpy (U, h, ibeg, iend, n)

  do i = ibeg, iend
    up(:) = U(:,i)

    !  ----  Get Lorentz Gamma factor  ---- 
   
    g2 = up(VELX_VAR)*up(VELX_VAR) + up(VELY_VAR)*up(VELY_VAR) + up(VELZ_VAR)*up(VELZ_VAR)

    if (up(PRES_VAR) .lt. 0.e0) then 
      print *,' ! Negative pressure in hy_rhd_primitiveToConserve'
      call Driver_abortFlash(" Stop in hy_rhd_primitiveToConserve")
    endif

    if (g2 .gt. 1.e0) then 
      print *,' v > 1 found in P2C',i
      call Driver_abortFlash(" Stop in hy_rhd_primitiveToConserve")
    endif

    g2   = 1.e0/(1.e0 - g2)
    g    = sqrt(g2)
    scrh = up(DENS_VAR)*h(i)*g2

    !  ----  convert U to conservative  ---- 

    U(DENS_VAR,i) = up(DENS_VAR)*g      ! Lab density 
    U(VELX_VAR,i) = up(VELX_VAR)*scrh   ! x-momentum
    U(VELY_VAR,i) = up(VELY_VAR)*scrh   ! y-momentum
    U(VELZ_VAR,i) = up(VELZ_VAR)*scrh   ! z-momentum
    U(ENER_VAR,i) = scrh - up(PRES_VAR) ! total energy density 

  end do

end subroutine hy_rhd_primitiveToConserve
