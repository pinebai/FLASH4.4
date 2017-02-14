!!****if* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_conserveToPrimitive
!!
!!
!! NAME
!!
!!   hy_rhd_conserveToPrimitive
!!
!! SYNOPSIS
!! 
!!  hy_rhd_conserveToPrimitive( real (INOUT) :: U(:,:),
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
!!   U    -   in input, an array of conserved quantities
!!            on output, an array of primitive quantities
!!   ibeg -   initial point for the conversion
!!   iend -   final   point for the conversion
!!   n    -   number of point in the specified direction
!!
!!
!!***

subroutine hy_rhd_conserveToPrimitive (U, ibeg, iend, n)

  use hy_rhd_interface, ONLY : hy_rhd_pressureFunc
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "RHD.h"

  integer, INTENT(in)                      :: ibeg, iend, n
  real, DIMENSION(NUNK_VARS,n), INTENT(inout)  :: U

  integer, PARAMETER :: max_iter = 50
  real, PARAMETER    :: acc      = 1.e-12
  real, PARAMETER    :: beta_max = 0.99d0

  integer               :: nv, i, iter, iflag
  real, DIMENSION(NUNK_VARS) :: uc
  real                  :: vx, vy, vz, p, m2
  real                  :: scrh, dp, yp, dyp

  do i = ibeg, iend

    uc(:) = U(:,i)

    if (uc(DENS_VAR) < 0.e0) then 
      print *,'Negative density in C2P'
      call Driver_abortFlash(" Stop in hy_rhd_conserveToPrimitive")
    end if

    if (uc(ENER_VAR) < 0.e0) then 
      print *,' !Negative energy in C2P ,',U(ENER_VAR,i)
      call Driver_abortFlash(" Stop in hy_rhd_conserveToPrimitive")
    end if
      
    m2 = uc(VELX_VAR)*uc(VELX_VAR) + uc(VELY_VAR)*uc(VELY_VAR) + uc(VELZ_VAR)*uc(VELZ_VAR)

    !  ------------------------------------------
    !   solve numerically the implicit equation 
    !   to recover pressure from conservative 
    !   quantities
    !  ------------------------------------------  

    p = sqrt(m2) - uc(ENER_VAR)  !     Provide initial guess   
    p = max(p, 1.e-18)

    iflag = 0
    do iter = 0, max_iter
    
      call hy_rhd_pressureFunc(uc, p, yp, dyp, iflag)
      if (iflag .eq. 1) then 
        print *," ! Bad iteration in C2P at ",i
      call Driver_abortFlash(" Stop in hy_rhd_conserveToPrimitive")
      end if

      dp = yp/dyp
      p  = p - dp

      if (abs (dp) .lt. acc) exit

    end do

    if(iter .eq. max_iter) then 
      !print *," ! Can not find pressure in C2P"
      p = 1.e-9
    end if

    if (p .lt. 0.e0) then
      !print *," ! Negative pressure found in C2P: ", p, i
      p = 1.e-12
    end if

    ! ----  Find velocity components  ----   

    scrh = 1.e0/(uc(ENER_VAR) + p)
    vx   = uc(VELX_VAR)*scrh
    vy   = uc(VELY_VAR)*scrh
    vz   = uc(VELZ_VAR)*scrh

    ! ----  check v^2 < 1 -----  

    scrh = vx*vx + vy*vy + vz*vz

    if (scrh .ge. 1.e0) then
        
      print *,'Vel > 1 in C2P',i
      call Driver_abortFlash(" Stop in hy_rhd_conserveToPrimitive")
      scrh = beta_max/sqrt(scrh)
      vx = vx*scrh
      vy = vy*scrh
      vz = vz*scrh
      scrh = beta_max*beta_max

    end if

    ! ---- obtain primitive quantities  ----  

    scrh = sqrt(1.e0 - scrh)  ! = 1/lorentz_factor

    ! ----  overwrite U  ----

    U(DENS_VAR,i) = uc(DENS_VAR)*scrh
    U(VELX_VAR,i) = vx
    U(VELY_VAR,i) = vy
    U(VELZ_VAR,i) = vz
    U(PRES_VAR,i) = p

  end do

end subroutine hy_rhd_conserveToPrimitive
