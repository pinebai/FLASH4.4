!!****if* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_checkBoundaryValues
!!
!!
!! NAME
!!
!!  hy_rhd_checkBoundaryValues
!!
!! SYNOPSIS
!!  hy_rhd_checkBoundaryValues( real   (IN) :: Vc(:),
!!                              integer(IN) :: n,
!!                              integer(IN) :: dir)
!!
!! DESCRIPTION
!!
!!  Check whether the boundary values for the current block 
!!  are set consistently.
!!
!!
!!  ARGUMENTS
!!
!!    Vc  -         1-D array of primitive values, center value
!!    n   -         number of points in a block including guardcells
!!    dir -         sweep direction
!!      
!!
!!***
subroutine hy_rhd_checkBoundaryValues (Vc, n, dir)

  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  
#include "Flash.h"
#include "constants.h"
#include "RHD.h"

  integer, INTENT(IN)                  :: n, dir
  real, DIMENSION(NUNK_VARS,n), INTENT(IN) :: Vc

  integer  :: i,k
  real     :: v2L, v2R

  k = n-NGUARD

  ! Check boundary values
  do i = 1,NGUARD

    v2L =   Vc(VELX_VAR,i)*Vc(VELX_VAR,i) &
          + Vc(VELY_VAR,i)*Vc(VELY_VAR,i) &
          + Vc(VELZ_VAR,i)*Vc(VELZ_VAR,i)
    if (v2L .ge. 1.e0) then
      print *," ! V > 1 in left boundary",v2L, i
      call Driver_abortFlash(" ! stop in hy_rhd_checkBoundaryValues")
    end if

    v2R =   Vc(VELX_VAR,i+k)*Vc(VELX_VAR,i+k) &
          + Vc(VELY_VAR,i+k)*Vc(VELY_VAR,i+k) &
          + Vc(VELZ_VAR,i+k)*Vc(VELZ_VAR,i+k)
    if (v2R .ge. 1.e0) then
      print *," V >1 in right boundary",v2R,i
      call Driver_abortFlash(" ! stop in hy_rhd_checkBoundaryValues")
    end if

    if (Vc(PRES_VAR,i) .lt. 0.e0) then
      print *," ! Negative pressure in left boundary !"
      call Driver_abortFlash(" ! stop in hy_rhd_checkBoundaryValues")
    end if
    if (Vc(PRES_VAR,i+k) .lt. 0.e0) then
      print *," ! Negative pressure in right boundary !"
      call Driver_abortFlash(" ! stop in hy_rhd_checkBoundaryValues")
    end if
    if (Vc(DENS_VAR,i) .lt. 0.e0) then
      print *," ! Negative density in left boundary !"
      call Driver_abortFlash(" ! stop in hy_rhd_checkBoundaryValues")
    end if
    if (Vc(DENS_VAR,i+k) .lt. 0.e0) then
      print *," ! Negative density in right boundary !"
      call Driver_abortFlash(" ! stop in hy_rhd_checkBoundaryValues")
    end if

  end do   
end subroutine hy_rhd_checkBoundaryValues
