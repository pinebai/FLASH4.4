!!****if* source/Driver/DriverMain/Driver_superTimeStep
!!
!! NAME
!!
!!  Driver_superTimeStep
!!
!! SYNOPSIS
!!
!!  call Driver_superTimeStep(
!!                            integer(in) :: nstepSTS)
!!
!! ARGUMENTS
!!
!!  nstepSTS  - the index i for which the substep tau_i is requested
!!              or, if negative two, the sum of all substeps is requested.
!!
!! DESCRIPTION
!!
!! This routine implements the super time steppping advancement algorithm
!! to overcome small diffusive time scales in explicit formulation.
!!
!! REFERENCES
!!
!! "Super-Time-Stepping Acceleration of Explicit Schemes for Parabolic Problems",
!! V. Alexiades, G. Amiez, and P. Gremaud, Com. Num. Meth. Eng, 1996
!! 
!!
!!***

subroutine Driver_superTimeStep(dt,nuSTS,nstepSTS,nstepTotalSTS,dt_subSTS)

  implicit none

#include "constants.h"
  !! Argument list -----------------------
  real, intent(IN)    :: dt,nuSTS
  integer, intent(IN) :: nstepSTS,nstepTotalSTS
  real, intent(OUT)   :: dt_subSTS
  !! -------------------------------------
  integer :: i

  if (nstepSTS .NE. -2) then
     !! Calculate a substep dt_subSTS (tau_i)
     dt_subSTS = dt/((nuSTS - 1.0)*cos((2.*nstepSTS - 1.0)/nstepTotalSTS * 0.5*PI) + 1.0 + nuSTS)
  else
     dt_subSTS = 0.0
     do i=1,nstepTotalSTS
        dt_subSTS = dt_subSTS + dt/((nuSTS - 1.0)*cos((2.*i - 1.0)/nstepTotalSTS * 0.5*PI) + 1.0 + nuSTS)
     end do
  end if

  return
end subroutine Driver_superTimeStep
