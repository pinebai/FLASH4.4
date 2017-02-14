!!****if* source/Grid/GridSolvers/Multipole/gr_mpoleEvalBin
!!
!! NAME
!!
!!  gr_mpoleEvalBin
!!
!! SYNOPSIS
!!
!!  gr_mpoleEvalBin    (real   ,intent(IN)   : iradius,
!!                      integer,intent(INOUT): q)
!!
!! DESCRIPTION
!!
!!   Evaluates the size of bin number (q) corresponding to a give
!!   location with respect to Center of mass (iradius)
!!
!!***

!!REORDER(4): solnData
subroutine gr_mpoleEvalBin (iradius, q)
  
  !==================================================================
  
  use gr_mpoleData, ONLY : r12, r23, dsinv, scaleType1, scaleType2, scaleType3, &
                           rscale1, rscale2, rscale3,rMax,q1,q2

  implicit none
  
#include "constants.h"
#include "Flash.h"

  real,intent(IN)     :: iradius
  real,intent(INOUT)  :: q

  integer :: func
  real    :: radius, fact 
  !=====================================================================

  if (iradius .le. r12) then
      func   = scaleType1
      fact   = rscale1
      radius = iradius
      q = 0.0
  else if (iradius .gt. r12 .and. iradius .le. r23) then
      func   = scaleType2
      fact   = rscale2
      radius = iradius - r12
      q = q1
  else
      func   = scaleType3
      fact   = rscale3
      radius = iradius - r23
      q = q1 + q2
  endif

  if (func == 1) then      ! Constant Scaling        
     q = q + (radius*dsinv) / fact
  else if (func == 2) then ! Exponential Scaling
     radius = (dsinv*radius)*(exp(fact)-1.0)
     q = q +  log(radius + 1.0)/fact
  else                     ! Finest mesh.
     q = q + radius*dsinv
  endif
  
  return
end subroutine gr_mpoleEvalBin
