!!****if* source/Grid/GridSolvers/Multipole/gr_mpoleEvalBinsize
!!
!! NAME
!!
!!  gr_mpoleEvalBinsize
!!
!! SYNOPSIS
!!
!!  gr_mpoleEvalBinsize (real,intent(IN)   : iradius,
!!                       real,intent(INOUT): dr)
!!
!! DESCRIPTION
!!
!!   Evaluates the size of a radial bin size (dr) corresponding to a give
!!   bin number (q)
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleEvalBinsize (q, dr)
  
  !==================================================================
  
  use gr_mpoleData, ONLY : r12, r23, dsinv, scaleType1, scaleType2, scaleType3, &
                           rscale1, rscale2, rscale3,rMax,q1,q2,q3

  implicit none
  
#include "constants.h"
#include "Flash.h"

  integer,intent(IN)     :: q
  real   ,intent(INOUT)  :: dr  

  integer :: func,qloc
  real    :: fact 
  !=====================================================================
  qloc = q

  ! Identify the zone (1,2,3)
  if  (qloc .le. int(q1)) then  
      func  = scaleType1
      fact  = rscale1      
  else if (qloc .le. int(q1+q2)) then
      func   = scaleType2
      fact   = rscale2
      qloc = q - int(q1)
  else
      func  = scaleType3
      fact  = rscale3
      qloc  = q -int (q1+q2)
  endif

  if (func == 1) then      ! Constant Scaling.        
      dr = (1.0/dsinv)*fact
  else if (func == 2) then ! Exponential Scaling.
      dr = ((1.0/(exp(fact)-1))/dsinv)*(exp(qloc*fact)-1.0-exp((qloc-1)*fact)+1.0)
  else                     ! Finest mesh (Default).
      dr = 1.0/dsinv     
  endif
  write (*,*) dr
  pause 
  return
end subroutine gr_mpoleEvalBinsize
