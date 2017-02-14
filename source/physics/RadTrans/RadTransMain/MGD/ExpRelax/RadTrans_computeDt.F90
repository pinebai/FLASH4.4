!!****if* source/physics/RadTrans/RadTransMain/MGD/ExpRelax/RadTrans_computeDt
!!
!!  NAME 
!!
!!  RadTrans
!!
!!  SYNOPSIS
!!
!!  call RadTrans_computeDt(integer(IN) :: blockID,
!!                          integer(IN) :: blkLimits(2,MDIM),
!!                          integer(IN) :: blkLimitsGC(2,MDIM),
!!                     real(IN),pointer::  solnData(:,:,:,:),   
!!                     real(OUT)   :: dt_radtrans, 
!!                     real(OUT)   :: dt_minloc(5)) 
!!  DESCRIPTION 
!!    Compute radiative transfer time step
!!
!!  ARGUMENTS
!!    blockID       --  local block ID
!!    blkLimits     --  the indices for the interior endpoints of the block
!!    blkLimitsGC   --  the indices for endpoints including the guardcells
!!    solnData      --  the physical, solution data from grid
!!    dt_radtrans   --  variable to hold timestep constraint
!!    dt_minloc(5)  --  array to hold limiting zone info:  zone indices
!!
!!
!!  J.E.Morel and R.G. McClarren, Stability of Explicit Radiation-Material
!!  Coupling in Radiative Transfer Calculations, Journal of Quantitative Spec-
!!  troscopy and Radiative Transfer, submitted July 2010.
!!
!!
!!***
subroutine RadTrans_computeDt(blockID,  blkLimits,blkLimitsGC, &
     solnData, dt_radtrans, dt_minloc)

#include "constants.h"

  use rt_data, ONLY: rt_useMGD
  use rt_data, ONLY: rt_computeDt
  use rt_data, ONLY: rt_precomputedDt
  use rt_data, ONLY: rt_precomputedMinLoc

  implicit none

  integer, intent(IN) :: blockID
  integer, intent(IN) :: blkLimits(2,MDIM)
  integer, intent(IN) :: blkLimitsGC(2,MDIM)
  real, pointer :: solnData(:,:,:,:) 
  real, intent(INOUT) :: dt_radtrans
  integer, intent(INOUT)  :: dt_minloc(5)
  
  if (.not. rt_useMGD .or. .not. rt_computeDt ) return

  ! The radtrans DT was already computed in RadTrans. Just use the
  ! value computed there...
  dt_radtrans = rt_precomputedDt
  dt_minloc = rt_precomputedMinLoc

  if(dt_radtrans <= 0.0) then
     call Driver_abortFlash("[ RadTrans]: computed dt is not positive! Aborting!")   
  end if

  return 
  
end subroutine RadTrans_computeDt
