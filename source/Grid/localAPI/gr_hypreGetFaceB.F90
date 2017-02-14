!!****if* source/Grid/localAPI/gr_hypreGetFaceB
!!
!!  NAME
!!
!!    gr_hypreGetFaceB
!!
!!  SYNOPSIS
!!
!!    call gr_hypreGetFaceB (integer, intent(IN) :: direction
!!                           integer, intent(IN) :: iFactorB
!!                           integer, intent(IN) :: blkLimits (2,MDIM)
!!                           integer, intent(IN) :: blkLimitsGC (2,MDIM)
!!                           real, intent(IN)    :: solnVec(NUNK_VARS,blkLimitsGC(HIGH,IAXIS),
!!                                                                    blkLimitsGC(HIGH,JAXIS), 
!!                                                                    blkLimitsGC(HIGH,KAXIS))
!!                           real, intent(INOUT) :: flux(NFLUXES,blkLimitsGC(HIGH,IAXIS), 
!!                                                               blkLimitsGC(HIGH,JAXIS), 
!!                                                               blkLimitsGC(HIGH,KAXIS)))
!!
!!  DESCRIPTION
!!   This routine helps to compute iFactorB*AREA on fine-coarse boundaries.
!!
!! ARGUMENTS
!!   direction    : IAXIS/JAXIS/KAXIS (direction along which flux is to be computed).
!!   iFactorB     : Conductivity/Opacitiy.
!!   blkLimitsGC  : an array that holds the lower and upper indices of the
!!                  section of block with the guard cells.
!!   blkLimits    : an array that holds the lower and upper indices of the
!!                  section of block without the guard cells.
!!   solnVec      : SolnVec(:,:,:,:) corresponding to a block.
!!   flux         : iFactorB*Area is stored in the flux vector (done only on the
!!                  outer cells, as they are the ones to be exchanged).
!!
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!***

!!REORDER(4): solnVec


subroutine gr_hypreGetFaceB (direction, iFactorB, blkLimits, blkLimitsGC, solnVec, flux, numVars)
  
  
  implicit none
  
#include "Flash.h"  
#include "constants.h"
  
  integer, intent(IN) :: direction
  integer, intent(IN) :: iFactorB
  integer, intent(IN) :: blkLimits (2,MDIM) 
  integer, intent(IN) :: blkLimitsGC (2,MDIM)
  
  real, intent(IN)    :: solnVec(NUNK_VARS, blkLimitsGC(HIGH,IAXIS), &
                                            blkLimitsGC(HIGH,JAXIS), &
                                            blkLimitsGC(HIGH,KAXIS))   
  
  real, intent(INOUT) :: flux(NFLUXES,blkLimitsGC(HIGH,IAXIS), &
                                      blkLimitsGC(HIGH,JAXIS), &
                                      blkLimitsGC(HIGH,KAXIS))  

  integer, intent(IN) :: numVars

  flux = 0.0
 
  return
  
end subroutine gr_hypreGetFaceB
