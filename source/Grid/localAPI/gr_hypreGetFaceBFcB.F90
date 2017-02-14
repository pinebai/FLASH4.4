!!****if* source/Grid/localAPI/gr_hypreGetFaceBFcB
!!
!!  NAME
!!
!!    gr_hypreGetFaceBFcB
!!
!!  SYNOPSIS
!!
!!    call gr_hypreGetFaceBFcB (integer(IN) :: direction,
!!                              integer(IN) :: blkLimits (2,MDIM),
!!                              integer(IN) :: blkLimitsGC (2,MDIM),
!!                    POINTER,  real(IN)    :: facBptr(:,:,:),
!!                              real(INOUT) :: flux(NFLUXES,blkLimitsGC(HIGH,IAXIS), &
!!                                                               blkLimitsGC(HIGH,JAXIS), &
!!                                                               blkLimitsGC(HIGH,KAXIS))
!!
!!  DESCRIPTION 
!!   This routine helps to compute iFactorB*AREA on fine-coarse boundaries. 
!!
!!   This FcB ("face-centered B") variant expects the diffusion coefficients to be
!!   given in the face-centered POINTER array facBptr.
!!
!! ARGUMENTS
!!   direction    : IAXIS/JAXIS/KAXIS (direction along which flux is to be computed).
!!   blkLimitsGC  : an array that holds the lower and upper indices of the
!!                  section of block with the guard cells.
!!   blkLimits    : an array that holds the lower and upper indices of the
!!                  section of block without the guard cells. 
!!   facBptr      : FacBptr(:,:,:) diffusion coefficient for corresponding to a block,
!!                  centered on faces in the given direction.
!!   flux         : iFactorB*Area is stored in the flux vector (done only on the
!!                  outer cells, as they are the ones to be exchanged).
!!   iVar         : if fluxes are being generated for more than one variable,
!!                  iVar is the 0-based number of the current variable.
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!
!!   With UG, there is only this stub.
!!
!! SEE ALSO
!!
!!  gr_hypreGetFaceB
!!***


subroutine gr_hypreGetFaceBFcB (direction, blkLimits, blkLimitsGC, facBptr, flux, iVar)
  
  implicit none
  
#include "constants.h"
#include "FortranLangFeatures.fh"
  
  integer, intent(IN) :: direction
  integer, intent(IN) :: blkLimits (2,MDIM) 
  integer, intent(IN) :: blkLimitsGC (2,MDIM)
  
  real, POINTER_INTENT_IN :: facBptr(:,:,:)
  real, intent(INOUT) :: flux(:,:,:,:)
  integer, intent(IN) :: iVar
end subroutine gr_hypreGetFaceBFcB
