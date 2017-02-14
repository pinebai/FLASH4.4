!!****if* source/Grid/GridSolvers/HYPRE/paramesh/gr_hypreGetFaceBFcB
!!
!!  NAME
!!
!!    gr_hypreGetFaceBFcB
!!
!!  SYNOPSIS
!!
!!    call gr_hypreGetFaceBFcB (integer, intent(IN) :: direction,
!!                           integer, intent(IN) :: blkLimits (2,MDIM),
!!                           integer, intent(IN) :: blkLimitsGC (2,MDIM),
!!                 POINTER,  real, intent(IN)    :: facBptr(:,:,:),
!!                           real, intent(INOUT) :: flux(NFLUXES,blkLimitsGC(HIGH,IAXIS), &
!!                                                               blkLimitsGC(HIGH,JAXIS), &
!!                                                               blkLimitsGC(HIGH,KAXIS))
!!
!!  DESCRIPTION 
!!   This routine helps to compute iFactorB*AREA on fine-coarse boundaries. 
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
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!
!!***


subroutine gr_hypreGetFaceBFcB (direction, blkLimits, blkLimitsGC, facBptr, flux, iVar)
  
  use Timers_interface, ONLY : Timers_start, Timers_stop  
  
  implicit none
  
#include "Flash.h"  
#include "constants.h"
#include "FortranLangFeatures.fh"
  
  integer, intent(IN) :: direction
  integer, intent(IN) :: blkLimits (2,MDIM) 
  integer, intent(IN) :: blkLimitsGC (2,MDIM)
  
  real, POINTER_INTENT_IN :: facBptr(:,:,:)
  
  real, intent(INOUT) :: flux(:,:,:,:)

  integer, intent(IN) :: iVar

  integer :: BoxLow  (MDIM)
  integer :: BoxHigh (MDIM)  
  
  integer :: i,j,k, ix(4), jy(4), kz(4), ii
  integer :: istep
  integer :: jstep
  integer :: kstep
  integer,parameter :: fpv = 2**(NDIM-1) ! fluxes per variable
  
  call Timers_start("gr_hypreGetFaceBFcB")      
  
  BoxLow (IAXIS)  = blkLimits(LOW,IAXIS)
  BoxLow (JAXIS)  = blkLimits(LOW,JAXIS)
  BoxLow (KAXIS)  = blkLimits(LOW,KAXIS)  
  BoxHigh(IAXIS)  = blkLimits(HIGH,IAXIS)
  BoxHigh(JAXIS)  = blkLimits(HIGH,JAXIS)
  BoxHigh(KAXIS)  = blkLimits(HIGH,KAXIS)
  
  istep = 2
  jstep = 2
  kstep = 2  

  ix = 0
  jy = 0
  kz = 0
  
  select case (direction)
     
  case (IAXIS)
     
     BoxHigh(IAXIS) = blkLimits(HIGH,IAXIS) + 1
     istep = BoxHigh(IAXIS)-BoxLow(IAXIS)               
     jy(2) = K2D
     jy(4) = K2D           
     kz(3) = K3D
     kz(4) = K3D     
     
  case (JAXIS)
     
     BoxHigh(JAXIS) = blkLimits(HIGH,JAXIS) + 1     
     jstep = BoxHigh(JAXIS)-BoxLow(JAXIS)     
     ix(2) = 1     
     ix(4) = 1
     kz(3) = K3D
     kz(4) = K3D     
     
  case (KAXIS)
     
     BoxHigh(KAXIS) = blkLimits(HIGH,KAXIS) + 1     
     kstep = BoxHigh(KAXIS)-BoxLow(KAXIS)          
     ix(2) = 1     
     ix(4) = 1          
     jy(3) = 1
     jy(4) = 1
     
  end select
  
  
  do k = BoxLow(KAXIS), BoxHigh(KAXIS), kstep
     do j = BoxLow(JAXIS), BoxHigh(JAXIS), jstep
        do i = BoxLow(IAXIS), BoxHigh(IAXIS), istep                                
           
           do ii = 1, fpv
              
                 flux(ii+iVar*fpv ,i+ix(ii), j+jy(ii), k+kz(ii))  =  &
                      facBptr(i+ix(ii),j+jy(ii),k+kz(ii))
           end do
           
        end do
     end do
  end do
  
  
  
  call Timers_stop("gr_hypreGetFaceBFcB") 
  
  return
  
end subroutine gr_hypreGetFaceBFcB
