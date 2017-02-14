!!****if* source/Grid/GridSolvers/HYPRE/paramesh/gr_hypreGetFaceB
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
  
  use Timers_interface, ONLY : Timers_start, Timers_stop  
  
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

  integer :: BoxLow  (MDIM)
  integer :: BoxHigh (MDIM)  
  
  integer :: i,j,k, ix(4), jy(4), kz(4), ii, iv
  integer :: ioffset, istep
  integer :: joffset, jstep
  integer :: koffset, kstep
  integer,parameter :: fpv = 2**(NDIM-1) ! fluxes per variable
  
  call Timers_start("gr_hypreGetFaceB")      
  
  BoxLow (IAXIS)  = blkLimits(LOW,IAXIS)
  BoxLow (JAXIS)  = blkLimits(LOW,JAXIS)
  BoxLow (KAXIS)  = blkLimits(LOW,KAXIS)  
  BoxHigh(IAXIS)  = blkLimits(HIGH,IAXIS)
  BoxHigh(JAXIS)  = blkLimits(HIGH,JAXIS)
  BoxHigh(KAXIS)  = blkLimits(HIGH,KAXIS)
  
  ioffset = 0
  joffset = 0
  koffset = 0

  istep = 2
  jstep = 2
  kstep = 2  

  ix = 0
  jy = 0
  kz = 0
  
  select case (direction)
     
  case (IAXIS)
     
     BoxHigh(IAXIS) = blkLimits(HIGH,IAXIS) + 1
     ioffset = -1     
     istep = BoxHigh(IAXIS)-BoxLow(IAXIS)               
     jy(2) = K2D
     jy(4) = K2D           
     kz(3) = K3D
     kz(4) = K3D     
     
  case (JAXIS)
     
     BoxHigh(JAXIS) = blkLimits(HIGH,JAXIS) + 1     
     joffset = -1     
     jstep = BoxHigh(JAXIS)-BoxLow(JAXIS)     
     ix(2) = 1     
     ix(4) = 1
     kz(3) = K3D
     kz(4) = K3D     
     
  case (KAXIS)
     
     BoxHigh(KAXIS) = blkLimits(HIGH,KAXIS) + 1     
     koffset = -1     
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
              
              do iv = 0,numVars-1
                 flux(ii+iv*fpv ,i+ix(ii), j+jy(ii), k+kz(ii))  =  &
                      0.5*(solnVec(iFactorB+iv,i+ix(ii),j+jy(ii),k+kz(ii)) &
                         + solnVec(iFactorB+iv,i+ix(ii)+ioffset,j+jy(ii)+joffset,k+kz(ii)+koffset))    
              end do
           end do
           
        end do
     end do
  end do
  
  
!!$  do k = BoxLow(KAXIS), BoxHigh(KAXIS), kstep
!!$     do j = BoxLow(JAXIS), BoxHigh(JAXIS), jstep
!!$        do i = BoxLow(IAXIS), BoxHigh(IAXIS), istep                                
!!$           
!!$           flux(1,i+ix(1),j+jy(1),k+kz(1)) = & 
!!$                0.5*(solnVec(iFactorB,i,j,k) + solnVec(iFactorB,i+ioffset,j+joffset,k+koffset))       
!!$           
!!$#if NDIM >= 2
!!$           
!!$           flux(2,i+ix(2),j+jy(2),k+kz(2))  =  &
!!$                0.5*(solnVec(iFactorB,i+ix(2),j+jy(2),k+kz(2)) &
!!$                + solnVec(iFactorB,i+ix(2)+ioffset,j+jy(2)+joffset,k+kz(2)+koffset))                  
!!$           
!!$#if NDIM == 3
!!$           
!!$           flux(3,i+ix(3),j+jy(3),k+kz(3))  =  &
!!$                0.5*(solnVec(iFactorB,i+ix(3),j+jy(3),k+kz(3)) &
!!$                + solnVec(iFactorB,i+ix(3)+ioffset,j+jy(3)+joffset,k+kz(3)+koffset))
!!$           
!!$           flux(4,i+ix(4),j+jy(4),k+kz(4))  =  &
!!$                0.5*(solnVec(iFactorB,i+ix(4),j+jy(4),k+kz(4)) &
!!$                + solnVec(iFactorB,i+ix(4)+ioffset,j+jy(4)+joffset,k+kz(4)+koffset))    
!!$
!!$
!!$           !!flux(3,i,j,k+kz)     = 0.5*(solnVec(iFactorB,i,j,k+kz)     + & 
!!$           !!     solnVec(iFactorB,i+ioffset,j+joffset,k+kz+koffset))           
!!$           !!flux(4,i,j+jy,k+kz)  = 0.5*(solnVec(iFactorB,i,j+jy,k+kz)  + &
!!$           !!     solnVec(iFactorB,i+ioffset,j+jy+joffset,k+kz+koffset))           
!!$#endif
!!$
!!$#endif
!!$           
!!$        end do
!!$     end do
!!$  end do
  
  call Timers_stop("gr_hypreGetFaceB") 
  
  return
  
end subroutine gr_hypreGetFaceB
