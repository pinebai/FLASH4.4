!!****if* source/physics/sourceTerms/Flame/FlameEffects/EIP/fl_effects
!!
!! NAME
!!
!!  fl_effects
!!
!! SYNOPSIS
!!
!!  call fl_effects(real, dimension(:,:,:,:),POINTER_INTENT_IN  :: solndata,
!!                  real,dimension(:,:,:)(in) :: flamdot,
!!                  real(in) :: dt,
!!                  integer(in) :: blockid)
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!!
!! ARGUMENTS
!!
!!   solndata : 
!!
!!   flamdot : 
!!
!!   dt : 
!!
!!   blockid : ID of block in current processor
!!
!!
!!
!!***


subroutine fl_effects( solnData, flamdot, dt, blockID)

  use fl_effData, only : fl_effDeltae
  use Timers_interface, only : Timers_start, Timers_stop
  use Grid_interface, only : Grid_getBlkIndexLimits
  use Eos_interface, only : Eos_wrapped
  implicit none

#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"
#include "Eos.h"
    
  real, dimension(:,:,:,:),POINTER_INTENT_IN  :: solnData
  real,dimension(:,:,:), intent(in)     :: flamdot
  real,intent(in)                       :: dt
  integer, intent(in)                   :: blockID

  integer                    :: i, j, k
  real                       :: qdot
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC


  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

  ! update interior cells
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
           qdot = flamdot(i,j,k) * fl_effDeltae
           solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) + qdot*dt
#ifdef EINT_VAR
           solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + qdot*dt
#endif

        enddo
     enddo
  enddo
     
  ! We changed internal energy so need to update other quantities
  call Timers_start("eos")
  call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
  call Timers_stop("eos")
     
  return
end subroutine
