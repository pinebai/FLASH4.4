!!****if* source/physics/sourceTerms/Heatexchange/HeatexchangeMain/Immediate/Heatexchange
!!
!! NAME
!!
!!  Heatexchange
!!
!!
!! SYNOPSIS
!!
!!   call Heatexchange ( integer(IN) :: blockCount, 
!!                       integer(IN) :: blockList(blockCount), 
!!                       real(IN)    ::  dt  )    
!!
!! DESCRIPTION
!!
!!  Apply thermal heat exchange among temperature components
!!  to all blocks in the specified list.
!!
!! ARGUMENTS
!!
!!   blockCount -- dimension of blockList
!!   blockList -- array of blocks where components should exchange
!!                internal energy
!!   dt  --       current time step (ignored in this implementation)
!!
!! PARAMETERS
!!
!!  useHeatexchange -- Boolean, True.  Turns on Heatexchange unit
!!  hx_applyToRadiation -- Boolean, False.  Does the Heatexchange unit handle radiation?
!!
!!
!!***


subroutine Heatexchange ( blockCount, blockList, dt )
  
  use Grid_interface, ONLY  : Grid_getBlkIndexLimits, &
       Grid_notifySolnDataUpdate
  use Eos_interface, ONLY   : Eos_wrapped
  use Timers_interface, ONLY : Timers_start, Timers_stop
  
  use Heatexchange_data,ONLY: hx_useHeatexchange, hx_applyToRadiation

  implicit none

#include "constants.h"

  !args
  integer, INTENT(in)                        :: blockCount
  integer, INTENT(in), DIMENSION(blockCount)  :: blockList
  real,    INTENT(in)                        :: dt

  ! locals
  integer                    :: blockID, thisBlock

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  ! Check useHeatexchange flag
  if (.not. hx_useHeatexchange) return

  call Grid_notifySolnDataUpdate(CENTER)

  ! START TIMERS
  call Timers_start("heatXchg")


  ! BEGIN LOOP OVER BLOCKS PASSED IN
  do thisBlock = 1, blockCount
     
     blockID = blockList(thisBlock)
     
     !GET DIMENSION AND COORD POSITIONS
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     

     if (hx_applyToRadiation) then
        ! We assume that the solution variables are thermodynamically
        ! consisten on entry. In particular, we assume that
        !      eint=eion+eele+erad
        ! is already true.
!!$        call Eos_wrapped(MODE_DENS_EI_GATHER,blkLimits,blockID)
        call Eos_wrapped(MODE_DENS_EI_SCATTER,blkLimits,blockID)
     else
        call Eos_wrapped(MODE_DENS_EI_MAT_EQUI,blkLimits,blockID)
        call Eos_wrapped(MODE_DENS_EI_GATHER,blkLimits,blockID)
     end if

  end do
  
  call Timers_stop("heatXchg")
  
  return
  
end subroutine Heatexchange
