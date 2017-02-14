!!****if* source/Grid/GridSolvers/HYPRE/paramesh/gr_hypreGridStatus
!!
!!  NAME 
!!
!!  gr_hypreGridStatus
!!
!!  SYNOPSIS
!!
!!  call gr_hypreGridStatus (integer(IN)  :: blockCount,
!!                           integer(IN), dimension(blockCount) :: blockList,
!!                  OPTIONAL,integer(IN)  :: nvars)
!!
!!  DESCRIPTION 
!!      With AMR mesh verifies if the grid has been modified, if
!!      modified it resets the HYPRE grid object. If called first time 
!!      it sets up the HYPRE grid. 
!!
!! ARGUMENTS
!!
!!   blockCount     : The number of blocks in the list.   
!!   blockList      : The list of blocks on which the solution must be updated.
!!   nvars          : Number of variables, also number of equations, for a
!!                    system. Default is 1.
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!   HYPRE grid is setup only once in UG. 
!!
!!***

#include "Flash.h"

subroutine gr_hypreGridStatus (blockCount, blockList, nvars)
  
  use gr_hypreData,     ONLY : gr_hypreGridIsSetUp, &
                               gr_hypreNVars,       &
                               gr_hypreNStep
  use gr_interface,     ONLY : gr_hypreSetupGrid
  use Driver_interface, ONLY : Driver_getNStep

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use tree, ONLY: grid_changed
  
  implicit none
  
#include "constants.h" 

  integer, intent(IN):: blockCount
  integer, dimension(blockCount),intent(IN):: blockList
  integer, intent(IN),OPTIONAL :: nvars

    integer :: NStep
  logical :: forceSetup

  call Timers_start ("gr_hypreGridStatus")

  
  if (.not. gr_hypreGridIsSetUp) then 
     call gr_hypreSetupGrid (blockCount, blockList, nvars)
     
     call Driver_getNStep(gr_hypreNStep)
  else
     if (present(nvars)) then
        forceSetup = ( gr_hypreNVars .NE. nvars )
     else
        forceSetup = ( gr_hypreNVars .NE. 1 )
     end if
     if (grid_changed == 1 .AND. .NOT. forceSetup) then
        call Driver_getNStep(NStep)
        forceSetup = (Nstep /= gr_hypreNStep)
     end if
     if (forceSetup) then
        call gr_hypreDestroyGrid ()
        call gr_hypreSetupGrid (blockCount, blockList, nvars)
        gr_hypreNStep = NStep
        
     end if
  end if
 
  
  call Timers_stop ("gr_hypreGridStatus")

  
  return
  
end subroutine gr_hypreGridStatus
