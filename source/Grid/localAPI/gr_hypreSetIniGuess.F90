!!****if* source/Grid/localAPI/gr_hypreSetIniGuess
!!
!!  NAME 
!!
!!  gr_hypreSetIniGuess
!!
!!  SYNOPSIS
!!
!!  call gr_hypreSetIniGuess (integer,intent(IN) :: iVar,
!!                            integer,intent(IN) :: blockCount,
!!                            integer,intent(IN) :: blockList (blockCount))
!!
!!  DESCRIPTION 
!!      This routine sets initial guess for AX=B.
!!
!! ARGUMENTS
!!   iVar       : Variable on which the diffusion operatorion is performed (e.g TEMP_VAR)
!!   blockCount : The number of blocks in the list.   
!!   blockList  : The list of blocks on which the solution must be updated.  
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!
!!   Uses HYPRE library.
!!
!!***

!!REORDER(4): solnVec


subroutine gr_hypreSetIniGuess (iVar, blockCount, blockList)
  
  implicit none
  
#include "Flash.h"  
#include "constants.h"
  
  integer,intent(IN) :: iVar
  integer,intent(IN) :: blockCount
  integer,intent(IN) :: blockList (blockCount)
  

  return
  
end subroutine gr_hypreSetIniGuess
