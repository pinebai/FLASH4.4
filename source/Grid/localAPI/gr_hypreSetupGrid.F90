!!****if* source/Grid/localAPI/gr_hypreSetupGrid
!!
!!  NAME 
!!
!! gr_hypreSetupGrid
!!
!!  SYNOPSIS
!!
!!  call gr_hypreSetupGrid (integer(IN)  :: blockCount,
!!                          integer(IN), dimension(blockCount) :: blockList,
!!                 OPTIONAL,integer(IN)  :: nvars)
!!
!!  DESCRIPTION 
!! This routine sets up the HYPRE Grid. 
!! Called only once in UG.
!! 
!!
!! ARGUMENTS
!!   blockCount     : The number of blocks in the list.   
!!   blockList      : The list of blocks on which the solution must be updated.
!!   nvars          : Number of variables, also number of equations, for a
!!                    system. Default is 1.
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!
!!
!!***

subroutine gr_hypreSetupGrid (blockCount, blockList, nvars)
  
  implicit none 
  
  integer,                      intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  integer,OPTIONAL,             intent(IN) :: nvars
  
  
  
end subroutine gr_hypreSetupGrid
