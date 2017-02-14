!!****if* source/Grid/localAPI/gr_hypreGridStatus
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
!!   Stub implementation.
!!
!!***


subroutine gr_hypreGridStatus (blockCount, blockList, nvars)
  
  implicit none  
  
  integer, intent(IN):: blockCount
  integer, dimension(blockCount),intent(IN):: blockList
  integer, intent(IN),OPTIONAL :: nvars

  
  return
  
end subroutine gr_hypreGridStatus
