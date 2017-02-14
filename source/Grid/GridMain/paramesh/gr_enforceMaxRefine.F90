!!***if* source/Grid/common/paramesh/gr_enforceMaxRefine
!!
!! NAME
!!  gr_enforceMaxRefine
!!
!! SYNOPSIS
!!
!!  call gr_enforceMaxRefine(integer(IN) :: lrefineUserMax)
!!  
!!  DESCRIPTION
!!
!!    Modify refinement flags in order to enforce a user-level global maximum
!!    refinement level if at all possible.
!!
!!    A user-level global maximum refinement level may be lower than the global
!!    lrefine_max known to the PARAMESH implementation.
!!
!!
!!  ARGUMENTS 
!!
!!    lrefineUserMax - the maximum current refinement level desired by the
!!                     application.
!! 
!!  NOTES
!!
!!   The current user-level global maximum refinement level is taken to be
!!   the variable gr_refineMax in the Grid_data module.
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_enforceMaxRefine(lrefineUserMax)


  use tree, ONLY: refine,derefine,lrefine_max,lnblocks,&
       lrefine_min,nodetype,lrefine

  implicit none

  integer, intent(in) :: lrefineUserMax

  integer :: i, maxLevel


!==============================================================================
  maxLevel = min(lrefine_max, lrefineUserMax)
  if (maxLevel < 1) maxLevel = 1

  do i = 1, lnblocks
     if(nodetype(i) == 1) then
        if (lrefine(i) > maxLevel) derefine(i) = .TRUE.
     end if

     if (lrefine(i) >= maxLevel) refine(i) = .FALSE.
  end do

#ifdef DEBUG_GRID
  print *, 'exiting gr_enforceMaxRefine ->', maxLevel
#endif
  return
end subroutine gr_enforceMaxRefine
