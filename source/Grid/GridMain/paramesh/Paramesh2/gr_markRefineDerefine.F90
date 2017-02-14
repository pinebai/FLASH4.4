!!***if* source/Grid/common/paramesh/Paramesh2/gr_markRefineDerefine
!!
!! NAME
!!  gr_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  gr_markRefineDerefine(integer(IN) :: iref,
!!                        real(IN) :: refine_cutoff.
!!                        real(IN) :: derefine_cutoff,
!!                        real(IN) :: refine_filter)
!!  
!!  DESCRIPTION
!!  
!!    Blocks are marked for refining or derefining.
!!    This version uses the second derivative calculations on the specified variable to 
!!    determine if the block needs more resoultion (refine) or less resolution (derefine)
!!    de/refine_cutoff are the thresholds for triggering the corresponding action.
!!    This version also takes lrefine_min and lrefine_max into account.
!!    Once the blocks have been marked, the control is passed to Paramesh to update refinement.
!!
!!    If no valid (i.e., positive) index into unk is specified, blocks are only marked for
!!    refinement/derefinement as far as necessary to make their refinement levels lie
!!    between lrefine_min and lrefine_max, inclusively.
!!
!!
!!  ARGUMENTS 
!!
!!
!!    iref - index of the refinement variable in data structure "unk".
!!           May be -1 (actually, any non-positive value) to indicate
!!           that the caller only wants to mark blocks so they will
!!           confirm to the lrefine_min and lrefine_max refinement limits.
!!
!!    refine_cutoff - the threshold value for triggering refinement 
!!
!!    derefine_cutoff - the threshold for triggereing derefinement
!!
!!    refine_filter - makes sure that error calculations to determine refinement
!!                    don't diverge numerically 
!! 
!!  NOTES
!!  
!! ***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_markRefineDerefine(&
                              iref,refine_cutoff,derefine_cutoff,refine_filter)


  use tree, ONLY: refine,derefine,lrefine_max,lnblocks,&
       lrefine_min,nodetype,lrefine
  use Grid_data, ONLY : gr_msgbuffer

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(in) :: iref
  real, intent(in) :: refine_cutoff,derefine_cutoff,refine_filter

  integer :: i

  if (iref > 0) then
     call ref_marking(refine_cutoff,derefine_cutoff,refine_filter,&
                   lrefine_max,gr_msgbuffer,iref)
  end if

!==============================================================================

       
  do i = 1, lnblocks
     if(nodetype(i) == 1) then
        if (lrefine(i) <= lrefine_min) derefine(i) = .FALSE.
         
        if (lrefine(i) < lrefine_min) refine(i) = .TRUE.


        if (lrefine(i) > lrefine_max) derefine(i) = .TRUE.

        if (lrefine(i) >= lrefine_max) refine(i) = .FALSE.
     end if
  end do

#ifdef DEBUG_GRID
  print *, 'exiting mark grid refinement'
#endif
  return
end subroutine gr_markRefineDerefine








