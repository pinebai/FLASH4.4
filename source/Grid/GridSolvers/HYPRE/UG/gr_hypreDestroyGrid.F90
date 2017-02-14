!!****if* source/Grid/GridSolvers/HYPRE/UG/gr_hypreDestroyGrid
!!
!!  NAME 
!!
!! gr_hypreDestroyGrid
!!
!!  SYNOPSIS
!!
!!  call gr_hypreDestroyGrid ()
!!
!!
!!  DESCRIPTION 
!! This routine destroys HYPRE Grid. 
!! 
!! 
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  
!! Uses HYPRE library.
!!
!!***

!!REORDER(4): solnVec

#include "Flash.h"

subroutine gr_hypreDestroyGrid ()
  
  use gr_hypreData, ONLY : gr_hypreLower, gr_hypreUpper, &
                           gr_hypreGrid, gr_hypreStencils, &
                           gr_hypreGraph, gr_hypreMatA, gr_hypreVecB, &
                           gr_hypreVecX,  &
                           gr_hypreGridIsSetUp
  
  implicit none
  
#include "HYPREf.h"     
  
  integer :: ierr, s
  
  if (gr_hypreGridIsSetUp) then
     
     deallocate (gr_hypreLower)
     deallocate (gr_hypreUpper)     
     !! HYPRE clean up
     call HYPRE_SStructGridDestroy(gr_hypreGrid, ierr)
     if (allocated(gr_hypreStencils)) then
        do s=0,size(gr_hypreStencils)-1
           call HYPRE_SStructStencilDestroy(gr_hypreStencils(s), ierr)
        end do
        deallocate(gr_hypreStencils)
     end if
     call HYPRE_SStructGraphDestroy(gr_hypreGraph, ierr)  
     call HYPRE_SStructMatrixDestroy(gr_hypreMatA, ierr)  
     call HYPRE_SStructVectorDestroy(gr_hypreVecB, ierr)
     call HYPRE_SStructVectorDestroy(gr_hypreVecX, ierr)     
     
     gr_hypreGridIsSetUp = .FALSE.
     
  end if
  
  
end subroutine gr_hypreDestroyGrid
