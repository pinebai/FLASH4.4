!!****f* source/Grid/Grid_updateRefinement
!!
!! NAME
!!
!!  Grid_updateRefinement
!!
!!
!! SYNOPSIS
!!  
!!  call Grid_updateRefinement(integer(IN) :: nstep,
!!                             real(IN)    :: time,
!!                    OPTIONAL,logical(OUT):: gridChanged)
!!
!!
!! DESCRIPTION
!!
!!  This routine is applicable only to mesh packages that use adaptive grid.
!!
!!  Applies user-defined refinment critera to determine which blocks need 
!!  to be refined and derefined.  Once the blocks are marked, calls
!!  amr_refine_derefine to actually carry out the refinements.  During this
!!  stage, the blocks are redistributed across processors (if needed).  
!!
!!  After the refinement, the newly created child blocks are filled via
!!  prolongation from the coarse parents.  This prolongation step can use
!!  prolongation routines from paramesh, or defined by the user
!!
!!  Once the prolongation is done, the guardcells are filled.  Finally, the
!!  EOS is called on the block interiors to make them thermodynamically
!!  consistent.
!!
!!
!! ARGUMENTS
!!
!!  nstep : current step number
!!  time  : current evolution time
!!  gridChanged : returns TRUE if grid may actually have changed.
!!
!!***


subroutine Grid_updateRefinement( nstep,time, gridChanged)

  implicit none

  integer, intent(in) :: nstep
  real, intent(in) :: time
  logical,intent(out),OPTIONAL :: gridChanged

  if (present(gridChanged)) gridChanged = .FALSE.
  return
end subroutine Grid_updateRefinement
