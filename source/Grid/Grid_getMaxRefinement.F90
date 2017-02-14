!!****f* source/Grid/Grid_getMaxRefinement
!!
!! NAME
!!
!!  Grid_getMaxRefinement
!!
!! SYNOPSIS
!!
!!  call Grid_getMaxRefinement(integer(OUT) :: maxRefinement,
!!                    OPTIONAL,integer(IN)  :: mode,
!!                    OPTIONAL,integer(IN)  :: scope,
!!                    OPTIONAL,integer(IN)  :: inputComm
!!                              )
!!
!! DESCRIPTION
!!
!!  This routine returns the maximum block refinement level in the grid.  
!!  Depending on the mode used (and modified by other optional arguments),
!!  the returned value represents the maximum refinement level that
!!  is either allowed, or is currently realized anywhere on the grid
!!  (or a subset of the grid).
!!
!!  We may have an AMR grid in which one portion of the grid is highly
!!  refined.  Here, it may be useful to determine globally the highest block
!!  refinement level that occurs, or can potentially occur during the
!!  simulation.
!!
!! ARGUMENTS
!!
!!  maxRefinement - Max common refinement level of blocks in the 
!!                  inputComm communicator.
!!  inputComm - Input MPI communicator, only used if mode=4 and scope=2.
!!              Default - none.
!!  mode      - 1 for lrefine_max,
!!              2 for gr_maxRefine (if used) or lrefine_max,
!!              3 for gr_maxRefine,
!!              4 for existing blocks.
!!              Default is 3.
!!  scope     - scope only used if mode=4;
!!              1 for local to this MPI task,
!!              2 for tasks in inputComm,
!!             [3 for mesh communicator,]
!!             [4 for MPI_COMM_WORLD].
!!              Default is 3.
!!
!! NOTES
!! 
!!  Communicator argument allows us to compare a subset of processes.
!!  It also makes it explicit to the user that this routine must be 
!!  called by all processes in the passed communicator ... otherwise
!!  deadlock. This applies only for modes that require communication.
!!
!!  For a uniform grid implementation, the returned level will always
!!  be 1.
!!
!!  This routine differs from Grid_getMaxCommonRefinement is several ways:
!!   1. Grid_getMaxCommonRefinement looks for existing LEAF blocks with the
!!      smallest refinement level (actual), while
!!      Grid_getMaxRefinement looks for blocks with the highest
!!      refinement level (either actual or potential).
!!   2. Grid_getMaxRefinement has additional optional arguments to select
!!      modes and task subsets.
!!
!!***

subroutine Grid_getMaxRefinement(maxRefinement, mode, scope, inputComm)

  implicit none
  integer, intent(IN), OPTIONAL :: mode, scope
  integer, intent(IN), OPTIONAL :: inputComm
  integer, intent(OUT) :: maxRefinement

  maxRefinement = 0

end subroutine Grid_getMaxRefinement
