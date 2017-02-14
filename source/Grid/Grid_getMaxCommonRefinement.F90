!!****f* source/Grid/Grid_getMaxCommonRefinement
!!
!! NAME
!!
!!  Grid_getMaxCommonRefinement
!!
!! SYNOPSIS
!!
!!  Grid_getMaxCommonRefinement(integer(IN) :: inputComm, &
!!                              integer(OUT) :: maxRefinement)
!!
!! DESCRIPTION
!!
!!  This is a simple routine to find the maximum common block refinement 
!!  level in the grid.  We may have an AMR grid in which one portion of 
!!  the grid is highly refined.  Here, it may be useful to determine the 
!!  highest block refinement level such that blocks of this level 
!!  completely cover the computational domain.
!!
!! ARGUMENTS
!!
!!  inputComm - Input MPI communicator.
!!  maxRefinement - Max common refinement level of blocks in the 
!!                  inputComm communicator.
!!
!! NOTES
!! 
!!  Communicator argument allows us to compare a subset of processes.
!!  It also makes it explicit to the user that this routine must be 
!!  called by all processes in the passed communicator ...otherwise
!!  deadlock.
!!
!!***

subroutine Grid_getMaxCommonRefinement(inputComm, maxRefinement)

  implicit none
  integer, intent(IN) :: inputComm
  integer, intent(OUT) :: maxRefinement

  maxRefinement = 1

end subroutine Grid_getMaxCommonRefinement
