!!****f* source/Driver/Driver_getNumProcs
!!
!! NAME
!!  Driver_getNumProcs
!!
!! SYNOPSIS
!!
!!  Driver_getNumProcs(integer(IN) :: communicatorType,
!!                 integer(OUT):: numProcs,
!!       optional, integer(IN) :: axis)
!!               
!!  
!! DESCRIPTION 
!!
!!  All units can use this interface to query the number of processors for
!!  communicators in use. At initialization all units call this interface.
!!  One can ask for either global communicator, or mesh communicators
!!  or axis communicators used in some implementations of Uniform mesh
!!  When there are no duplicate copies of the mesh, the mesh communicator
!!  defaults to MPI_COMM_WORLD
!! 
!!
!! ARGUMENTS
!!
!!  communicatorType - Input argument indicating whether to return the global 
!!                     communicator (GLOBAL_COMM), or mesh communicator
!!                     (MESH_COMM) that allows duplication of mesh or
!!                     it is communicator that handles rows and colums of
!!                     processors (AXIS_COMM)
!! numProcs          - output argument in which to return the count of procs
!! axis              - this optional argument is included only if 
!!                     communicatorType is AXIS_COMM, the valid 
!!                     values are IAXIS, JAXIS or KAXIS
!!
!!***

#include "constants.h"

subroutine Driver_getNumProcs(communicatorType, numProcs, axis)
  implicit none
  integer, INTENT(IN) :: communicatorType
  integer, INTENT(OUT) :: numProcs
  integer, optional, INTENT(IN) :: axis
  numProcs = 1
end subroutine Driver_getNumProcs
