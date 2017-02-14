!!****if* source/Driver/DriverMain/Driver_getNumProcs
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

  use Driver_data, ONLY : dr_meshNumProcs, dr_meshAcrossNumProcs, &
       dr_globalNumProcs,&
       dr_axisNumProcs
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, INTENT(IN) :: communicatorType
  integer, INTENT(OUT) :: numProcs
  integer, optional, intent(IN) :: axis
  logical :: alright

  select case(communicatorType)
  case(GLOBAL_COMM)
     numProcs=dr_globalNumProcs
  case(MESH_COMM)
     numProcs=dr_meshNumProcs
  case(MESH_ACROSS_COMM)
     numProcs=dr_meshAcrossNumProcs
  case(AXIS_COMM)
     alright=present(axis)
     if(alright) alright=(axis.le.MDIM) 
     if(alright) then
        numProcs=dr_axisNumProcs(axis)
     else
        call Driver_abortFlash("Driver_getNumProcs : for directional comm, a valid axis value is needed")
     end if
  end select
end subroutine Driver_getNumProcs
