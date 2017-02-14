!!****if* source/Driver/DriverMain/Driver_getComm
!!
!! NAME
!!  Driver_getComm
!!
!! SYNOPSIS
!!
!!  Driver_getComm(integer(IN) :: communicatorType,
!!                 integer(OUT):: communicator,
!!       optional, integer(IN) :: axis)
!!               
!!  
!! DESCRIPTION 
!!
!!  All units can use this interface to query the 
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
!! communicator      - output argument in which to return the communicator
!! axis              - this optional argument is included only if 
!!                     communicatorType is AXIS_COMM, the valid 
!!                     values are IAXIS, JAXIS or KAXIS
!!
!!***

#include "constants.h"

subroutine Driver_getComm(communicatorType, communicator, axis)

  use Driver_data, ONLY : dr_meshComm, dr_meshAcrossComm, dr_globalComm,&
       dr_axisComm
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, INTENT(IN) :: communicatorType
  integer, INTENT(OUT) :: communicator
  integer, optional, intent(IN) :: axis
  logical :: alright

  select case(communicatorType)
  case(GLOBAL_COMM)
     communicator=dr_globalComm
#ifdef DEBUG_ALL
     print*,'Driver_getComm: returning communicator=',communicator
#endif
  case(MESH_COMM)
     communicator=dr_meshComm
  case(MESH_ACROSS_COMM)
     communicator=dr_meshAcrossComm
  case(AXIS_COMM)
     alright=present(axis)
     if(alright) alright=(axis.le.MDIM) 
     if(alright) then
        communicator=dr_axisComm(axis)
     else
        call Driver_abortFlash("Driver_getComm : for directional comm, right axis value is needed")
     end if
  case default
     call Driver_abortFlash("Driver_getComm : unrecognized communicatorType")
  end select
end subroutine Driver_getComm
