!!****if* source/Driver/DriverMain/Driver_getMype
!!
!! NAME
!!  Driver_getMype
!!
!! SYNOPSIS
!!
!!  Driver_getMype(integer(IN) :: communicatorType,
!!                 integer(OUT):: mype,
!!       optional, integer(IN) :: axis)
!!               
!!  
!! DESCRIPTION 
!!
!!  All units can use this interface to query the value of 
!!  mype for communicators in use. At initialization all units call this interface
!!  to get the global mype. Some units may call to get the meshComm mype, 
!!  only units such as UG and Grid solvers for UG are likely to want axisComm
!! 
!!
!! ARGUMENTS
!!
!!  communicatorType - Input argument indicating whether to return the global 
!!                     communicator (GLOBAL_COMM), or mesh communicator
!!                     (MESH_COMM) that allows duplication of mesh or
!!                     it is communicator that handles rows and colums of
!!                     processors (AXIS_COMM)
!! mype              - output argument in which to return the mype value
!! axis              - this optional argument is included only if 
!!                     communicatorType is AXIS_COMM, the valid 
!!                     values are IAXIS, JAXIS or KAXIS
!!
!!***

#include "constants.h"

subroutine Driver_getMype(communicatorType, mype, axis)

  use Driver_data, ONLY : dr_meshMe, dr_meshAcrossMe, dr_globalMe,&
       dr_axisMe
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, INTENT(IN) :: communicatorType
  integer, INTENT(OUT) :: mype
  integer, optional, intent(IN) :: axis
  logical :: alright

  select case(communicatorType)
  case(GLOBAL_COMM)
     mype=dr_globalMe
  case(MESH_COMM)
     mype=dr_meshMe
  case(MESH_ACROSS_COMM)
     mype=dr_meshAcrossMe
  case(AXIS_COMM)
     alright=present(axis)
     if(alright) alright=(axis.le.MDIM) 
     if(alright) then
        mype=dr_axisMe(axis)
     else
        call Driver_abortFlash("Driver_getMype : for directional comm, a valid axis value is needed")
     end if
  end select
end subroutine Driver_getMype
