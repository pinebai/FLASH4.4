!!****f* source/Driver/Driver_getMype
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

  implicit none
  integer, INTENT(IN) :: communicatorType
  integer, INTENT(OUT) :: mype
  integer, optional, INTENT(IN) :: axis
  mype=0
end subroutine Driver_getMype
