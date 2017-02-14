!!****if* source/physics/Gravity/GravityMain/UserDefined/Gravity_init
!!
!! NAME
!!
!!  Gravity_init
!!
!! SYNOPSIS
!!
!!  Gravity_init()
!!
!! DESCRIPTION
!!
!!  This is the initialization routine for user-defined gravity. 
!!  This implementation of gravity unit is a placeholder for
!!  user-defined functions for the gravitational potential and/or
!!  gravitational acceleration.
!!
!! ARGUMENTS
!!
!!
!!
!!
!!***

#include "constants.h"
subroutine Gravity_init ()

  use Gravity_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs

  implicit none

!==============================================================================

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,grv_meshMe)
  call Driver_getNumProcs(MESH_COMM,grv_meshNumProcs)
  call RuntimeParameters_get("useGravity", useGravity)

!==============================================================================

  return

end subroutine Gravity_init
