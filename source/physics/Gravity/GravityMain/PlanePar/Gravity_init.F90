!!****if* source/physics/Gravity/GravityMain/PlanePar/Gravity_init
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
!!  This routine initializes the planePar gravitational physics module.
!!
!! ARGUMENTS
!!
!!  
!!
!!***

subroutine Gravity_init ()

  use Gravity_data
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"

  ! Everybody should know these....
  call Driver_getMype(MESH_COMM,grv_meshMe)
  call Driver_getNumProcs(MESH_COMM,grv_meshNumProcs)

  call RuntimeParameters_get("ptxpos", grv_ptxpos)
  call RuntimeParameters_get("ptmass", grv_ptmass)
  call RuntimeParameters_get("ptdirn", grv_ptdirn)
  call PhysicalConstants_get('Newton', grv_newton)
  call RuntimeParameters_get("useGravity", useGravity)

!==============================================================================

!==============================================================================

  return
end subroutine Gravity_init
