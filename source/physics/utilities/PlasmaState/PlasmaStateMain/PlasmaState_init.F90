!!****if* source/physics/utilities/PlasmaState/PlasmaStateMain/PlasmaState_init
!!
!! NAME
!!  
!!  PlasmaState_init
!!
!! SYNOPSIS
!! 
!!  call PlasmaState_init ()
!!
!! DESCRIPTION
!!
!!  Perform various initializations for the PlasmaState unit.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

#include "constants.h"

subroutine PlasmaState_init ()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface,            ONLY : Driver_abortFlash, &
                                          Driver_getMype
  use pls_interface,                ONLY : pls_initComposition

  use PlasmaState_data, ONLY: pls_usePlasmaState, &
                              pls_globalMe, pls_meshMe

  implicit none
    
  call Driver_getMype(MESH_COMM,  pls_meshMe)
  call Driver_getMype(GLOBAL_COMM,pls_globalMe)

  call RuntimeParameters_get ("usePlasmaState",   pls_usePlasmaState)

  if (.not.pls_usePlasmaState) then
       return
  end if

  call pls_initComposition()

  return
end subroutine PlasmaState_init
