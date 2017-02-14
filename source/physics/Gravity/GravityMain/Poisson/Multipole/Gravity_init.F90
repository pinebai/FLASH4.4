!!****if* source/physics/Gravity/GravityMain/Poisson/Multipole/Gravity_init
!!
!! NAME
!!
!!  Gravity_init
!!
!! 
!! SYNOPSIS
!!
!!  Gravity_init()
!!
!!
!! DESCRIPTION
!!
!!  Initialize the multipole Poisson solver.  Read in any of the
!!  runtime parameters for this solver.  All solver common data
!!  is stored in the mpole_common module
!!
!!  ARGUMENTS
!!
!!  
!!
!!***

subroutine Gravity_init()

  use Gravity_data
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype,&
      Driver_getComm, Driver_getNumProcs

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none
  real,save :: newton

#include "constants.h"

  
  character(len=MAX_STRING_LENGTH) :: strGeometry

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,grv_meshMe)
  call Driver_getComm(MESH_COMM,grv_meshComm)
  call Driver_getNumProcs(MESH_COMM,grv_meshNumProcs)


  call RuntimeParameters_get("geometry", strGeometry)

  call RuntimeParameters_mapStrToInt(strGeometry, grav_geometry)
  !! DEV testing for invalid gravity geometries?  Perhaps it is done in Grid

  
  call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)

  ! Can't use periodic b.c. with Multipole
  select case (grav_boundary_type)
     case ("periodic")
        call Driver_abortFlash('[Gravity_init] No periodic gravity boundary conditions with Multipole.')
     case ("isolated")
        !! Life is good here.
        !! the following variable is not used in multipole implementations, setting it here anyway
        grav_boundary = ISOLATED
     case default
        call Driver_abortFlash('[Gravity_init] Unsupported gravity boundary conditions, only isolated allowed.')
  end select

  call RuntimeParameters_get("useGravity", useGravity)
  call RuntimeParameters_get("updateGravity", updateGravity)
  call RuntimeParameters_get("grav_unjunkPden", grav_unjunkPden)
  call PhysicalConstants_get("Newton", newton)

  grav_poisfact = 4. * PI * Newton


  return
end subroutine Gravity_init
