!!****if* source/physics/Gravity/GravityMain/Poisson/Pfft/Gravity_init
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
!! ARGUMENTS
!!
!!  
!!
!! DESCRIPTION
!!
!!  Initialize the Poisson solver.  Read in any of the
!!  runtime parameters for this solver.  All solver common data
!!  is stored in the mpole_common module
!!
!!***

subroutine Gravity_init()

  use Gravity_data 

  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype,&
       Driver_getNumProcs, Driver_getComm
       
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"

  

  character(len=MAX_STRING_LENGTH) :: strGeometry
  real                             :: newton


  ! Everybody should know these two
  call Driver_getMype(MESH_COMM,grv_meshMe)
  call Driver_getNumProcs(MESH_COMM,grv_meshNumProcs)
  call Driver_getComm(MESH_COMM,grv_meshComm)

  call RuntimeParameters_get("geometry", strGeometry)
  call RuntimeParameters_mapStrToInt(strGeometry, grav_geometry)


  call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)

  call RuntimeParameters_get("useGravity", useGravity)
  call RuntimeParameters_get("updateGravity", updateGravity)
  call RuntimeParameters_get("grav_unjunkPden", grav_unjunkPden)
  call PhysicalConstants_get("Newton", newton)
  grav_poisfact = 4.* PI * newton

  ! abort if we don't have the correct type of BCs.
  select case (grav_boundary_type)

    case ("isolated")
       call Driver_abortFlash('Gravity_init: Isolated boundary conditions are not supported.')
       grav_boundary = ISOLATED
    case ("dirichlet")
       call Driver_abortFlash('Gravity_init: Dirichlet boundary conditions are not supported.')
       grav_boundary = DIRICHLET
    case ("outflow")
       call Driver_abortFlash('Gravity_init: Outflow/Neumann boundary conditions are not supported.')
       grav_boundary = OUTFLOW
    case ("hydrostatic")
       call Driver_abortFlash('Gravity_init: Hydrostatic boundary conditions are not supported.')
       grav_boundary = HYDROSTATIC
    case ("periodic")
       grav_boundary = PERIODIC
    case default
       call Driver_abortFlash('Gravity_init: unrecognized or unsupported gravity boundary type')
  end select

  return
end subroutine Gravity_init
