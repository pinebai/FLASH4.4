!!****if* source/physics/Gravity/GravityMain/Poisson/Multigrid/Gravity_init
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
!!  Initialize the multipole Poisson solver.  Read in any of the
!!  runtime parameters for this solver.  All solver common data
!!  is stored in the mpole_common module
!!
!!***

subroutine Gravity_init()

  use Gravity_data 

  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype,&
       Driver_getComm, Driver_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use Logfile_interface, ONLY: Logfile_stampMessage
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Multigrid.h"
#include "Flash.h"

  

  character(len=MAX_STRING_LENGTH) :: strGeometry
  real                             :: newton
  integer, save                    :: iimg_mass, iimg_pot

  ! Everybody should know these two
  call Driver_getMype(MESH_COMM,grv_meshMe)
  call Driver_getComm(MESH_COMM,grv_meshComm)
  call Driver_getNumProcs(MESH_COMM,grv_meshNumProcs)

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
       grav_boundary = MG_BND_ISOLATED
    case ("periodic")
       grav_boundary = MG_BND_PERIODIC
    case default
       call Driver_abortFlash('Gravity_init: unrecognized or unsupported gravity boundary type')
  end select


  !specific init routine for the multigrid solver is called gr_hInit, located
  !  in Grid/GridSolvers/Multigrid
  !  it is called from Grid/Grid_initDomain via gr_solversInit()
  

#ifdef IMGM_VAR
  iimg_mass = IMGM_VAR
#else
  iimg_mass = NONEXISTENT
#endif
#ifdef IMGP_VAR
  iimg_pot  = IMGP_VAR
#else
  iimg_pot = NONEXISTENT
#endif

  if (grav_boundary == MG_BND_ISOLATED) then
     ! Check for existence of variables
     if ( (iimg_mass .EQ. NONEXISTENT) .or. (iimg_pot .EQ. NONEXISTENT)) then
        call Logfile_stampMessage(&
             "[Gravity_init] must setup with unit implementation GravityMain/Poisson/Multigrid")
        call Driver_abortFlash( &
             "[Gravity_init] must include isolated boundary subunit implementation IsoBndMultipole")
     end if
  endif

  return
end subroutine Gravity_init
