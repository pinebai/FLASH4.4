!!****if* source/Grid/GridMain/paramesh/paramesh4/Paramesh_init
!!
!! NAME
!!  Paramesh_init
!! 
!! SYNOPSIS
!!  call Paramesh_init()
!!
!! DESCRIPTION
!!  Simple interface to mesh initialization routine.
!!  Each mesh package will have a different copy of this function.
!!
!!  This Paramesh3 / Paramesh4 implementation:
!!  does some geometry-related checks and (if possible) readjusts
!!  geometry-related paramesh variables to the geometry requested by
!!  Flash runtime parameter.
!!
!!  Paramesh runtime parameters (these are different from Flash runtime
!!  parameters!) will be dumped to a file.
!!
!! USED BY
!!  Grid_init
!!***

#include "Flash.h"
! Undefine these symbols if defined by Flash.h, we want the ones from constants.h! - KW
#ifdef CARTESIAN
#undef CARTESIAN
#endif
#ifdef SPHERICAL
#undef SPHERICAL
#endif
#ifdef CYLINDRICAL
#undef CYLINDRICAL
#endif
#ifdef POLAR
#undef POLAR
#endif
#include "constants.h"

subroutine Paramesh_init()
  use paramesh_interfaces , ONLY : amr_initialize
  use Driver_interface, ONLY: Driver_abortFlash
  use paramesh_dimensions, ONLY: ndim, l2p5d, nxb, nyb, nzb
  use physicaldata, ONLY: curvilinear, cartesian_pm, cylindrical_pm, &
       spherical_pm, polar_pm, gcell_on_fc, lsingular_line

#ifdef FLASH_GRID_PARAMESH4DEV_SURR_BLKS_OPTION
  use physicaldata, ONLY: use_flash_surr_blks_fill, use_reduced_orrery
#endif

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
!  use physicaldata, ONLY: interp_mask_unk
!  use workspace, ONLY: interp_mask_work
  use Grid_data, ONLY: gr_geometry, gr_kmin, gr_kmax, gr_meshMe


  implicit none
  logical :: testvar


  call amr_initialize()

#ifdef LIBRARY
#define LIBRARY_OR_PM4DEV
#else
#ifdef FLASH_GRID_PARAMESH4DEV
#define LIBRARY_OR_PM4DEV
#endif
#endif

#ifdef LIBRARY_OR_PM4DEV
  if (nxb /= NXB .or. nyb /= NYB .or.  nzb /= NZB) then     
     if(gr_meshMe == MASTER_PE) then
        print*,'PARAMESH  nxb,nyb,nzb are ',nxb, nyb, nzb
        print*,'FLASH NXB,NYB,NZB are ',NXB, NYB, NZB
     end if
     
     call Driver_abortFlash('[Paramesh_init] The value of the nxb, nyb, nzb runtime parameter is  &
          & incompatible with the NXB, NYB, NZB compiled into FLASH.')
  end if
#endif

#ifdef LIBRARY_OR_PM4DEV
  select case (gr_geometry)
     case (CARTESIAN)
        if (curvilinear) then
           cartesian_pm = .true.
        else
           cartesian_pm = .false.
        end if
        spherical_pm = .false.
        cylindrical_pm = .false.
        polar_pm = .false.
     case (SPHERICAL)
        cartesian_pm = .false.
        spherical_pm = .true.
        cylindrical_pm = .false.
        polar_pm = .false.
     case (CYLINDRICAL)
        cartesian_pm = .false.
        spherical_pm = .false.
        cylindrical_pm = .true.
        polar_pm = .false.
     case (POLAR)
        cartesian_pm = .false.
        spherical_pm = .false.
        cylindrical_pm = .false.
        polar_pm = .true.
     end select
  call gr_amr_dump_runtime_parameters
#else
  call gr_amr_dump_runtime_parameters
  select case (gr_geometry)
     case (CARTESIAN)
        testvar = (cartesian_pm .or. .not.curvilinear)
     case (SPHERICAL)
        testvar = spherical_pm
     case (CYLINDRICAL)
        testvar = cylindrical_pm
     case (POLAR)
        testvar = polar_pm
     case default
        testvar = .false.
     end select
     if (.not. testvar) then
        call Driver_abortFlash('[Paramesh_init] The value of the geometry runtime parameter is incompatible    &
             & with the geometry compiled into the AMR code.')
     end if
#endif


  if (cylindrical_pm .or. spherical_pm .or. polar_pm) then
     if (.not. curvilinear) then
        call Driver_abortFlash('[Paramesh_init] Curvilinear support must be enabled in order to use the requested geometry!')
     end if
  end if


  if (ndim==2 .AND. l2p5d==1) then
     gcell_on_fc(KAXIS,:) = .false.

     if (gr_kmax == gr_kmin) then
        if(gr_meshMe == MASTER_PE) print*,&
             'Grid_init: Forcing zmin > zmax to prevent trouble. It should not matter anyway since NDIM=2.'
        if (cylindrical_pm .OR. spherical_pm) then
           gr_kmax = gr_kmin + 180 ! these are in degrees at this point.
           gr_kmin = gr_kmin - 180
        else
           gr_kmax = gr_kmin + 1
        end if
     end if
  end if

#ifdef FLASH_GRID_PARAMESH4DEV_SURR_BLKS_OPTION
  call RuntimeParameters_get("use_flash_surr_blks_fill", &
       use_flash_surr_blks_fill)
  
  if (use_flash_surr_blks_fill .eqv. .true.) then
     !Check the setting is consistent with the geometry.
     if ( lsingular_line .and. NDIM > 1 .and. &
          (spherical_pm .or. polar_pm) ) then
        if (gr_meshMe == 0) then
           print *, "[Paramesh_init]: use_flash_surr_blks_fill "//&
                "option has been reset to false because of the geometry"
        end if
        use_flash_surr_blks_fill = .false.
     end if
  end if

  call RuntimeParameters_get("use_reduced_orrery", &
       use_reduced_orrery)
#endif


!! NOTE: Maybe this is the place to put stuff like this?
!  interp_mask_unk = 2
!  interp_mask_work = 2
!!  call init_flash_physicaldata()
end subroutine Paramesh_init
