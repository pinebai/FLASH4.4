!!****if* source/physics/Hydro/HydroMain/split/RHD/Hydro_init
!!
!! NAME
!!
!!  Hydro_init
!!
!!
!! SYNOPSIS
!!
!!  Hydro_init()
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!  PLUTO would be a better choice for doing Relativistic Hydro simulations than FLASH.
!!
!!***

subroutine Hydro_init ()

  use Hydro_data
  use IO_interface,                ONLY : IO_getScalar
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Grid_interface, ONLY:  Grid_setFluxHandling
!!$  use Logfile_interface, ONLY : Logfile_stampMessage
!!$  use Driver_interface,  ONLY : Driver_abortFlash
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs
  implicit none


#include "constants.h"
#include "Flash.h"


  character(len=MAX_STRING_LENGTH) :: eosModeString
  character(len=MAX_STRING_LENGTH) :: hy_str_geometry

  logical   :: restart
!!$  character(len=11) ::  msgGeom
!!$  character(len=100)::  internalFile
!!$  integer   :: imsgDim, hy_meshMe
!!$  real      :: Rconst

!!$  ! icube   = dBaseKeyNumber("xyzCube")
!!$  ! DEV Need to check that cylindrical coordinates can only be in 2D.
!!$  ! DEV  and spherical is only in 1D.  Maybe already checked in geometry init...

  ! Everybody should know these

  call Driver_getMype(MESH_COMM,hy_meshMe)
  call Driver_getNumProcs(MESH_COMM,hy_meshNumProcs)

  call RuntimeParameters_get("geometry", hy_str_geometry)
  call RuntimeParameters_mapStrToInt(hy_str_geometry, hy_meshGeom)

  call RuntimeParameters_get("eosMode", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString, hy_eosMode)

!!$  ! DEV doesn't appear to be ever used
!!$  call RuntimeParameters_get('UnitSystem',units)
  call RuntimeParameters_get('gamma', hy_gamma)
!!$  call PhysicalConstants_get("ideal gas constant", Rconst)
  call RuntimeParameters_get('cfl', hy_cfl)
  call RuntimeParameters_get('irenorm', hy_renorm)
  call RuntimeParameters_get('reconType', hy_reconType)
  call RuntimeParameters_get('flux_correct',hy_fluxCorrect)

  hy_useGravity = .false.
#ifdef GRAVITY
  call RuntimeParameters_get("useGravity", hy_useGravity)
#endif

  call RuntimeParameters_get("restart", restart)

  if (.not.restart) then
     hy_dtmin = HUGE(hy_dtmin)
  else
     call IO_getScalar("dt", hy_dtmin)
  endif
  
  hy_xref = 1.0
  hy_vref = 1.0
  hy_dref = 1.0

  hy_tref = hy_xref/hy_vref
  hy_eref = hy_vref*hy_vref
  hy_nref = hy_dref*hy_vref*hy_xref
  hy_pref = hy_dref*hy_vref*hy_vref
  hy_gref = hy_vref*hy_vref/hy_xref

  !! We assume here that the database knows which system
  !! of units is being used. This is not the case yet.

!!$   hy_qref = hy_vref*hy_vref/Rconst
!!$   hy_kref = hy_dref*hy_vref*hy_xref*Rconst


!!$  call Grid_getGeometry(hy_meshGeom)
!!$  geometryOK = .true.
!!$  if (hy_meshGeom == CYLINDRICAL) then
!!$     if (NDIM > 2) then
!!$        geometryOK = .false.
!!$        msgGeom = "CYLINDRICAL"
!!$        imsgDim = 2
!!$     end if
!!$  else if (hy_meshGeom == SPHERICAL) then
!!$     if (NDIM > 1) then
!!$        geometryOK = .false.
!!$        msgGeom = "SPHERICAL  "
!!$        imsgDim = 1
!!$     end if
!!$     if (.NOT. geometryOK) then
!!$        write(internalFile,900)msgGeom,imsgDim
!!$        call Grid_gethy_meshMe(hy_meshMe)
!!$        call Logfile_stampMessage(internalFile)
!!$        write(*,*) internalFile
!!$        call Driver_abortFlash("[Hydro_init] Bad geometry initialization")
!!$     end if
!!$  end if
!!$900  format('[Hydro_init] Cannot have ',A11,' geometry with dimension greater than ',I3)

  if (NDIM > 1) then
     if (hy_fluxCorrect) then
        call Grid_setFluxHandling('consv_flux_densities')
     end if
  end if

end subroutine Hydro_init
