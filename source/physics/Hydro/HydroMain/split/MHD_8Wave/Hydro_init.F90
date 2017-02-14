!!****if* source/physics/Hydro/HydroMain/split/MHD_8Wave/Hydro_init
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
!!  Initialize the 8Wave MHD unit
!!
!!
!! ARGUMENTS
!!
!!
!!***

subroutine Hydro_init()

  use Hydro_data, ONLY : hy_units, hy_killdivb, hy_cfl, hy_Rconst, &
                         hy_dtmin, hy_xref, hy_vref, hy_dref,      &
                         hy_mref,  hy_tref, hy_eref, hy_nref,      &
                         hy_pref,  hy_gref, hy_qref, hy_kref,      &
                         hy_bref,  hy_units,hy_useGravity,         &
                         hy_smallpres,   hy_smalldens,             &
                         hy_fluxCorrect, hy_irenorm,               &
                         hy_meshMe, hy_meshNumProcs,&
                         hy_useDiffuse,hy_useViscosity,&
                         hy_useConductivity,hy_useMagneticResistivity,&
                         hy_maxMagDiff

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs
  use Grid_interface, ONLY:   Grid_setFluxHandling
  use Driver_data, ONLY : dr_dt, dr_restart

  implicit none

#include "constants.h"
#include "Flash.h"

  real :: magResistivity,thermConductivity,kinViscosity

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,hy_meshMe)
  call Driver_getNumProcs(MESH_COMM,hy_meshNumProcs)

  call RuntimeParameters_get('UnitSystem', hy_units)
  call RuntimeParameters_get('killdivb', hy_killdivb)
  call RuntimeParameters_get('smallP', hy_smallpres)
  call RuntimeParameters_get('cfl', hy_cfl)
  call RuntimeParameters_get('flux_correct', hy_fluxCorrect)
  if (NDIM > 1) then
     if (hy_fluxCorrect) then
        call Grid_setFluxHandling('consv_flux_densities')
     end if
  end if

  hy_useGravity = .false.
#ifdef GRAVITY
  call RuntimeParameters_get("useGravity", hy_useGravity)
#endif


  call RuntimeParameters_get("useMagneticResistivity", hy_useMagneticResistivity)
  if (hy_useMagneticResistivity) then
     call RuntimeParameters_get("resistivity", magResistivity)
     hy_maxMagDiff=magResistivity
  endif
  
  call RuntimeParameters_get("useViscosity",    hy_useViscosity)
  if (hy_useViscosity) then
     call RuntimeParameters_get("diff_visc_nu", kinViscosity )
     hy_maxMagDiff=max(hy_maxMagDiff,kinViscosity)
  endif

  call RuntimeParameters_get("useConductivity", hy_useConductivity)
  if (hy_useConductivity) then
     call RuntimeParameters_get("conductivity_constant", thermConductivity)
     hy_maxMagDiff=max(hy_maxMagDiff,thermConductivity)
  endif


  hy_useDiffuse = .false.
  if (hy_useMagneticResistivity .or. hy_useViscosity .or. hy_useConductivity) then
     hy_useDiffuse = .true.
  endif


!!$  ! Hall MHD
!!$  ! It is not supported yet!
!!$  call RuntimeParameters_get("hall_parameter", hy_hall_parameter)
!!$  call RuntimeParameters_get("hyperResistivity", hy_hyperResistivity)

  call RuntimeParameters_get('irenorm', hy_irenorm)
  call PhysicalConstants_get("ideal gas constant", hy_Rconst)

  if (.not.dr_restart) then
     hy_dtmin = HUGE(hy_dtmin)
  else
     hy_dtmin = dr_dt
  endif

  hy_xref = 1.0
  hy_vref = 1.0
  hy_dref = 1.0

  hy_mref = hy_xref*hy_vref
  hy_tref = hy_xref/hy_vref
  hy_eref = hy_vref*hy_vref
  hy_nref = hy_dref*hy_vref*hy_xref
  hy_pref = hy_dref*hy_vref*hy_vref
  hy_gref = hy_vref*hy_vref/hy_xref

  hy_qref = hy_vref*hy_vref/hy_Rconst
  hy_kref = hy_dref*hy_vref*hy_xref*hy_Rconst

  if ( hy_units == "SI" .or. hy_units == "si" ) then
    hy_bref = hy_vref*sqrt(4.0*PI*hy_dref*1.e-7)
  else if ( hy_units == "CGS" .or. hy_units == "cgs" ) then
    hy_bref = hy_vref*sqrt(4.0*PI*hy_dref)
  else
    hy_bref = hy_vref*sqrt(hy_dref)
  end if

end subroutine Hydro_init
