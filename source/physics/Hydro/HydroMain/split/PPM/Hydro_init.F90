!!****if* source/physics/Hydro/HydroMain/split/PPM/Hydro_init
!!
!! NAME
!!
!!  Hydro_init
!!
!!
!! SYNOPSIS
!!
!!  call Hydro_init()
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
!!  
!!
!! PARAMETERS
!!
!!   These are the runtime parameters used in the split PPM Hydro 
!!   implementation.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have overwritten these values with the flash.par values
!!   for your specific run.  
!!
!!    useHydro [BOOLEAN]
!!        Should any Hydro calculations be performed?
!!    updateHydroFluxes [BOOLEAN]
!!        whether fluxes computed by Hydro should be used to update the solution
!!    cfl [REAL]
!!    geometry [STRING]
!!        Grid geometry
!!    eosMode [STRING]
!!        the default Eos mode, usually MODE_DENS_EI, 
!!        where density and energy are provided to 
!!        calculate pressure and temperature. This is used
!!        for calling Eos on guard cells initially (where necessary).
!!    hy_eosModeAfter [STRING]
!!        the Eos mode, usually MODE_DENS_EI,
!!        where density and energy are provided to 
!!        calculate pressure and temperature. This is used
!!        for calling Eos on interior cells of all leaf blocks
!!        at the end of each Hydro sweep.
!!    irenorm
!!    flux_correct [BOOLEAN]
!!    hybrid_riemann [BOOLEAN]
!!    smlrho [REAL]
!!        Cutoff value for density
!!    smallp [REAL]
!!        Cutoff value for pressure
!!    eintSwitch [REAL]  Defined in Eos Unit
!!    cvisc [REAL]
!!        Artificial viscosity constant
!!    dp_sh_md
!!    epsiln [REAL]
!!    nriem [INTEGER]
!!    omg1
!!        PPM dissipation parameter omega1
!!    omg2
!!        PPM dissipation parameter omega2
!!    rieman_tol [REAL]
!!    small [REAL]
!!       Generic small value that can be used as floor where needed
!!    smallu [REAL]
!!       Cutoff value for velocity
!!    smallx [REAL]
!!       Cutoff value for abundances
!!    vgrid
!!    ppm_modifystates
!!    leveque
!!    igodu [INTEGER]
!!       Enable Guodunov method instead of PPM if set to 1
!!    iplm [INTEGER]
!!       Enable piecewise linear method instead of PPM if set to 1
!!    use_steepening [BOOLEAN]
!!    use_cma_flattening [BOOLEAN]
!!    use_cma_advection [BOOLEAN]
!!    ppmEnerFluxConstructionMeth [INTEGER]
!!    ppmEintFluxConstructionMeth [INTEGER]
!!    ppmEnerCompFluxConstructionMeth [INTEGER]
!!    ppmEnerCompFluxConstructionMeth [INTEGER]
!!***

subroutine Hydro_init()

  !!These are all the runtime parameters.  First the logicals, then the
  !! integers, then the reals    

  use Hydro_data
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype, &
    Driver_getNumProcs
  use Logfile_interface, ONLY : Logfile_stampMessage, &
    Logfile_stampVarMask, Logfile_stamp
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use Grid_interface, ONLY:  Grid_setFluxHandling

  implicit none


  logical :: useSpecialFluxVar
  integer :: istat

#include "constants.h"
#include "Flash.h"  

  character(len=MAX_STRING_LENGTH) :: str_geometry,eosModeString
  logical :: threadBlockListBuild, threadWithinBlockBuild
  
  ! Everybody should know these
  call Driver_getMype(MESH_COMM,hy_meshMe)
  call Driver_getNumProcs(MESH_COMM,hy_meshNumProcs)

  
  call RuntimeParameters_get("useHydro", hy_useHydro)

  !!hydro_timestep
  call RuntimeParameters_get ("cfl", hy_cfl)

  !!**Hydro_sweep RuntimeParameters

  call RuntimeParameters_get ("geometry", str_geometry)
  call RuntimeParameters_mapStrToInt(str_geometry, hy_geometry)

  call RuntimeParameters_get ("eosMode", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString, hy_eosMode)
  if(hy_useHydro) then
    if(hy_eosMode/=MODE_DENS_EI)&
       call Driver_abortFlash("Hydro : Wrong Eos mode for PPM")
  endif
  hy_useGravity = .false.
#ifdef GRAVITY
  call RuntimeParameters_get("useGravity", hy_useGravity)
#endif
  call RuntimeParameters_get("useDiffuse", hy_useDiffuse)
  
  call RuntimeParameters_get("irenorm", hy_irenorm)
  call RuntimeParameters_get("flux_correct",   hy_fluxCorrect)
  call RuntimeParameters_get("hybrid_riemann", hy_hybridRiemann)

  !!**Hydro_updateSolution
  call RuntimeParameters_get("updateHydroFluxes", hy_updateHydroFluxes)
  call RuntimeParameters_get("smlrho",hy_smlrho)
  call RuntimeParameters_get("smallp",hy_smallp)  
  call RuntimeParameters_get("eintSwitch",hy_eintSwitch)  
 
  !!**Hydro_1d needs
  call RuntimeParameters_get("cvisc", hy_cvisc)
  call RuntimeParameters_get("dp_sh_md", hy_dp_sh_md )
  call RuntimeParameters_get("epsiln", hy_epsiln)

  call RuntimeParameters_get("nriem", hy_nriem)
  call RuntimeParameters_get("omg1", hy_omg1)
  call RuntimeParameters_get("omg2", hy_omg2)
  call RuntimeParameters_get("rieman_tol", hy_riemanTol)
  call RuntimeParameters_get("small", hy_small)
  !!smallp from update_soln
  call RuntimeParameters_get("smallu",hy_smallu)
  call RuntimeParameters_get("smallx",hy_smallx)
  !!smlrho from update_soln
  call RuntimeParameters_get("vgrid", hy_vgrid)

  !!**PPM inputs
  call RuntimeParameters_get("ppm_modifystates", hy_ppmModifystates)
  call RuntimeParameters_get("leveque",hy_leveque)
  call RuntimeParameters_get("igodu", hy_igodu)
  call RuntimeParameters_get("iplm", hy_iplm)
  call RuntimeParameters_get("ppmEnerFluxConstructionMeth", hy_ppmEnerFluxConstructionMeth)
  call RuntimeParameters_get("ppmEintFluxConstructionMeth", hy_ppmEintFluxConstructionMeth)
  if (hy_ppmEintFluxConstructionMeth==-1) hy_ppmEintFluxConstructionMeth = hy_ppmEnerFluxConstructionMeth
  call RuntimeParameters_get("ppmEnerCompFluxConstructionMeth", hy_ppmEnerCFluxConstructionMeth)
  call RuntimeParameters_get("ppmEintCompFluxConstructionMeth", hy_ppmEintCFluxConstructionMeth)
  if (hy_ppmEintCFluxConstructionMeth==-1) hy_ppmEintCFluxConstructionMeth = hy_ppmEnerCFluxConstructionMeth
  call RuntimeParameters_get("use_steepening", hy_useSteepening)
  call RuntimeParameters_get("use_cma_flattening", hy_useCmaFlattening)
  call RuntimeParameters_get("use_cma_advection", hy_useCmaAdvection)
  call RuntimeParameters_get("charLimiting", hy_charLimiting) ! new characteristic limiting - DL

  call RuntimeParameters_get("hy_fluxRepresentation", hy_fluxRepresentation)
  if (trim(hy_fluxRepresentation) == "auto") then
#ifdef FLASH_GRID_PARAMESH2
     hy_fluxRepresentation = "hybrid" !Paramesh2 always assumes this - KW
#else
     if (hy_geometry == CARTESIAN) then
        hy_fluxRepresentation = "hybrid"
     else
        hy_fluxRepresentation = "fluxes"
     end if
#endif
  end if

  if (.NOT. hy_useHydro) return ! If Hydto is turned off; return here before anything serious gets done.

  if (hy_fluxRepresentation == "hybrid") then
     hy_useCellAreasForFluxes = .FALSE.
     call Grid_setFluxHandling('consv_flux_densities',status=istat)
  else if (hy_fluxRepresentation == "fluxes") then
     hy_useCellAreasForFluxes = .TRUE.
     call Grid_setFluxHandling('consv_fluxes',status=istat)
  else
     call Driver_abortFlash('Hydro_init: Runtime Parameter hy_fluxRepresentation must be either '//&
          '"hybrid" or "fluxes".')
  end if
#if NDIM > 1
  if (istat .NE. 0) then
     if (hy_fluxCorrect) then
        if (hy_meshMe .EQ. MASTER_PE) print*,'WARNING from Hydro_init: hy_fluxRepresentation '//&
             'was requested as "',hy_fluxRepresentation,'",',' but the Grid unit does not support this,'//&
             ' using the handling supported by Grid.'
        call Logfile_stampMessage('WARNING from Hydro_init: hy_fluxRepresentation '//&
             'was requested as "'//hy_fluxRepresentation//'", but the Grid unit does not support this,'//&
             ' using the handling supported by Grid.')
        hy_useCellAreasForFluxes = .NOT. hy_useCellAreasForFluxes
     end if
  end if
#endif

  useSpecialFluxVar = hy_useCellAreasForFluxes

  if (useSpecialFluxVar) then
     hy_specialFluxVars(1) = P_FLUX
  end if
  if (hy_meshMe .EQ. MASTER_PE) then
     print*,'Info: Hydro_init has set hy_specialFluxVars to ',hy_specialFluxVars
  end if
  

!! Determine some unit-wide variables that are directly derived from
!! Runtime parameters
  if (hy_cvisc == 0.0) then
     hy_transverseStencilWidth = 0
  else
     hy_transverseStencilWidth = 1
  end if

!! Determine the geometries of the individual dimensions

  if (hy_geometry == CARTESIAN)then
     hy_dirGeom(IAXIS) = XYZ
     hy_dirGeom(JAXIS) = XYZ
     hy_dirGeom(KAXIS) = XYZ
  elseif(hy_geometry == POLAR)then
     hy_dirGeom(IAXIS) = RAD_CYL
     hy_dirGeom(JAXIS) = PHI_CYL
     hy_dirGeom(KAXIS) = XYZ
  elseif(hy_geometry == CYLINDRICAL) then
     hy_dirGeom(IAXIS) = RAD_CYL
     hy_dirGeom(JAXIS) = XYZ
     hy_dirGeom(KAXIS) = PHI_CYL
  elseif(hy_geometry == SPHERICAL) then
     hy_dirGeom(IAXIS) = RAD_SPH
     hy_dirGeom(JAXIS) = THETA
     hy_dirGeom(KAXIS) = PHI_SPH
  else
     call Driver_abortFlash("unsupported geometry ")
  end if
  
  if (hy_vgrid /= 0.e0) then
     hy_movingGrid = .true.
  else
     hy_movingGrid = .false.
  endif

!! Now initialize the GC Mask

  hy_gcMask = .FALSE.

  hy_gcMask(PRES_VAR) = .TRUE.

  hy_gcMask(DENS_VAR) = .TRUE.
  hy_gcMask(ENER_VAR) = .TRUE.
#ifdef EINT_VAR
  hy_gcMask(EINT_VAR) = .TRUE.
#endif
  hy_gcMask(TEMP_VAR) = .TRUE.     !for now - only used for initial guess by Helmholtz Eos - KW

#ifdef GAMC_VAR
  hy_gcMask(GAMC_VAR) = .TRUE.
#endif
#ifdef GAME_VAR
  hy_gcMask(GAME_VAR) = .TRUE.
#endif

#ifdef GPOT_VAR
  hy_gcMask(GPOT_VAR) = .TRUE.
#endif
#ifdef GPOL_VAR
  hy_gcMask(GPOL_VAR) = .TRUE.
#endif

#ifdef SGAX_VAR
  hy_gcMask(SGAX_VAR) = .TRUE.  !special meaning for sink particles
#endif
#ifdef SGAY_VAR
  hy_gcMask(SGAY_VAR) = .TRUE.  !special meaning for sink particles
#endif
#ifdef SGAZ_VAR
  hy_gcMask(SGAZ_VAR) = .TRUE.  !special meaning for sink particles
#endif
#ifdef SGXO_VAR
  hy_gcMask(SGXO_VAR) = .TRUE.  !special meaning for sink particles
#endif
#ifdef SGYO_VAR
  hy_gcMask(SGYO_VAR) = .TRUE.  !special meaning for sink particles
#endif
#ifdef SGZO_VAR
  hy_gcMask(SGZO_VAR) = .TRUE.  !special meaning for sink particles
#endif
  
  hy_gcMask(VELX_VAR) = .TRUE.
#if NDIM >= 2
#ifdef VELY_VAR
  hy_gcMask(VELY_VAR) = .TRUE.
#endif
#endif
#if NDIM == 3
#ifdef VELZ_VAR
  hy_gcMask(VELZ_VAR) = .TRUE.
#endif
#endif
#if SPECIES_BEGIN <= UNK_VARS_END
  hy_gcMask(SPECIES_BEGIN:UNK_VARS_END) = .TRUE.
#endif

  hy_gcMaskSize=NUNK_VARS

  call Logfile_stampVarMask(hy_gcMask, .FALSE., '[Hydro_init]', 'gcNeed')

  call RuntimeParameters_get ("hy_eosModeAfter", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString, hy_eosModeAfter)
  if(hy_eosModeAfter/=MODE_DENS_EI .AND. &
       hy_eosModeAfter/=hy_eosMode)&
       call Driver_abortFlash("Hydro : Wrong Eos mode for After PPM Sweep")

  call RuntimeParameters_get("threadBlockListBuild", threadBlockListBuild)
  call RuntimeParameters_get("threadHydroBlockList", hy_threadBlockList)

  call RuntimeParameters_get("threadWithinBlockBuild", threadWithinBlockBuild)
  call RuntimeParameters_get("threadHydroWithinBlock", hy_threadWithinBlock)

  if (hy_threadBlockList .and. .not. threadBlockListBuild) then
     call Logfile_stamp('WARNING! Turning off block list threading '//&
          'because FLASH is not built appropriately','[Hydro_init]')
     hy_threadBlockList = .false.
  end if
  if (hy_threadWithinBlock .and. .not. threadWithinBlockBuild) then
     call Logfile_stamp('WARNING! Turning off within block threading '//&
          'because FLASH is not built appropriately','[Hydro_init]')
     hy_threadWithinBlock = .false.
  end if

end subroutine Hydro_init
