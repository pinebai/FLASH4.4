!!****if* source/physics/Hydro/HydroMain/split/PPM/multiTemp/Hydro_init
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
!!        For 3T Hydro, MODE_DENS_EI_RECAL_GATHER is recommended.
!!    hy_eosModeAfter [STRING]
!!        the Eos mode, usually MODE_DENS_EI,
!!        where density and energy are provided to 
!!        calculate pressure and temperature. This is used
!!        for calling Eos on interior cells of all leaf blocks
!!        at the end of each Hydro sweep.
!!        For 3T Hydro, MODE_DENS_EI_GATHER is recommended.
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
!!    hy_dbgReconstConsvSele [BOOLEAN]
!!    eos_smallEion [REAL]
!!    eos_smallEele [REAL]
!!    eos_smallErad [REAL]
!!***

subroutine Hydro_init()

  !!These are all the runtime parameters.  First the logicals, then the
  !! integers, then the reals    

  use Hydro_data
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype, &
    Driver_getNumProcs
  use Logfile_interface, ONLY : Logfile_stampMessage, &
    Logfile_stampVarMask
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Eos_interface, ONLY               : Eos_getParameters
  use Grid_interface, ONLY:  Grid_setFluxHandling

  implicit none

  logical :: useSpecialFluxVar
  integer :: istat

#include "constants.h"
#include "Flash.h"  
#include "Hydro_components.h"

  character(len=MAX_STRING_LENGTH) :: str_geometry,eosModeString
  
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
     if(hy_eosMode/=MODE_DENS_EI .AND. hy_eosMode/=MODE_DENS_EI_ALL .AND. &
          hy_eosMode/=MODE_DENS_EI_SCATTER .AND. &
          hy_eosMode/=MODE_DENS_EI_RECAL_GATHER .AND. hy_eosMode/=MODE_DENS_EI_GATHER)&
          call Driver_abortFlash("Hydro : Wrong Eos mode for PPM")
  end if
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
  call Eos_getParameters(smallE1=hy_smallEion, smallE2=hy_smallEele, smallE3=hy_smallErad)
  call RuntimeParameters_get("eintSwitch",hy_eintSwitch)  
  call RuntimeParameters_get("eint1Switch",hy_eint1Switch)  
  call RuntimeParameters_get("eint2Switch",hy_eint2Switch)  
  call RuntimeParameters_get("eint3Switch",hy_eint3Switch)  
  if (hy_eint1Switch==-1.0) hy_eint1Switch = hy_eintSwitch
  if (hy_eint2Switch==-1.0) hy_eint2Switch = hy_eintSwitch
  if (hy_eint3Switch==-1.0) hy_eint3Switch = hy_eintSwitch
 
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

  call PhysicalConstants_get("electron mass",hy_eMass)
  call PhysicalConstants_get("proton mass",hy_pMass)
  call PhysicalConstants_get("electron mass",hy_eMassInUAmu,unitMass="amu")

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
!!$     print*,'size(hy_specialFluxVars,1) is',size(hy_specialFluxVars,1)
!!$     print*,'A', '(size(hy_specialFluxVars,1).LE.', HYDROCOMP_ION,'+1) -->',(size(hy_specialFluxVars,1).LE. HYDROCOMP_ION+1) 
!!$     print*,'B', '(size(hy_specialFluxVars,1).LE.', HYDROCOMP_ELE,'+1) -->',(size(hy_specialFluxVars,1).LE. HYDROCOMP_ELE+1) 
!!$     print*,'C', '(size(hy_specialFluxVars,1).LE.', HYDROCOMP_RAD,'+1) -->',(size(hy_specialFluxVars,1).LE. HYDROCOMP_RAD+1) 
     if (size(hy_specialFluxVars,1).GE. HYDROCOMP_ION+1) hy_specialFluxVars(HYDROCOMP_ION+1) = PION_FLUX
     if (size(hy_specialFluxVars,1).GE. HYDROCOMP_ELE+1) hy_specialFluxVars(HYDROCOMP_ELE+1) = PELE_FLUX
     if (size(hy_specialFluxVars,1).GE. HYDROCOMP_RAD+1) hy_specialFluxVars(HYDROCOMP_RAD+1) = PRAD_FLUX
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
#ifdef EION_VAR
  hy_gcMask(EION_VAR) = .TRUE.
#endif
#ifdef EELE_VAR
  hy_gcMask(EELE_VAR) = .TRUE.
#endif
#ifdef ERAD_VAR
  hy_gcMask(ERAD_VAR) = .TRUE.
#endif
#ifdef PION_VAR
  hy_gcMask(PION_VAR) = .TRUE.
#endif
#ifdef PELE_VAR
  hy_gcMask(PELE_VAR) = .TRUE.
#endif
#ifdef PRAD_VAR
  hy_gcMask(PRAD_VAR) = .TRUE.
#endif
  hy_gcMask(TEMP_VAR) = .TRUE.     !for now - used for initial guess by Eos Newton-Raphson loop - KW
#ifdef TION_VAR
  hy_gcMask(TION_VAR) = .TRUE.     !for now - only used for initial guess by Helmholtz Eos - KW
#endif
#ifdef TELE_VAR
  hy_gcMask(TELE_VAR) = .TRUE.     !for now - only used for initial guess by Helmholtz Eos - KW
#endif
#ifdef TRAD_VAR
  hy_gcMask(TRAD_VAR) = .TRUE.     !Now too - only used for initial guess by Helmholtz Eos - KW
#endif

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

  ! For testing ways to advect components and handle shock heating
  
  call RuntimeParameters_get ("hy_eosModeAfter", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString, hy_eosModeAfter)
  if(hy_eosModeAfter/=MODE_DENS_EI_SELE_GATHER .AND. &
       hy_eosModeAfter/=MODE_DENS_EI_SHOCKSELE_GATHER .AND. &
       hy_eosModeAfter/=MODE_DENS_EI_SCATTER .AND. &
       hy_eosModeAfter/=MODE_DENS_EI_GATHER .AND. &
       hy_eosModeAfter/=hy_eosMode)&
       call Driver_abortFlash("Hydro : Wrong Eos mode for After PPM Sweep")

#ifndef SELE_MSCALAR
  if(hy_eosModeAfter == MODE_DENS_EI_SELE_GATHER) then
     call Logfile_stampMessage('[Hydro_init] ERROR:')
     call Logfile_stampMessage('[Hydro_init] You set hy_eosModeAfter to "dens_ie_sele_gather"')
     call Logfile_stampMessage('[Hydro_init] To use this mode, the electron entropy mass scalar')
     call Logfile_stampMessage('[Hydro_init] must exist. Make sure the following line exists in')
     call Logfile_stampMessage('[Hydro_init] your simulation Config file:')
     call Logfile_stampMessage('[Hydro_init] MASS_SCALAR sele EOSMAP: SELE')
     call Driver_abortFlash("[Hydro_init] Must create SELE_MSCALAR, see log file for details")
  end if
#endif

  call RuntimeParameters_get ("hy_3Ttry_Arelated", hy_3Ttry_Arelated)
  call RuntimeParameters_get ("hy_3Ttry_useShockDetect", hy_3Ttry_useShockDetect)
  call RuntimeParameters_get ("hy_3Ttry_B", hy_3Ttry_B)
  call RuntimeParameters_get ("hy_3Ttry_B_rad", hy_3Ttry_B_rad)
  if (hy_3Ttry_B_rad < 0) hy_3Ttry_B_rad = hy_3Ttry_B
  call RuntimeParameters_get ("hy_3Ttry_D", hy_3Ttry_D)
  call RuntimeParameters_get ("hy_3Ttry_E", hy_3Ttry_E)
  call RuntimeParameters_get ("hy_3Ttry_F", hy_3Ttry_F)
  call RuntimeParameters_get ("hy_3Ttry_G", hy_3Ttry_G)
  call RuntimeParameters_get ("hy_3Ttry_Q", hy_3Ttry_Q)

  call RuntimeParameters_get ("hy_dbgReconstConsvSele", hy_dbgReconstConsvSele)
end subroutine Hydro_init
