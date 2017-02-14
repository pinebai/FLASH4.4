!!****if* source/Grid/GridMain/Chombo/AMR/Grid_init
!!
!! NAME
!!  Grid_init
!!
!! SYNOPSIS
!!
!!  Grid_init()
!!           
!!
!! DESCRIPTION
!!  Initialize Grid_data
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS 
!!
!!  lrefine_max [INTEGER] 
!!      maximum AMR refinement level
!!  lrefine_min [INTEGER] 
!!      minimum AMR refinement level
!!  nrefs [INTEGER] 
!!      refine/derefine AMR grid every nrefs timesteps
!!
!!  refine_var_1 [INTEGER] 
!!     indicates first variable on which to refine
!!  refine_cutoff_1 [REAL] 
!!      threshold value to trigger refinement for refine_var_1
!!  derefine_cutoff_1 [REAL]
!!      threshold value to trigger derefinement for refine_var_1
!!  refine_filter_1 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_1
!!
!!  refine_var_2 [INTEGER] 
!!     indicates second variable on which to refine
!!  refine_cutoff_2 [REAL] 
!!      threshold value to trigger refinement for refine_var_2
!!  derefine_cutoff_2 [REAL]
!!      threshold value to trigger derefinement for refine_var_2
!!  refine_filter_2 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_2
!!
!!  refine_var_3 [INTEGER] 
!!     indicates third variable on which to refine (if needed)
!!  refine_cutoff_3 [REAL] 
!!      threshold value to trigger refinement for refine_var_3
!!  derefine_cutoff_3 [REAL]
!!      threshold value to trigger derefinement for refine_var_3
!!  refine_filter_3 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_3
!!
!!  refine_var_4 [INTEGER] 
!!     indicates fourth variable on which to refine (if needed)
!!  refine_cutoff_4 [REAL] 
!!      threshold value to trigger refinement for refine_var_4
!!  derefine_cutoff_4 [REAL]
!!      threshold value to trigger derefinement for refine_var_4
!!  refine_filter_4 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_4
!!
!!  flux_correct [BOOLEAN]
!!     turns on or off flux correction
!! small  [REAL]
!!   Generic small value that can be used as floor where needed
!! smlrho [REAL]  
!!   Cutoff value for density    
!! smallp [REAL]  
!!   Cutoff value for pressure
!! smalle [REAL]  
!!   Cutoff value for energy
!! smallt [REAL]  
!!   Cutoff value for temperature
!! smallu [REAL]  
!!   Cutoff value for velocity
!! smallx [REAL]  
!!   Cutoff value for abundances
!! eosMode[STRING]
!!   determines which variables to calculate from the ones
!!   defined. Possible values are "dens_ie", "dens_pres" and "dens_temp"
!! interpol_order [INTEGER]
!!   Order of interpolation, used in Paramesh2 "monotonic" interpolation
!!   for mesh prolongation
!! grid_monotone_hack [BOOLEAN]
!!   If .true., apply radical monotonicity constraints to interpolants,
!!   i.e., completely flatten them if they violate monotonicity.  Used
!!   in Paramesh2 "quadratic_cartesian" interpolation for mesh prolongation.
!! earlyBlockDistAdjustment [BOOLEAN]
!!   If .true., let Paramesh redistribute blocks
!!   across processors early, so that the block distribution chosen by
!!   Paramesh will be in effect when time evolution begins after restart.
!!   If earlyBlockDistAdjustment is .false., the block distribution enacted
!!   by the IO unit when it read a checkpoint file will normally still be
!!   in effect when time evolution begins after a restart.
!!   This flag is ignored if not restarting from a checkpoint.
!! 
!! lrefine_del [INTEGER]
!! gr_lrefineMaxRedDoByTime [BOOLEAN]
!! gr_lrefineMaxRedTRef [REAL]
!! gr_lrefineMaxRedTimeScale [REAL]
!! gr_lrefineMaxRedLogBase [REAL]
!!
!! gr_lrefineMaxRedDoByLogR [BOOLEAN]
!! gr_lrefineMaxRedRadiusFact [REAL]
!! x_refine_center [REAL]
!! y_refine_center [REAL]
!! z_refine_center [REAL]
!!
!! gr_restrictAllMethod [INTEGER]
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_init()

  use Grid_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype, &
    Driver_getNumProcs, Driver_getComm
  use Logfile_interface, ONLY : Logfile_stampMessage, Logfile_stamp
  use Simulation_interface, ONLY : Simulation_mapStrToInt, &
    Simulation_getVarnameType
  use Grid_interface, ONLY : Grid_setFluxHandling
  implicit none

  include "Flash_mpi.h"

  integer :: i, j, k, localNumBlocks, ii, numLeafBlks

  character(len=MAX_STRING_LENGTH),save :: refVarname,refVarString,paramString
  character(len=MAX_STRING_LENGTH),save :: refCutoffName,refCutOffString
  character(len=MAX_STRING_LENGTH),save :: derefCutoffName,derefCutOffString
  character(len=MAX_STRING_LENGTH),save :: refFiltername,refFilterString
  character(len=MAX_STRING_LENGTH) :: xl_bcString,xr_bcString
  character(len=MAX_STRING_LENGTH) :: yl_bcString,yr_bcString
  character(len=MAX_STRING_LENGTH) :: zl_bcString,zr_bcString
  character(len=MAX_STRING_LENGTH) :: eosModeString, grav_boundary_type
  real :: dx, dy, dz
  real, dimension(NDIM) :: rnb
  integer,save :: refVar
  integer :: istat
  character(len=MAX_STRING_LENGTH) :: fluxRepresentation

  character(len=*), parameter :: meshInterpConsvdMsg = "WARNING!!! "//&
       "convertToConsvdInMeshInterp is not yet implemented for Chombo"
  character(len=*), parameter :: meshCallsConsvdMsg = "WARNING!!! "//&
       "convertToConsvdForMeshCalls is set to .false.. "//&
       "Our experience shows that it needs to be .true. for good conservation"
  character(len=*), parameter :: meshConsvdExclusiveMsg = "WARNING!!! "//&
       "convertToConsvdForMeshCalls and convertToConsvdInMeshInterp are "//&
       "exclusive parameters.  We have only implemented primitive "//&
       "to conservative conversions for convertToConsvdForMeshCalls, so it "//&
       "is selected over convertToConsvdInMeshInterp"
  character(len=*), parameter :: meshRefRatioMsg = "WARNING!!! "//&
       "We have barely tested refinement jumps greater than 2"
  character(len=*), parameter :: meshTagRadiusMsg = "WARNING!!! "//&
       "Small values of tagRadius have led to non-convergence in Riemann "//&
       "solver.  We think this is because it can lead to fine-coarse interfaces "//&
       "in regions of the domain where there are steep gradients"
  character(len=*), parameter :: meshBlockSizeMsg = "WARNING!!! "//&
       "Small values of BRMeshRefineBlockFactor have led to worse than "//&
       "expected conservation of mass and total energy"


!----------------------------------------------------------------------------------
! mesh geometry - moved here so Paramesh_init can use gr_geometry for some checking
!----------------------------------------------------------------------------------
  call RuntimeParameters_get("geometry",gr_str_geometry)
  call RuntimeParameters_mapStrToInt(gr_str_geometry, gr_geometry)
  call RuntimeParameters_get("geometryOverride",gr_geometryOverride)

  call Driver_getMype(GLOBAL_COMM, gr_globalMe)
  call Driver_getNumProcs(GLOBAL_COMM, gr_globalNumProcs)
  call Driver_getComm(GLOBAL_COMM, gr_globalComm)

  call Driver_getMype(MESH_COMM, gr_meshMe)
  call Driver_getNumProcs(MESH_COMM, gr_meshNumProcs)
  call Driver_getComm(MESH_COMM, gr_meshComm)

  call Driver_getMype(MESH_ACROSS_COMM, gr_meshAcrossMe)
  call Driver_getNumProcs(MESH_ACROSS_COMM, gr_meshAcrossNumProcs)
  call Driver_getComm(MESH_ACROSS_COMM, gr_meshAcrossComm)


#ifdef GRID_WITH_MONOTONIC
  if (NGUARD < 4) then
     if (myPE==MASTER_PE) then
        print*,'Grid_init: Monotonic grid interpolation requires at least 4 layers of guard cells.'
        print*,' However, NGUARD is only ', NGUARD
        print*," Maybe you want to setup with '-gridinterpolation=native',"
        print*," or make sure that NGUARD is set correctly in Config file."
        call Driver_abortFlash("Please setup with '-gridinterpolation=native', or change NGUARD.")
     end if
  endif
#endif


! The following renaming was done: "conserved_var" -> "convertToConsvdForMeshCalls". - KW
  call RuntimeParameters_get("convertToConsvdForMeshCalls", gr_convertToConsvdForMeshCalls)
  call RuntimeParameters_get("convertToConsvdInMeshInterp", gr_convertToConsvdInMeshInterp)
  call RuntimeParameters_get("enableMaskedGCFill", gr_enableMaskedGCFill)

  call RuntimeParameters_get("nrefs", gr_nrefs)
  call RuntimeParameters_get('lrefine_min', lrefine_min)
  call RuntimeParameters_get('lrefine_max', lrefine_max)

  call RuntimeParameters_get("smalle",gr_smalle)
  call RuntimeParameters_get("smlrho",gr_smallrho)
  call RuntimeParameters_get("smallx",gr_smallx) !
!  call RuntimeParameters_get("grid_monotone_hack", gr_monotone) ! for "quadratic_cartesian" interpolation
  call RuntimeParameters_get("interpol_order",gr_intpol) ! for "monotonic" interpolation
#ifdef GRID_WITH_MONOTONIC
  gr_intpolStencilWidth = 2     !Could possibly be less if gr_intpol < 2  - KW
#endif


  !get the boundary conditions stored as strings in the flash.par file
  call RuntimeParameters_get("xl_boundary_type", xl_bcString)
  call RuntimeParameters_get("xr_boundary_type", xr_bcString)
  call RuntimeParameters_get("yl_boundary_type", yl_bcString)
  call RuntimeParameters_get("yr_boundary_type", yr_bcString)
  call RuntimeParameters_get("zl_boundary_type", zl_bcString)
  call RuntimeParameters_get("zr_boundary_type", zr_bcString)

  !map the string boundary conditions to integer constants defined in constants.h
  call RuntimeParameters_mapStrToInt(xl_bcString,gr_domainBC(LOW,IAXIS))
  call RuntimeParameters_mapStrToInt(xr_bcString,gr_domainBC(HIGH,IAXIS))
  call RuntimeParameters_mapStrToInt(yl_bcString,gr_domainBC(LOW,JAXIS))
  call RuntimeParameters_mapStrToInt(yr_bcString,gr_domainBC(HIGH,JAXIS))
  call RuntimeParameters_mapStrToInt(zl_bcString,gr_domainBC(LOW,KAXIS))
  call RuntimeParameters_mapStrToInt(zr_bcString,gr_domainBC(HIGH,KAXIS))

  call RuntimeParameters_get("bndPriorityOne",gr_bndOrder(1))
  call RuntimeParameters_get("bndPriorityTwo",gr_bndOrder(2))
  call RuntimeParameters_get("bndPriorityThree",gr_bndOrder(3))

  !get the initial grid layout
  call RuntimeParameters_get("refine_on_particle_count",gr_refineOnParticleCount)

  call RuntimeParameters_get("min_particles_per_blk",gr_minParticlesPerBlk)
  call RuntimeParameters_get("max_particles_per_blk",gr_maxParticlesPerBlk)

!------------------------------------------------------------------------------
! mesh geometry       (gr_geometry was already set above)
!------------------------------------------------------------------------------
  !get the physical domain limits
  call RuntimeParameters_get('xmin', gr_imin)
  call RuntimeParameters_get('xmax', gr_imax)
  call RuntimeParameters_get('ymin', gr_jmin)
  call RuntimeParameters_get('ymax', gr_jmax)
  call RuntimeParameters_get('zmin', gr_kmin)
  call RuntimeParameters_get('zmax', gr_kmax)

  !Store computational domain limits in a convenient array.  Used later in Grid_getBlkBC.
  gr_globalDomain(LOW,IAXIS) = gr_imin
  gr_globalDomain(LOW,JAXIS) = gr_jmin
  gr_globalDomain(LOW,KAXIS) = gr_kmin
  gr_globalDomain(HIGH,IAXIS) = gr_imax
  gr_globalDomain(HIGH,JAXIS) = gr_jmax
  gr_globalDomain(HIGH,KAXIS) = gr_kmax


! Determine the geometries of the individual dimensions, and scale
! angle value parameters that are expressed in degrees to radians.
! This call must be made after gr_geometry, gr_domainBC, and gr_{j,k}{min,max}
! have been set based on the corresponding runtime parameters.
  call gr_initGeometry()


  call RuntimeParameters_get("eosMode", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString, gr_eosMode)

  call RuntimeParameters_get("eosModeInit", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString, gr_eosModeInit)

  gr_eosModeNow = gr_eosModeInit ! may change after initialization is done

  call RuntimeParameters_get("earlyBlockDistAdjustment", gr_earlyBlockDistAdjustment)
  gr_justExchangedGC = .false.


  !! This section of the code identifies the variables to used in
  !! the refinement criterion. If a variable is a refinement variable
  !! then the corresponding refinement/derefinement cutoff and filter
  !! values also have to be fetched. The config file defines
  !! refinement variables as strings, names as "refine_var_1",
  !! "refine_var_2" etc, with the current maximum being 4. The
  !! general utility routine takes the base "refine_var_" and appends
  !! the index at the end of the string to generate the parameter
  !! name and the routine Simulation_mapStrToInt finds its index into UNK.

  call RuntimeParameters_get("refine_var_count",gr_numRefineVarsMax)
  gr_refine_var = NONEXISTENT
  gr_numRefineVars=0

  refVarName='refine_var_'
  refCutoffName='refine_cutoff_'
  derefCutoffName='derefine_cutoff_'
  refFilterName='refine_filter_'

  do i = 1,gr_numRefineVarsMax
     call concatStringWithInt(refVarName,i,refVarString)
     call RuntimeParameters_get( refVarString, paramString)
     if(paramString /= "none") then
        call Simulation_mapStrToInt(paramString, refVar, MAPBLOCK_UNK)
        if (refVar > 0) then
           gr_numRefineVars=gr_numRefineVars+1
           gr_refine_var(gr_numRefineVars)=refVar
           call concatStringWithInt(refCutoffName,gr_numRefineVars,refCutoffString)
           call concatStringWithInt(derefCutoffName,gr_numRefineVars,derefCutOffString)
           call concatStringWithInt(refFilterName,gr_numRefineVars,refFilterString)
           call RuntimeParameters_get( refCutoffString, gr_refine_cutoff(gr_numRefineVars)  )
           call RuntimeParameters_get( derefCutoffString, gr_derefine_cutoff(gr_numRefineVars) )
           call RuntimeParameters_get( refFilterString,  gr_refine_filter(gr_numRefineVars) )
        else
           if (gr_meshMe==MASTER_PE) print*, 'WARNING: Unrecognized variable name in refine_var_',i,' treating it as "none"'
           call Logfile_stampMessage('WARNING: Unrecognized variable name in refine_var, treating it as "none"')
           
        end if
     end if
  end do

  gr_enforceMaxRefinement = .FALSE.

  call RuntimeParameters_get("lrefine_del", gr_lrefineDel)
  gr_maxRefine=lrefine_max

  call RuntimeParameters_get("gr_lrefineMaxRedDoByLogR", gr_lrefineMaxRedDoByLogR)
  call RuntimeParameters_get("gr_lrefineMaxRedRadiusFact", gr_lrefineMaxRedRadiusSq)
  gr_lrefineMaxRedRadiusSq = gr_lrefineMaxRedRadiusSq * gr_lrefineMaxRedRadiusSq
  call RuntimeParameters_get("x_refine_center", gr_lrefineCenterI)
  call RuntimeParameters_get("y_refine_center", gr_lrefineCenterJ)
  call RuntimeParameters_get("z_refine_center", gr_lrefineCenterK)

  call RuntimeParameters_get("gr_lrefineMaxRedDoByTime", gr_lrefineMaxRedDoByTime)
  if (gr_lrefineMaxRedDoByTime) gr_enforceMaxRefinement = .TRUE.
  call RuntimeParameters_get("gr_lrefineMaxRedTimeScale", gr_lrefineMaxRedTimeScale)
  call RuntimeParameters_get("gr_lrefineMaxRedLogBase", gr_lrefineMaxRedLogBase)
  call RuntimeParameters_get("gr_lrefineMaxRedTRef", gr_lrefineMaxRedTRef)

  if(gr_numRefineVars==0)then
     if(gr_meshMe == MASTER_PE) print*,'WARNING : Adaptive Grid did not find any refinement variables'
     call Logfile_stampMessage("WARNING : Adaptive Grid did not find any variable to refine")
  end if

#ifdef FLASH_PARTICLES
  call RuntimeParameters_get('useParticles',gr_useParticles)
  call RuntimeParameters_get('pt_maxPerProc',gr_maxParticlesPerProc)
#else
  gr_useParticles=.false.
#endif

  gr_guard = NGUARD
  gr_guard(JAXIS) = gr_guard(JAXIS)*K2D
  gr_guard(KAXIS) = gr_guard(KAXIS)*K3D

  gr_allPeriodic = .true.
  do i = 1,NDIM
     if(gr_domainBC(LOW,i)/=PERIODIC)gr_allPeriodic=.false.
     if(gr_domainBC(HIGH,i)/=PERIODIC)gr_allPeriodic=.false.
  end do

  !Check if there are gravitational isolated boundary conditions
  !in order to determine which solvers to intialize.
  call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
  gr_isolatedBoundaries = (grav_boundary_type=="isolated")

  gr_anyVarToConvert = .FALSE.
  do i = UNK_VARS_BEGIN,UNK_VARS_END
     gr_vars(i)=i
     call Simulation_getVarnameType(i, gr_vartypes(i))
     if (gr_vartypes(i) .eq. VARTYPE_PER_MASS) gr_anyVarToConvert = .TRUE.
  end do


  call gr_setDataStructInfo()
  call gr_bcInit()


  !Only call the particle initialization routines when
  !we are using particles.
  if (gr_useParticles .eqv. .true. ) then
     call gr_ptInit()
     call gr_ptMapInit()
  endif

  call RuntimeParameters_get("iGridSize", gr_gIndexSize(IAXIS))
  call RuntimeParameters_get("jGridSize", gr_gIndexSize(JAXIS))
  call RuntimeParameters_get("kGridSize", gr_gIndexSize(KAXIS))
  call RuntimeParameters_get("maxBlockSize", gr_maxBlockSize)
  call RuntimeParameters_get("BRMeshRefineFillRatio", gr_BRMeshRefineFillRatio)
  call RuntimeParameters_get("BRMeshRefineBufferSize", gr_BRMeshRefineBufferSize)
  call RuntimeParameters_get("BRMeshRefineBlockFactor", gr_BRMeshRefineBlockFactor)
  call RuntimeParameters_get("tagRadius", gr_tagRadius)
  call RuntimeParameters_get("verbosity", gr_verbosity)
  call RuntimeParameters_get("QuadCFInterp", gr_useQuadCFInterp)
  call RuntimeParameters_get("flux_correct", gr_useFluxCorrect)
  call RuntimeParameters_get("refRatio", gr_refRatio)
  call RuntimeParameters_get("restrictBeforeGhostExchange", &
       gr_restrictBeforeGhostExchange)

  
  !------------------------- Set the flux handling ---------------------------
  !DEV CD: Code needs to be re-evaluated so that Grid_setFluxHandling is not called
  !both here and in Hydro_init.  I added the call here because Grid_initDomain
  !is called before Hydro_init.  My Chombo initialization function needs to know
  !how we wish to handle the fluxes (gr_scaleFineFluxes).

  !Testing for the existance of property fluxes is a hack so that
  !EOS unit test still works.
#if NPROP_FLUX > 0 && defined(FLASH_HYDRO_PPM)
  call RuntimeParameters_get("hy_fluxRepresentation", fluxRepresentation)
#else
  fluxRepresentation = "fluxes"  !Just so gr_scaleFineFluxes gets initialized.
#endif

  if (trim(fluxRepresentation) == "auto") then
     if (gr_geometry == CARTESIAN) then
        fluxRepresentation = "hybrid"
     else
        fluxRepresentation = "fluxes"
     end if
  end if

  !After calling Grid_setFluxHandling we obtain gr_scaleFineFluxes.
  if (fluxRepresentation == "hybrid") then
     call Grid_setFluxHandling('consv_flux_densities',status=istat)
  else if (fluxRepresentation == "fluxes") then
     call Grid_setFluxHandling('consv_fluxes',status=istat)
  else
     call Driver_abortFlash('Grid_init: Runtime Parameter hy_fluxRepresentation must be either '//&
          '"hybrid" or "fluxes".')
  end if
  !------------------------- Set the flux handling ---------------------------


  !Ensure that all unused dimensions have a domain size of 1.
  do i = NDIM+1, MDIM
     gr_gIndexSize(i) = 1
  end do


  !! Now find the minimum cell size (gr_minCellSize)
  !! - only considering non-angle coordinates

  gr_minCellSizes(IAXIS) = (gr_imax - gr_imin) / &
       (gr_gIndexSize(IAXIS) * gr_refRatio**(lrefine_max-1))
  gr_minCellSize = gr_minCellSizes(IAXIS)

  if (NDIM >= 2) then
     gr_minCellSizes(JAXIS) = (gr_jmax - gr_jmin) / &
          (gr_gIndexSize(JAXIS) * gr_refRatio**(lrefine_max-1))
     if (.not.gr_dirIsAngular(JAXIS)) then
        gr_minCellSize = min(gr_minCellSize,gr_minCellSizes(JAXIS))
     end if
  end if

  if (NDIM == 3) then
     gr_minCellSizes(KAXIS) = (gr_kmax - gr_kmin) / &
          (gr_gIndexSize(KAXIS) * gr_refRatio**(lrefine_max-1))
     if (.not. gr_dirIsAngular(KAXIS)) then
        gr_minCellSize = min(gr_minCellSize,gr_minCellSizes(KAXIS))
     end if
  end if



  ! Display warning messages for features that are either 
  ! 1) not yet implemented
  ! 2) implemented but lead to bad things, e.g. worse than expected mass and
  !    energy conservation.
  ! For the time being write to the screen and the FLASH logfile.
  if (lrefine_max > 1) then
     if (gr_useFluxCorrect) then
        if (gr_convertToConsvdInMeshInterp) then
           call Logfile_stamp(meshInterpConsvdMsg, '[Grid_init]')
           if (gr_meshMe==MASTER_PE) print *, meshInterpConsvdMsg

           if (gr_convertToConsvdForMeshCalls) then
              gr_convertToConsvdInMeshInterp = .false.
              call Logfile_stamp(meshConsvdExclusiveMsg, '[Grid_init]')
              if (gr_meshMe==MASTER_PE) print *, meshConsvdExclusiveMsg
           end if
        end if

        if (.not.gr_convertToConsvdForMeshCalls) then
           call Logfile_stamp(meshCallsConsvdMsg, '[Grid_init]')
           if (gr_meshMe==MASTER_PE) print *, meshCallsConsvdMsg
        end if
     end if

     if (gr_tagRadius < 2) then
        call Logfile_stamp(meshTagRadiusMsg, '[Grid_init]')
        if (gr_meshMe==MASTER_PE) print *, meshTagRadiusMsg
     end if

     if (gr_refRatio > 2) then
        call Logfile_stamp(meshRefRatioMsg, '[Grid_init]')
        if (gr_meshMe==MASTER_PE) print *, meshRefRatioMsg
     end if

     if (gr_BRMeshRefineBlockFactor < 8) then
        call Logfile_stamp(meshBlockSizeMsg, '[Grid_init]')
        if (gr_meshMe==MASTER_PE) print *, meshBlockSizeMsg
     end if
  end if

  gr_region=0.0
end subroutine Grid_init
