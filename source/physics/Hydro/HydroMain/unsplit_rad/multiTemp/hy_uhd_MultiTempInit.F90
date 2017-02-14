!!****if* source/physics/Hydro/HydroMain/unsplit_rad/multiTemp/hy_uhd_MultiTempInit
!!
!! NAME
!!
!!  hy_uhd_MultiTempInit
!!
!!
!! SYNOPSIS
!!
!!  call hy_uhd_MultiTempInit()
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize some unit scope variables which are typically taken
!!  directly from the runtime parameters and which are specific
!!  to the multiTemp variant of the unsplit Hydro implementation.
!!  This must be called once by Hydro_init.F90.  Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!
!!   These are the runtime parameters used in the multiTemp variant of the
!!   the split PPM Hydro inmplementation, in addition to those used in the
!!   standard implementation.
!!
!!DEV:no    ppmEnerFluxConstructionMeth [INTEGER]
!!DEV:no    ppmEintFluxConstructionMeth [INTEGER]
!!DEV:no    ppmEnerCompFluxConstructionMeth [INTEGER]
!!DEV:no    ppmEnerCompFluxConstructionMeth [INTEGER]
!!    eos_smallEion [REAL]
!!    eos_smallEele [REAL]
!!    eos_smallErad [REAL]
!! DEV: List of PARAMETERS is out of date / incomplete.
!!***

subroutine hy_uhd_MultiTempInit()

  !!These are all the runtime parameters.  First the logicals, then the
  !! integers, then the reals    

  use hy_uhd_MultiTempData, ONLY: hy_3Ttry_B, hy_3Ttry_D, hy_3Ttry_E, hy_3Ttry_F, hy_3Ttry_G, &
       hy_3Ttry_B_rad
!!$  use hy_uhd_MultiTempData, ONLY: hy_ppmEnerFluxConstructionMeth, &
!!$       hy_ppmEintFluxConstructionMeth, &
!!$       hy_ppmEnerCFluxConstructionMeth, &
!!$       hy_ppmEintCFluxConstructionMeth
  use hy_uhd_MultiTempData, ONLY: hy_smallEion,hy_smallEele,hy_smallErad
  use hy_uhd_MultiTempData, ONLY: hy_mtScaleAccel, hy_mtScaleWork, hy_mtScaleLorentz
  use hy_uhd_MultiTempData, ONLY: hy_mtFactorB0,hy_mtFactorB1,hy_mtFactorB2
  use hy_uhd_MultiTempData, ONLY: hy_mtFactorBrad0,hy_mtFactorBrad1,hy_mtFactorBrad2
  use hy_uhd_MultiTempData, ONLY: hy_mtPresRatLambda3Min
  use hy_uhd_MultiTempData, ONLY : hy_smoothIterations, &
                                   hy_smoothFluxLim0,   &
                                   hy_smoothFluxLim1,   &
                                   hy_smoothMethod,     &
                                   hy_smoothCoeff,      &
                                   hy_useMinSmoothVarVal, hy_useMaxSmoothVarVal, &
                                   hy_minSmoothVarVal, hy_maxSmoothVarVal
  use Hydro_data, ONLY: hy_useHydro, hy_unsplitEosMode, hy_eosModeAfter
  use Hydro_data, ONLY: hy_3TMode
  use Hydro_data, ONLY: hy_lam3ScaleFactor
  use Hydro_data, ONLY: hy_computeRadFlBefore, hy_doUnsplitLoop0
       
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Eos_interface, ONLY               : Eos_getParameters
  use Grid_interface, ONLY: Grid_setFluxHandling

  implicit none


#include "constants.h"
#include "Flash.h"  
!!$#include "Hydro_components.h"

  character(len=MAX_STRING_LENGTH) :: str
  
  !!**Hydro_sweep RuntimeParameters

  call RuntimeParameters_get ("eosMode", str)
  call RuntimeParameters_mapStrToInt(str, hy_unsplitEosMode)
  if(hy_unsplitEosMode/=MODE_DENS_EI .AND. hy_unsplitEosMode/=MODE_DENS_EI_ALL .AND. &
       hy_unsplitEosMode/=MODE_DENS_EI_SCATTER .AND. hy_unsplitEosMode/=MODE_DENS_EI_GATHER &
       .AND. hy_unsplitEosMode/=MODE_DENS_EI_RECAL_GATHER &
       .AND. hy_unsplitEosMode/=MODE_DENS_EI_MAT_GATHER_PRADSCALE)&
       call Driver_abortFlash("Hydro : Wrong Eos mode for multiTemp unsplit Hydro")

  call RuntimeParameters_get("hy_lam3ScaleFactor",  hy_lam3ScaleFactor)

#ifdef FLLM_VAR
  hy_computeRadFlBefore = .TRUE.
  if (.NOT. hy_doUnsplitLoop0) then
     hy_doUnsplitLoop0 = hy_computeRadFlBefore
  end if
#endif


  !!**PPM inputs
!!$  call RuntimeParameters_get("ppmEnerFluxConstructionMeth", hy_ppmEnerFluxConstructionMeth)
!!$  call RuntimeParameters_get("ppmEintFluxConstructionMeth", hy_ppmEintFluxConstructionMeth)
!!$  if (hy_ppmEintFluxConstructionMeth==-1) hy_ppmEintFluxConstructionMeth = hy_ppmEnerFluxConstructionMeth
!!$  call RuntimeParameters_get("ppmEnerCompFluxConstructionMeth", hy_ppmEnerCFluxConstructionMeth)
!!$  call RuntimeParameters_get("ppmEintCompFluxConstructionMeth", hy_ppmEintCFluxConstructionMeth)
!!$  if (hy_ppmEintCFluxConstructionMeth==-1) hy_ppmEintCFluxConstructionMeth = hy_ppmEnerCFluxConstructionMeth

!!$  call RuntimeParameters_get("charLimiting", hy_charLimiting) ! new characteristic limiting - DL
!!$
!!$  call PhysicalConstants_get("electron mass",hy_eMass)
!!$  call PhysicalConstants_get("proton mass",hy_pMass)
!!$  call PhysicalConstants_get("electron mass",hy_eMassInUAmu,unitMass="amu")


  if (.NOT. hy_useHydro) return ! If Hydro is turned off; return here before anything serious gets done.

  

  ! For testing ways to advect components and handle shock heating
  
  call RuntimeParameters_get ("hy_eosModeAfter", str)
  call RuntimeParameters_mapStrToInt(str, hy_eosModeAfter)
  if(hy_eosModeAfter/=MODE_DENS_EI_SELE_GATHER .AND. &
       hy_eosModeAfter/=MODE_DENS_EI_SCATTER .AND. &
       hy_eosModeAfter/=MODE_DENS_EI_GATHER .AND. &
       hy_eosModeAfter/=MODE_DENS_EI_RECAL_GATHER .AND. &
       hy_eosModeAfter/=hy_unsplitEosMode)&
       call Driver_abortFlash("Hydro : Wrong Eos mode for hy_eosModeAfter")
  call Eos_getParameters(smallE1=hy_smallEion, smallE2=hy_smallEele, smallE3=hy_smallErad)
  call RuntimeParameters_get ("hy_3Ttry_B", hy_3Ttry_B)
  call RuntimeParameters_get ("hy_3Ttry_D", hy_3Ttry_D)
  call RuntimeParameters_get ("hy_3Ttry_E", hy_3Ttry_E)
  call RuntimeParameters_get ("hy_3Ttry_F", hy_3Ttry_F)
  call RuntimeParameters_get ("hy_3Ttry_G", hy_3Ttry_G)

  call RuntimeParameters_get ("hy_3Ttry_B_rad", hy_3Ttry_B_rad)
  if (hy_3Ttry_B_rad < 0) hy_3Ttry_B_rad = hy_3Ttry_B

  hy_mtFactorB0 = 1.0; hy_mtFactorBrad0 = 0.0
  hy_mtFactorB1 = 0.0; hy_mtFactorBrad1 = 0.0
  hy_mtFactorB2 = 0.0; hy_mtFactorBrad2 = 1.0

  if (hy_3Ttry_B==0) then
     hy_mtFactorB0 = 1.0
     hy_mtFactorB1 = 0.0
     hy_mtFactorB2 = 0.0
  else if (hy_3Ttry_B==1) then
     hy_mtFactorB0 = 0.0
     hy_mtFactorB1 = 1.0
     hy_mtFactorB2 = 0.0
  else if (hy_3Ttry_B==2) then
     hy_mtFactorB0 = 0.0
     hy_mtFactorB1 = 0.0
     hy_mtFactorB2 = 1.0
  end if
  if (hy_3Ttry_B_rad==0) then
     hy_mtFactorBrad0 = 1.0
     hy_mtFactorBrad1 = 0.0
     hy_mtFactorBrad2 = 0.0
  else if (hy_3Ttry_B_rad==1) then
     hy_mtFactorBrad0 = 0.0
     hy_mtFactorBrad1 = 1.0
     hy_mtFactorBrad2 = 0.0
  else if (hy_3Ttry_B_rad==2) then
     hy_mtFactorBrad0 = 0.0
     hy_mtFactorBrad1 = 0.0
     hy_mtFactorBrad2 = 1.0
  end if


#ifndef SELE_MSCALAR
  if(hy_eosModeAfter == MODE_DENS_EI_SELE_GATHER) then
     call Logfile_stampMessage('[hy_uhd_MultiTempInit] ERROR:')
     call Logfile_stampMessage('[hy_uhd_MultiTempInit] You set hy_eosModeAfter to "dens_ie_sele_gather"')
     call Logfile_stampMessage('[hy_uhd_MultiTempInit] To use this mode, the electron entropy mass scalar')
     call Logfile_stampMessage('[hy_uhd_MultiTempInit] must exist. Make sure the following line exists in')
     call Logfile_stampMessage('[hy_uhd_MultiTempInit] your simulation Config file:')
     call Logfile_stampMessage('[hy_uhd_MultiTempInit] MASS_SCALAR sele EOSMAP: SELE')
     call Driver_abortFlash("[hy_uhd_MultiTempInit] Must create SELE_MSCALAR, see log file for details")
  end if
#endif

  call RuntimeParameters_get("hy_3TMode", str)
  if(trim(str) == "ragelike") then
     hy_3TMode = HY3T_RAGELIKE
!     if(hy_eosModeAfter /= MODE_DENS_EI_GATHER) &
!          call Driver_abortFlash( &
!          "[hy_uhd_MultiTempInit] ERROR: hy_eosModeAfter must be dens_ie_gather")

  else if(trim(str) == "crashlike") then
     hy_3TMode = HY3T_CRASHLIKE
     if(hy_eosModeAfter /= MODE_DENS_EI_GATHER) &
          call Driver_abortFlash( &
          "[hy_uhd_MultiTempInit] ERROR: hy_eosModeAfter must be dens_ie_gather")

  else if(trim(str) == "entropy") then
     hy_3TMode = HY3T_ENTROPY
     if(hy_eosModeAfter /= MODE_DENS_EI_SELE_GATHER) &
          call Driver_abortFlash( &
          "[hy_uhd_MultiTempInit] ERROR: hy_eosModeAfter must be dens_ie_sele_gather")

  else if(trim(str) == "castrolike") then
     hy_3TMode = HY3T_CASTROLIKE
     call Logfile_stampMessage('[hy_uhd_MultiTempInit] Using the experimental CASTROLIKE 3T mode!!!')

  else
     call Driver_abortFlash("[hy_uhd_MultiTempInit] ERROR: Unknown hy_3TMode")
  end if

  call RuntimeParameters_get("hy_mtScaleWork", hy_mtScaleWork)
  call RuntimeParameters_get("hy_mtScaleAccel", hy_mtScaleAccel)
  call RuntimeParameters_get("hy_mtScaleLorentz", hy_mtScaleLorentz)
  call RuntimeParameters_get("hy_mtPresRatLambda3Min", hy_mtPresRatLambda3Min)

  call RuntimeParameters_get("hy_smoothIterations", hy_smoothIterations)
  call RuntimeParameters_get("hy_smoothMethod", str)
  call RuntimeParameters_mapStrToInt(str, hy_smoothMethod)
  call RuntimeParameters_get("hy_smoothCoeff", hy_smoothCoeff)
  call RuntimeParameters_get("hy_useMinSmoothVarVal", hy_useMinSmoothVarVal)
  call RuntimeParameters_get("hy_useMaxSmoothVarVal", hy_useMaxSmoothVarVal)
  call RuntimeParameters_get("hy_minSmoothVarVal", hy_minSmoothVarVal)
  call RuntimeParameters_get("hy_maxSmoothVarVal", hy_maxSmoothVarVal)



end subroutine hy_uhd_MultiTempInit
