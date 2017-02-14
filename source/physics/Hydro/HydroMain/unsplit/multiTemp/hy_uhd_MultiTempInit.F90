!!****if* source/physics/Hydro/HydroMain/unsplit/multiTemp/hy_uhd_MultiTempInit
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
  use Hydro_data, ONLY: hy_useHydro, hy_unsplitEosMode, hy_eosModeAfter
  use Hydro_data, ONLY: hy_3TMode
       
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
       .AND. hy_unsplitEosMode/=MODE_DENS_EI_RECAL_GATHER)&
       call Driver_abortFlash("Hydro : Wrong Eos mode for multiTemp unsplit Hydro")

 

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
     if(hy_eosModeAfter /= MODE_DENS_EI_GATHER) &
          call Driver_abortFlash( &
          "[hy_uhd_MultiTempInit] ERROR: hy_eosModeAfter must be dens_ie_gather")

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
  else
     call Driver_abortFlash("[hy_uhd_MultiTempInit] ERROR: Unknown hy_3TMode")
  end if

end subroutine hy_uhd_MultiTempInit
