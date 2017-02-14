!!****if* source/physics/RadTrans/RadTransMain/MGD/Unified/rt_initUnifiedMGD
!!
!!  NAME 
!!
!!  rt_initUnifiedMGD
!!
!!  SYNOPSIS
!!
!!  call rt_initUnifiedMGD()
!!
!!  DESCRIPTION 
!!    Initialize additional data for using the Unified MGD solver for radiative transfer
!!
!!***
subroutine rt_initUnifiedMGD
  use rt_data, ONLY: rt_useMGD, rt_mgdthetaC, rt_mgdthetaD, rt_mgdthetaImplct, &
       rt_tightIonCoupling, &
       rt_gcMask
  use RadTrans_data, ONLY: rt_useRadTrans, rt_meshCopyCount, &
       rt_useRadTrans,                                       &
       rt_meshMe, rt_acrossMe
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY: Logfile_stampMessage, Logfile_stampVarMask
  implicit none

#include "Flash.h"
#include "constants.h"

  
  call RuntimeParameters_get('rt_useMGD', rt_useMGD)
  if (rt_useMGD .AND. .NOT. rt_useRadTrans) then
     if (rt_meshME==MASTER_PE .AND. &
          rt_acrossMe==MASTER_PE) then
        print *,' Forcing rt_useMGD to FALSE because useRadTrans is FALSE!'
     end if
     rt_useMGD = .FALSE.
  end if
  if(.not. rt_useRadTrans) return

  if(.not. rt_useMGD) return

  call RuntimeParameters_get("rt_mgdthetaC", rt_mgdthetaC)
  call RuntimeParameters_get("rt_mgdthetaD", rt_mgdthetaD)
  if (rt_mgdthetaC == -1.0) rt_mgdthetaC = rt_mgdthetaImplct
  if (rt_mgdthetaD == -1.0) rt_mgdthetaD = rt_mgdthetaImplct

  call RuntimeParameters_get("rt_tightIonCoupling", rt_tightIonCoupling)

!!$  rt_gcMask = .TRUE.

#if NSPECIES == 1
#ifdef SPECIES_BEGIN
  rt_gcMask(SPECIES_BEGIN) = .FALSE.
#endif
#endif
#ifdef MGDC_VAR
  rt_gcMask(MGDC_VAR) = .FALSE.
#endif
#ifdef COND_VAR
  rt_gcMask(COND_VAR) = .FALSE.
#endif
#ifdef DFCF_VAR
  rt_gcMask(DFCF_VAR) = .FALSE.
#endif
#ifdef ABSR_VAR
  rt_gcMask(ABSR_VAR) = .FALSE.
#endif
#ifdef EMIS_VAR
  rt_gcMask(EMIS_VAR) = .FALSE.
#endif
#ifdef FLLM_VAR
  rt_gcMask(FLLM_VAR) = .FALSE.
#endif
#ifdef GAME_VAR
  rt_gcMask(GAME_VAR) = .FALSE.
#endif
#ifdef PRES_VAR
  rt_gcMask(PRES_VAR) = .FALSE.
#endif
#ifdef PELE_VAR
  rt_gcMask(PELE_VAR) = .FALSE.
#endif
#ifdef PION_VAR
  rt_gcMask(PION_VAR) = .FALSE.
#endif
#ifdef PIPE_VAR
  rt_gcMask(PIPE_VAR) = .FALSE.
#endif
#ifdef TITE_VAR
  rt_gcMask(TITE_VAR) = .FALSE.
#endif
#ifdef SHKS_VAR
  rt_gcMask(SHKS_VAR) = .FALSE.
#endif
#ifdef DBGS_VAR
  rt_gcMask(DBGS_VAR) = .FALSE.
#endif
#ifdef SHOK_VAR
  rt_gcMask(SHOK_VAR) = .FALSE.
#endif


  call Logfile_stampVarMask(rt_gcMask, .FALSE., '[rt_initUnifiedMGD]', 'gcNeed')

end subroutine rt_initUnifiedMGD
