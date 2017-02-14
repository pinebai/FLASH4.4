!!****if* source/physics/RadTrans/RadTransMain/NeutrinoLeakage/rt_data
!!
!!  NAME 
!!    rt_data
!!
!!  SYNOPSIS
!!    use rt_data
!!
!!  DESCRIPTION 
!!    Stores data for Leakage
!!
!!  NOTES
!!      This unit implements ray-by-ray multispecies neutrino leakage.
!!      Parts of this unit are released under a different license than the
!!      usual FLASH license.  Specifically, some subroutines in rt_calcLeak.F90 and 
!!      rt_calcTau.F90 are released under the Creative Commons 
!!      attribution-noncommercial-share alike license.  Basically, if you use this
!!      unit in your work, the license requires that you cite the two articles 
!!      mentioned below.  More details may be found here:  stellarcollapse.org/codes.html.
!!
!!      * O'Connor, E.P., & Ott, C.D. 2010, CQGra, 27, 114103
!!      * Couch, S.M., & O'Connor, E.P. 2013, arXiv:1310.5728
!!
!!***
module rt_data
#include "constants.h"
#include "Flash.h"

  implicit none

  real, save :: rt_leakRadLog, rt_leakRadMax, rt_leakThtMax, rt_leakDx, rt_leakPhiMax
  real, pointer, save :: rt_tauRuff(:,:,:), rt_chiRoss(:,:,:), rt_heatFlx(:,:,:)
  real, pointer, save :: rt_heatErms(:,:), rt_heatEave(:,:)
#ifdef LEAK_STATIC
  integer, parameter :: rt_leakNumRad = 1000, rt_leakNumTht = 37, rt_leakNumPhi = 75, rt_leakNumRays = 2775
  real, save :: rt_leakRadii(rt_leakNumRad), rt_leakTheta(rt_leakNumTht), rt_leakPhi(rt_leakNumPhi)
  real, save :: rt_leakX(rt_leakNumRad,rt_leakNumTht,rt_leakNumPhi), &
       rt_leakY(rt_leakNumRad,rt_leakNumTht,rt_leakNumPhi), &
       rt_leakZ(rt_leakNumRad,rt_leakNumTht,rt_leakNumPhi)
  real, save :: rt_leakArr(3,rt_leakNumRad,rt_leakNumRays), rt_leakSrc(4,rt_leakNumRad,rt_leakNumRays)
  real, save :: rt_dr(rt_leakNumRad) 
  real, target, save :: rt_rayData(5,rt_leakNumRad,3,rt_leakNumRays)
#else
  integer, save :: rt_leakNumRad, rt_leakNumTht, rt_leakNumPhi, rt_leakNumRays
  real, allocatable, save :: rt_leakRadii(:), rt_leakTheta(:), rt_leakPhi(:)
  real, allocatable, save :: rt_leakX(:,:,:), rt_leakY(:,:,:), rt_leakZ(:,:,:)
  real, allocatable, save :: rt_leakArr(:,:,:), rt_leakSrc(:,:,:)
  real, allocatable, save :: rt_dr(:) 
  real, allocatable, target, save :: rt_rayData(:,:,:,:)
#endif
  integer, allocatable :: rt_recvCnt(:), rt_dsplCnt(:)

  integer, save :: rt_globalMe, rt_meshComm, rt_meshNumProcs, rt_subMeshMe, rt_subCommSize, rt_subMeshComm
  logical, save :: rt_leakDoHeat
  real, save :: rt_leakHeatFac
  integer, save :: rt_leakRay
  real, save :: rt_dTht, rt_dPhi
  integer, save :: rt_arraySize, rt_arraySmall, rt_arraySizeLocal
  
  real, save :: rt_leakLumTot(3), rt_leakHeatTot(3), rt_leakNetTot(3), rt_leakEaveTot(3)

  real,parameter :: temp_mev_to_kelvin = 1.16045221d10

  integer, save :: rt_minTht, rt_maxTht

  real, save :: rt_pnsCoord(MDIM), rt_pnsDens

  logical, save :: rt_threadWithinBlock, rt_threadBlockList
  integer, save :: rt_istart, rt_iend
  integer, parameter :: rt_itau=1, rt_ichi=2, rt_iflx=3, rt_ierm=4, rt_ieav=5
  real, save :: rt_radNu(3)
  integer, save :: rt_nstepStart, rt_nstep
  real, save :: rt_bounceTime
  real, save :: rt_reducedTime
  integer, save :: rt_reducedSteps

end module rt_data
