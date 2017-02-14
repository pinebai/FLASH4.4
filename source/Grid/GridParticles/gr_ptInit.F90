!!****if* source/Grid/GridParticles/gr_ptInit
!!
!! NAME
!!
!!  gr_ptInit
!!
!! SYNOPSIS
!!
!!  gr_ptInit()
!!
!! DESCRIPTION
!!
!!  Initialize values for all data in the module gr_ptData,
!!  and allocate the scratch buffers
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!         !!! IMPORTANT !!!
!!
!!  The allocations that are done here right now are EXTREMELY dangerous!
!!  They fix the dimensions of the grid particle destination and source
!!  buffers assuming a FIXED 1st dimension size. This has the potential
!!  to give catastrophic results when using the grid particles unit for
!!  many different kind of particles with different particle attribute
!!  sizes. This works well for those particles whose attribute sizes
!!  match the fixed 1st dimension of the buffers, but fails to work for
!!  those whose attribute sizes do not match the fixed 1st dimension of
!!  the buffers.
!!
!!
!!***

#include "Flash.h"
#ifdef FLASH_PIMG
#include "ProtonImaging.h"
#endif
#ifdef FLASH_PEMI
#include "ProtonEmission.h"
#endif
#ifdef FLASH_EDEP
#include "EnergyDeposition.h"
#endif

subroutine gr_ptInit()
  use gr_ptData, ONLY :   gr_ptDestBuf,gr_ptSourceBuf,gr_ptMaxPerProc,&
                          gr_ptRemove,gr_ptRemoveAlgo,&
                          gr_ptNumToReduce,gr_ptSieveCheckFreq,&
                          gr_ptLogLevel,gr_ptKeepLostParticles,gr_ptMaxVirtualCount
#ifndef FLASH_GRID_UG
  use gr_ptData, ONLY : gr_ptMaxPerProcUpperThresh, gr_ptMaxPerProcLowerThresh, &
       gr_ptMaxPerProcBlockFactor,gr_ptMaxPerProcBlockNoFuzz, gr_ptRefineOnPtMaxPerProc
  use Grid_data, ONLY : gr_meshMe, gr_refineOnParticleCount
#endif  
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stampMessage, Logfile_stamp
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none 
  integer, save :: maxPerProc, propCount
 

  call RuntimeParameters_get("gr_ptRemove",gr_ptRemove)
  call RuntimeParameters_get("gr_ptRemoveAlgo",gr_ptRemoveAlgo)
  call RuntimeParameters_get("gr_ptNumToReduce",gr_ptNumToReduce)
  call RuntimeParameters_get("gr_ptSieveCheckFreq",gr_ptSieveCheckFreq)
  call RuntimeParameters_get("keepLostParticles",gr_ptKeepLostParticles)

!  call RuntimeParameters_get("pt_logLevel",gr_ptLogLevel)


  gr_ptMaxPerProc=1
  propCount=1

#ifdef NPART_PROPS
  if (NPART_PROPS > 1) then
     call RuntimeParameters_get("pt_maxPerProc",maxPerProc)
     propCount=max(NPART_PROPS,propCount)                       ! DANGER
     gr_ptMaxPerProc = max(maxPerProc,gr_ptMaxPerProc)          ! OK -> 2nd buffer dimension
  endif
#endif
  gr_ptMaxVirtualCount=gr_ptMaxPerProc

#ifdef RAY_ATTR_COUNT
  call RuntimeParameters_get("ed_maxRayCount",maxPerProc)
  propCount=max(RAY_ATTR_COUNT,propCount)                       ! DANGER
  gr_ptMaxPerProc = max(maxPerProc,gr_ptMaxPerProc)             ! OK -> 2nd buffer dimension
#endif

#ifdef PROTON_ATTRCOUNT
  call RuntimeParameters_get("pi_maxProtonCount",maxPerProc)
  propCount=max(PROTON_ATTRCOUNT,propCount)                     ! DANGER
  gr_ptMaxPerProc = max(maxPerProc,gr_ptMaxPerProc)             ! OK -> 2nd buffer dimension
#endif

#ifdef EMPROTON_ATTRCOUNT
  call RuntimeParameters_get("pem_maxProtonCount",maxPerProc)
  propCount=max(EMPROTON_ATTRCOUNT,propCount)                   ! DANGER
  gr_ptMaxPerProc = max(maxPerProc,gr_ptMaxPerProc)             ! OK -> 2nd buffer dimension
#endif
!
!
!    ...The following check has been added, such that the FLASH application
!       stops smoothly with an informative error massage, rather than an
!       uncontrolled code crash.
!
!
#ifdef NPART_PROPS
  if (NPART_PROPS > 1 .and. NPART_PROPS /= propCount) then
      call Logfile_stamp (NPART_PROPS,"[gr_ptInit]: Value of NPART_PROPS is = ")
      call Logfile_stamp (propCount,"[gr_ptInit]: Value of propCount is = ")
      call Logfile_stampMessage ("Different 1st dimensions for Destination and Source Buffer!")
      call Driver_abortFlash ("[gr_ptInit]: NPART_PROPS must match propCount (see Logfile).")
  endif
#endif

#ifdef RAY_ATTR_COUNT
  if (RAY_ATTR_COUNT /= propCount) then
      call Logfile_stamp (RAY_ATTR_COUNT,"[gr_ptInit]: Value of RAY_ATTR_COUNT is = ")
      call Logfile_stamp (propCount,"[gr_ptInit]: Value of propCount is = ")
      call Logfile_stampMessage ("Different 1st dimensions for Destination and Source Buffer!")
      call Driver_abortFlash ("[gr_ptInit]: RAY_ATTR_COUNT must match propCount (see Logfile).")
  endif
#endif

#ifdef PROTON_ATTRCOUNT
  if (PROTON_ATTRCOUNT /= propCount) then
      call Logfile_stamp (PROTON_ATTRCOUNT,"[gr_ptInit]: Value of PROTON_ATTRCOUNT is = ")
      call Logfile_stamp (propCount,"[gr_ptInit]: Value of propCount is = ")
      call Logfile_stampMessage ("Different 1st dimensions for Destination and Source Buffer!")
      call Driver_abortFlash ("[gr_ptInit]: PROTON_ATTRCOUNT must match propCount (see Logfile).")
  endif
#endif


#ifdef EMPROTON_ATTRCOUNT
  if (EMPROTON_ATTRCOUNT /= propCount) then
      call Logfile_stamp (EMPROTON_ATTRCOUNT,"[gr_ptInit]: Value of EMPROTON_ATTRCOUNT is = ")
      call Logfile_stamp (propCount,"[gr_ptInit]: Value of propCount is = ")
      call Logfile_stampMessage ("Different 1st dimensions for Destination and Source Buffer!")
      call Driver_abortFlash ("[gr_ptInit]: EMPROTON_ATTRCOUNT must match propCount (see Logfile).")
  endif
#endif


  ! These are deallocated in gr_ptFinalize
  allocate(gr_ptDestBuf(propCount,gr_ptMaxPerProc))             ! DANGER -> propCount fixed
  allocate(gr_ptSourceBuf(propCount,gr_ptMaxPerProc))           ! DANGER -> propCount fixed

#ifndef FLASH_GRID_UG
  ! additional RPs for maxPerProc-based refinement criteria - KW
  call RuntimeParameters_get("gr_ptMaxPerProcUpperThresh",gr_ptMaxPerProcUpperThresh)
  call RuntimeParameters_get("gr_ptMaxPerProcLowerThresh",gr_ptMaxPerProcLowerThresh)
  call RuntimeParameters_get("gr_ptMaxPerProcBlockFactor",gr_ptMaxPerProcBlockFactor)
  call RuntimeParameters_get("gr_ptMaxPerProcBlockNoFuzz",gr_ptMaxPerProcBlockNoFuzz)
  call RuntimeParameters_get("gr_ptRefineOnPtMaxPerProc", gr_ptRefineOnPtMaxPerProc)
  if (gr_ptRefineOnPtMaxPerProc .AND. .NOT. gr_refineOnParticleCount) then
     call Logfile_stampMessage( &
          "WARNING: Ignoring gr_ptRefineOnPtMaxPerProc because RP refine_on_particle_count is .FALSE.")
     gr_ptRefineOnPtMaxPerProc = .FALSE.
  end if
#endif
  return
end subroutine gr_ptInit
