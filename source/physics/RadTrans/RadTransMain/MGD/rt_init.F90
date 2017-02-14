!!****if* source/physics/RadTrans/RadTransMain/MGD/rt_init
!!
!!  NAME 
!!
!!  rt_init
!!
!!  SYNOPSIS
!!
!!  call rt_init()
!!
!!  DESCRIPTION 
!!    Initialize local data for each radiative transfer model
!!
!!***
subroutine rt_init
  use rt_data
  use RadTrans_data, ONLY: rt_useRadTrans, rt_meshCopyCount, &
       rt_useRadTrans,                                       &
       rt_meshMe, rt_acrossMe, rt_radconst, rt_boltz
  use RadTrans_interface, ONLY: RadTrans_mgdSetBound
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
       RuntimeParameters_mapStrToInt
  use Logfile_interface, ONLY: Logfile_stampMessage, Logfile_stampVarMask
  implicit none

#include "Flash.h"
#include "constants.h"

  character(len=MAX_STRING_LENGTH) :: xl_bcString,xr_bcString
  character(len=MAX_STRING_LENGTH) :: yl_bcString,yr_bcString
  character(len=MAX_STRING_LENGTH) :: zl_bcString,zr_bcString
  character(len=MAX_STRING_LENGTH) :: flmode_str
  character(len=MAX_STRING_LENGTH) :: entry_str, basename, fullname
  integer :: entry_type
  integer :: g, f
  real, parameter :: EV_TO_ERGS = 1.60217653e-12
  real :: group_bound
  real :: bcValues(6)
  real :: eg, egp1, xg, xgp1, KB, pxg, pxgp1, A
  
  ! Load flux limiter options:
  call RuntimeParameters_get('rt_mgdFlMode', flmode_str)
  call makeLowercase(flmode_str)
  call RuntimeParameters_mapStrToInt(flmode_str, rt_mgdFlMode)
  call RuntimeParameters_get('rt_mgdFlCoef', rt_mgdFlCoef)

  ! Multigroup Diffusion Number of Groups:
  call RuntimeParameters_get('rt_mgdNumGroups', rt_mgdNumGroups)

  call RuntimeParameters_get('rt_useMGD', rt_useMGD)
  if (rt_useMGD .AND. .NOT. rt_useRadTrans) then
     if (rt_meshME==MASTER_PE .AND. &
          rt_acrossMe==MASTER_PE) then
        print *,' Forcing rt_useMGD to FALSE because useRadTrans is FALSE!'
     end if
     rt_useMGD = .FALSE.
  end if
  if(.not. rt_useRadTrans) return

  ! Do some error checking to detect configuration problems:
#ifndef MGDR_NONREP
  call Driver_abortFlash("MGD Error: MGDR_NONREP not defined")
#endif

  ! Get options relating to timing:
  call RuntimeParameters_get('rt_timeGroups', rt_timeGroups)
  call RuntimeParameters_get('rt_groupBarrier', rt_groupBarrier)

  if(rt_timeGroups .and. rt_meshCopyCount > 1) then
     call Driver_abortFlash("MGD Error: rt_timeGroups cannot be true when meshCopyCount > 1")
  end if

  ! Read flag which will determine whether a RadTrans time step
  ! associated with MGD is computed:
  call RuntimeParameters_get('rt_computeDt', rt_computeDt)
  
  ! Tolerance for relative temp change in a time step, used by SOME
  ! RadTrans implementations(s)
  call RuntimeParameters_get('rt_tempChangeRelTol', rt_tempChangeRelTol)
  
  if(rt_mgdNumGroups <= 0) then
     call Driver_abortFlash("MGD Error: rt_mgdNumGroups must be greater than zero to use MGD!")
  end if

  if(rt_mgdNumGroups > MGDR_NONREP_MAXLOCS*rt_meshCopyCount) then
     ! If you get this error, you have to increase the amount of
     ! memory available for storing radiation groups. You do this by
     ! specifying a larger number for the mgd_meshgroups setup option.
     ! You must meet the condition:
     !   rt_mgdNumGroups > mgd_meshgroups*rt_meshCopyCount
     call Driver_abortFlash("MGD Error: Not enough unk space to store group energies!")
  end if

  ! Allocate space for energy group boundaries:
  allocate(rt_mgdBounds(rt_mgdNumGroups + 1))

  ! Set MGD Group Boundaries:
  call RuntimeParameters_get('rt_mgdBoundEntry', entry_str)
  call RuntimeParameters_mapStrToInt(entry_str, entry_type)

  ! Load the group boundaries
  select case(entry_type)
  case (GRBD_MANUAL)
     if(rt_mgdNumGroups > MGD_MAXGROUPS) then
        ! There are only enough runtime parameters defined to store
        ! the group boundaries for MGD_MAXGROUPS groups. If this isn't
        ! enough, just set mgd_maxgroups to a larger number on the
        ! setup line.
        call Driver_abortFlash("MGD Error: MGD_MAXGROUPS is too small")
     end if

     ! Loop through each group boundary and store it. Note that the
     ! group boundaries are input in units of eV. These have to be
     ! converted to ergs.
     basename = 'rt_mgdBounds_'
     do g = 1, rt_mgdNumGroups+1
        call concatStringWithInt(basename, g, fullname)
        call RuntimeParameters_get(fullname, group_bound)

        if (group_bound < 0.0) then
           call Logfile_stampMessage( &
                "[rt_init] Error: One of the group boundaries was negative.")
           call Logfile_stampMessage( &
                "This probably means that a group boundary was not specified,")
           call Logfile_stampMessage( &
                "check your paramter file to make sure you specified boundaries! ")
           call Driver_abortFlash("MGD Error: Bad group boundary, see logfile")
        end if

        call RadTrans_mgdSetBound(g, group_bound * EV_TO_ERGS)
     end do
     
  case DEFAULT
     call Driver_abortFlash("MGD Error: Bad value for rt_mgdBoundEntry")
  end select

  if(.not. rt_useMGD) return

  ! Load MGD boundary conditions:
  allocate(rt_mgdDomainBC(rt_mgdNumGroups, 6))
  allocate(rt_mgdBcVals(rt_mgdNumGroups, 6))

  call RuntimeParameters_get("rt_mgdXlBoundaryType", xl_bcString)
  call RuntimeParameters_get("rt_mgdXrBoundaryType", xr_bcString)
  call RuntimeParameters_get("rt_mgdYlBoundaryType", yl_bcString)
  call RuntimeParameters_get("rt_mgdYrBoundaryType", yr_bcString)
  call RuntimeParameters_get("rt_mgdZlBoundaryType", zl_bcString)
  call RuntimeParameters_get("rt_mgdZrBoundaryType", zr_bcString)

  call RuntimeParameters_get("rt_mgdXlBoundaryTemp", bcValues(1))
  call RuntimeParameters_get("rt_mgdXrBoundaryTemp", bcValues(2))
  call RuntimeParameters_get("rt_mgdYlBoundaryTemp", bcValues(3))
  call RuntimeParameters_get("rt_mgdYrBoundaryTemp", bcValues(4))
  call RuntimeParameters_get("rt_mgdZlBoundaryTemp", bcValues(5))
  call RuntimeParameters_get("rt_mgdZrBoundaryTemp", bcValues(6))

  KB = rt_boltz
  A  = rt_radconst

  do g = 1, rt_mgdNumGroups
     call RuntimeParameters_mapStrToInt(xl_bcString,rt_mgdDomainBC(g,1))
     call RuntimeParameters_mapStrToInt(xr_bcString,rt_mgdDomainBC(g,2))
     call RuntimeParameters_mapStrToInt(yl_bcString,rt_mgdDomainBC(g,3))
     call RuntimeParameters_mapStrToInt(yr_bcString,rt_mgdDomainBC(g,4))
     call RuntimeParameters_mapStrToInt(zl_bcString,rt_mgdDomainBC(g,5))
     call RuntimeParameters_mapStrToInt(zr_bcString,rt_mgdDomainBC(g,6))

     do f = 1, 6
        if(bcValues(f) > 0.0) then
           call RadTrans_mgdGetBound(g, eg)
           call RadTrans_mgdGetBound(g+1, egp1)
           xg   = eg / (KB*bcValues(f))
           xgp1 = egp1 / (KB*bcValues(f))
           call RadTrans_planckInt(xg, pxg)
           call RadTrans_planckInt(xgp1, pxgp1)
           pxg = pxg * 15.0/PI**4
           pxgp1 = pxgp1 * 15.0/PI**4
           rt_mgdBcVals(g,f) = A * bcValues(f)**4 * (pxgp1 - pxg)
        else
           rt_mgdBcVals(g,f) = 0.0
        end if
     end do
  end do  
  
  call RuntimeParameters_get("rt_mgdthetaImplct", rt_mgdthetaImplct)

  rt_gcMask = .TRUE.

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


  call Logfile_stampVarMask(rt_gcMask, .FALSE., '[rt_init]', 'gcNeed')

  call rt_initUnifiedMGD
  call rt_initExpRelax

end subroutine rt_init
