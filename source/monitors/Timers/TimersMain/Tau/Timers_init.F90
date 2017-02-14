!!****if* source/monitors/Timers/TimersMain/Tau/Timers_init
!!
!! NAME
!!  Timers_init
!!
!! SYNOPSIS
!!
!!  Timers_init(real(OUT) :: initialWCTime)
!!  
!! DESCRIPTION 
!!  
!!  Initialize the timer data structures.  This will
!!  essentially delete all information previously gathered by all timers
!!  and make it safe to start timers from scratch.  In the middle of
!!  a run, for instance, this could be called once per timestep along with
!!  Timers_getSummary to get timer summary information for each timestep.
!!  
!! ARGUMENTS 
!!
!!  initialWCTime -- the initial wall clock time when this was called. 
!! 
!!***

subroutine Timers_init(initialWCTime)

  use TauMetadata_interface, ONLY: &
    write_tau_metadata_int, write_tau_metadata_real, &
    write_tau_metadata_str, write_tau_metadata_log
  use Timers_data, ONLY : tmr_freeSlot, tmr_globalMe, tmr_globalNumProcs, &
       tmr_MAX_CUSTOM_TIMERS, tmr_customPrefix, tmr_prefixLen
  use Driver_interface, ONLY : Driver_abortFlash
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, Driver_getComm
  use RuntimeParameters_interface, ONLY : &
    RuntimeParameters_getNumReal, RuntimeParameters_getNumInt, &
    RuntimeParameters_getNumIgn, RuntimeParameters_getNumStr, &
    RuntimeParameters_getNumLog, RuntimeParameters_getAllInt, &
    RuntimeParameters_getAllReal, RuntimeParameters_getAllStr, &
    RuntimeParameters_getAllLog
  use tmr_interface, ONLY : tmr_etime

  implicit none

#include "constants.h"
#include "Flash.h"

  real, intent(out) :: initialWCTime

  ! Everybody should know this

  character(len=MAX_STRING_LENGTH)   :: names, names2
  integer                            :: i

  character(len=80)                  :: setup_flashRelease
  character (len=MAX_STRING_LENGTH)  :: buildDate, buildDir, buildMachine
  character (len=400)                :: cflags, fflags, setupCall

  integer :: numRealParms, numIntParms, numStrParms, numLogParms

  character (len=MAX_STRING_LENGTH), allocatable, save :: intParmNames(:)
  integer, allocatable, save :: intParmValues(:)
  logical, allocatable, save :: intParmChanged(:)

  character (len=MAX_STRING_LENGTH), allocatable, save :: realParmNames(:)
  real, allocatable, save :: realParmValues(:)
  logical, allocatable, save :: realParmChanged(:)

  character (len=MAX_STRING_LENGTH), allocatable, save :: strParmNames(:)
  character (len=MAX_STRING_LENGTH), allocatable, save :: strParmValues(:)
  logical, allocatable, save :: strParmChanged(:)

  character (len=MAX_STRING_LENGTH), allocatable, save :: logParmNames(:)
  logical, allocatable, save :: logParmValues(:)
  logical, allocatable, save :: logParmChanged(:)

  character(len=MAX_STRING_LENGTH), allocatable, save :: s1(:)

  integer :: numUnits

  character(len=7) :: unitstr = "unit   "

!==============================================================================

  call Driver_getMype(GLOBAL_COMM, tmr_globalMe)
  call Driver_getNumProcs(GLOBAL_COMM, tmr_globalNumProcs)


  tmr_freeSlot = 1   !The first place in timer_tauList to store data.
  if (tmr_freeSlot > tmr_MAX_CUSTOM_TIMERS) then
     call Driver_abortFlash("[Timers_init]: No space for timers")
  end if
  tmr_prefixLen = len(tmr_customPrefix)

  call tmr_etime(initialWCTime)

!==============================================================================

! Add metadata regarding this run to the Tau profile.  Include everything that
! is normally written to the log file header.

  call write_tau_metadata_int("Number of processors",     tmr_globalNumProcs)
  call write_tau_metadata_int("Dimensionality",           NDIM)
  call write_tau_metadata_int("Max blocks per processor", MAXBLOCKS)

#ifdef FIXEDBLOCKSIZE
  call write_tau_metadata_int("Number x zones per block", NXB)
  call write_tau_metadata_int("Number y zones per block", NYB)
  call write_tau_metadata_int("Number z zones per block", NZB)
#else
  call write_tau_metadata_str("Number x zones per block", "variable")
  call write_tau_metadata_str("Number y zones per block", "variable")
  call write_tau_metadata_str("Number z zones per block", "variable")
#endif

  call setup_buildstamp (names, names2, MAX_STRING_LENGTH)
  call write_tau_metadata_str("Setup stamp", trim(names))
  call write_tau_metadata_str("Build stamp", trim(names2))

  call setup_systemInfo (names, MAX_STRING_LENGTH)
  call write_tau_metadata_str("System info", trim(names))
  call write_tau_metadata_str("Flash version", setup_flashRelease())

  call setup_buildstats(buildDate, buildDir, buildMachine, setupCall, &
                        cflags, fflags)
  call write_tau_metadata_str("Build directory",        buildDir)
  call write_tau_metadata_str("Setup syntax",           setupCall)
  call write_tau_metadata_str("Fortran compiler flags", fflags)
  call write_tau_metadata_str("C compiler flags",       cflags)

  call setup_getNumFlashUnits(numUnits)
  allocate(s1(numUnits))
  call log_getUnitsArr(numUnits, s1)
  do i = 1, numUnits
    write(unitstr(5:7),'(I3.3)') i
    call write_tau_metadata_str(unitstr, trim(s1(i)))
  enddo
  deallocate(s1)

  call RuntimeParameters_getNumReal(numRealParms)
  call RuntimeParameters_getNumInt(numIntParms)
  call RuntimeParameters_getNumStr(numStrParms)
  call RuntimeParameters_getNumLog(numLogParms)

  allocate (realParmNames (numRealParms))
  allocate (realParmValues (numRealParms))
  allocate (realParmChanged (numRealParms))

  allocate (intParmNames (numIntParms))
  allocate (intParmValues (numIntParms))
  allocate (intParmChanged (numIntParms))

  allocate (strParmNames (numStrParms))
  allocate (strParmValues (numStrParms))
  allocate (strParmChanged (numStrParms))

  allocate (logParmNames (numLogParms))
  allocate (logParmValues (numLogParms))
  allocate (logParmChanged (numLogParms))

  call RuntimeParameters_getAllInt(numIntParms, &
                                 intParmNames, intParmValues, intParmChanged)

  call RuntimeParameters_getAllReal(numRealParms, &
                                 realParmNames, realParmValues, realParmChanged)

  call RuntimeParameters_getAllStr(numStrParms, &
                                 strParmNames, strParmValues, strParmChanged)

  call RuntimeParameters_getAllLog(numLogParms, &
                                 logParmNames, logParmValues, logParmChanged)

  do i = 1, numIntParms
    call write_tau_metadata_int("parameter " // intParmNames(i), &
                                 intParmValues(i))
  enddo

  do i = 1, numRealParms
    call write_tau_metadata_real("parameter " // realParmNames(i), &
                                 realParmValues(i))
  enddo

  do i = 1, numStrParms
    call write_tau_metadata_str("parameter " // strParmNames(i), &
                                 strParmValues(i))
  enddo

  do i = 1, numLogParms
    call write_tau_metadata_log("parameter " // logParmNames(i), &
                                 logParmValues(i))
  enddo

  deallocate(intParmNames)
  deallocate(intParmValues)
  deallocate(intParmChanged)

  deallocate(realParmNames)
  deallocate(realParmValues)
  deallocate(realParmChanged)

  deallocate(strParmNames)
  deallocate(strParmValues)
  deallocate(strParmChanged)

  deallocate(logParmNames)
  deallocate(logParmValues)
  deallocate(logParmChanged)

!==============================================================================

end subroutine Timers_init
