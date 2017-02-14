!!****if* source/monitors/Logfile/LogfileMain/Logfile_create
!!
!! NAME
!!     Logfile_create
!! 
!! SYNOPSIS
!!     Logfile_create()
!!
!! DESCRIPTION
!!
!!     Creates the named log file and writes some header information
!!     to it, including the included units, runtime parameters,
!!     physical constants.  Meta data about the run is also stored
!!     like a time stamp, the run dimensionality, compiler flags and
!!     more.  The logfile can be stamped from any unit to store the
!!     simulation's progress.
!!
!!     The name of the parameter file is taken as an input; it is
!!     echoed to the log file.  Only the master processor actually
!!     writes anything.  In order to avoid accidentally overwriting an
!!     important logfile during a science run, the logfile is always
!!     opened in append mode.
!!
!!     
!!
!! ARGUMENTS
!!     
!!
!! NOTES
!!  variables that begin with "log_" are defined in the fortran 
!!  module Logfile_data.  The prefix "log_" is meant to indicate
!!  that these variables have Logfile unit scope.  Other variables
!!  are local to the individual subroutines
!!
!!***

subroutine Logfile_create ()

  !$ use omp_lib
  use Logfile_data, ONLY : log_fileOpen, log_lun, log_fileName, &
       log_runComment, log_runNum, log_globalMe, log_globalNumProcs
  use Driver_interface, ONLY : Driver_abortFlash, Driver_putTimeStamp, &
       Driver_mpiThreadSupport
  use RuntimeParameters_interface, ONLY : &
    RuntimeParameters_getNumReal, RuntimeParameters_getNumInt, &
    RuntimeParameters_getNumIgn, RuntimeParameters_getNumStr, RuntimeParameters_getNumLog, &
    RuntimeParameters_getAllInt, RuntimeParameters_getAllReal, RuntimeParameters_getAllStr, & 
    RuntimeParameters_getAllLog, RuntimeParameters_stampIgnored
  use Logfile_interface, ONLY : Logfile_close, Logfile_break, Logfile_stampMessage
  use PhysicalConstants_interface, ONLY : &
    PhysicalConstants_listUnits, PhysicalConstants_list
  use Multispecies_interface, ONLY : Multispecies_list

  implicit none

#include "constants.h"
#include "Flash.h"
  include "Flash_mpi.h"

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
  character(len=40), save                :: dateStr
  character(len=5) ,save                 :: parmNameFormat

  integer :: ioStat, numUnits, numIgnoredParms
  integer :: mpiVersion, mpiSubversion, ierr
  logical :: mpiThreadSupport, isOpenmpMacroDefined


  !Only the master processor creates and initializes the log file.
    
  if (log_globalMe .eq. MASTER_PE) then
     
     ! The timestamp in the log file should be echoed to the .dat file
     call current_date_time(dateStr)
     call Driver_putTimeStamp(dateStr)  
    
     call Logfile_close()    ! in case it's already open
     ! open the file -- use status=unknown so if we are restarting
     ! or re-running, we can append to the end of it
     open (log_lun, file=log_fileName, iostat=ioStat,  & 
          status='unknown', position='append')

     if (ioStat == 0) then
        log_fileOpen = .true.
        write (log_lun,*) 'FLASH log file:  ', & 
             dateStr(1:len_trim(dateStr)), '    ', &  ! 
             'Run number:  ', log_runNum(1:len_trim(log_runNum))

        call Logfile_break("=")

        write (log_lun,*) 'Number of MPI tasks:       ', log_globalNumProcs
        call MPI_Get_version(mpiVersion, mpiSubversion, ierr)
        write (log_lun,*) 'MPI version:               ', mpiVersion
        write (log_lun,*) 'MPI subversion:            ', mpiSubversion

        !Write these stamps only when we are compiling with OpenMP support.
        !$ call Driver_mpiThreadSupport(mpiThreadSupport)
        !$ write (log_lun,*) 'MPI thread support:                  ', &
        !$  mpiThreadSupport

        !$omp parallel
        !$ if (omp_get_thread_num() == 0) then
        !$   write (log_lun,*) 'OpenMP threads/MPI task:   ', &
        !$    omp_get_num_threads()
        !$ end if
        !$omp end parallel
        
        !openmp_version is defined in omp_lib module (see OpenMP standard).
        !$ write (log_lun,*) 'OpenMP version:                 ', &
        !$  openmp_version

        !I have found that Absoft 64-bit Pro Fortran 11.1.3 for Linux x86_64
        !does not define _OPENMP when compiling Fortran files.  We have
        !some code that is conditionally compiled in FLASH based on _OPENMP       
        !If it is not defined you can do something like this in Makefile.h
        !OPENMP_FORTRAN = -openmp -D_OPENMP=200805
#ifdef _OPENMP
        isOpenmpMacroDefined = .true.
#else
        isOpenmpMacroDefined = .false.
#endif
        !$ write (log_lun,*) 'Is "_OPENMP" macro defined:          ', &
        !$  isOpenmpMacroDefined

        write (log_lun,*) 'Dimensionality:            ', NDIM
        write (log_lun,*) 'Max Number of Blocks/Proc: ', MAXBLOCKS
#ifdef FIXEDBLOCKSIZE
        write (log_lun,*) 'Number x zones:            ', NXB
        write (log_lun,*) 'Number y zones:            ', NYB
        write (log_lun,*) 'Number z zones:            ', NZB
#else
        write (log_lun,*) 'Number x zones:            (igridsize / iprocs)'
        if (NDIM > 1) then
           write (log_lun,*) 'Number y zones:            (jgridsize / jprocs)'
        else
           write (log_lun,*) 'Number y zones:            ', NYB !should be 1
        end if
        if (NDIM > 2) then
           write (log_lun,*) 'Number z zones:            (kgridsize / kprocs)'
        else
           write (log_lun,*) 'Number z zones:            ', NZB !should be 1
        end if
#endif

        call setup_buildstamp (names, names2, MAX_STRING_LENGTH)
        write (log_lun,*) "Setup stamp:     ", trim(names)
        write (log_lun,*) "Build stamp:     ", trim(names2)
        
        call setup_systemInfo (names, MAX_STRING_LENGTH)
        write (log_lun,*) "System info:     ", trim(names)
        
        write (log_lun,*) 'Version:         ', setup_flashRelease()
        
        call setup_buildstats(buildDate, buildDir, buildMachine, setupCall, &
             cflags, fflags)


        call removeNullChar(buildDate)
        call removeNullChar(buildDir)
        call removeNullChar(buildMachine)
        call removeNullChar(setupCall)
        call removeNullChar(cflags)
        call removeNullChar(fflags)
        

    
        write (log_lun,*) 'Build directory: ', trim(buildDir)
        write (log_lun,*) 'Setup syntax:    ', trim(setupCall)
        
        
        write (log_lun,*) 'f compiler flags: ', trim(fflags)
        write (log_lun,*) 'c compiler flags:    ', trim(cflags)
        
        
        call Logfile_break("=")
        write (log_lun,*) 'Comment:  ', log_runComment(1:len_trim(log_runComment))
        call Logfile_break("=")
        
        write (log_lun,*) 'FLASH Units used:'
        

        call setup_getNumFlashUnits(numUnits)
        allocate(s1(numUnits))
        call log_getUnitsArr(numUnits, s1)
        
        do i = 1, numUnits
           write (log_lun,*) '  ', trim(s1(i))
        enddo
        deallocate(s1)
        
        call Logfile_break("=")
        

        write (log_lun,*) 'RuntimeParameters:'
        call Logfile_break(" ")
        call Logfile_break("=")

        if (.true.) then

           ! get number of runtime parameters        
           call RuntimeParameters_getNumReal(numRealParms)
           call RuntimeParameters_getNumInt(numIntParms)
           call RuntimeParameters_getNumStr(numStrParms)
           call RuntimeParameters_getNumLog(numLogParms)

           !! allocate the space for the Runtime Parameters
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

           ! get runtime parameters
           call RuntimeParameters_getAllInt(numIntParms, &
                intParmNames, intParmValues, intParmChanged)
           
           call RuntimeParameters_getAllReal(numRealParms, &
                realParmNames, realParmValues, realParmChanged)
           
           call RuntimeParameters_getAllStr(numStrParms, &
                strParmNames, strParmValues, strParmChanged)
           
           call RuntimeParameters_getAllLog(numLogParms, &
                logParmNames, logParmValues, logParmChanged)


           !write runtime parameters
999        format('(A',I2,',')
           do i=1, numIntParms
              if (len_trim(intParmNames(i)) < 27) then
                 write(parmNameFormat,999) 27
              else
                 write(parmNameFormat,999) len_trim(intParmNames(i))
              endif
              if (intParmChanged(i)) then
                 write (log_lun, parmNameFormat // "' = ', I10, ' [CHANGED]')") intParmNames(i), intParmValues(i) 
              else if (intParmValues(i) == -HUGE(1)) then
                 write (log_lun, parmNameFormat // "' = -HUGE(1)')") intParmNames(i)
              else
                 write (log_lun, parmNameFormat // "' = ', I10)") intParmNames(i), intParmValues(i) 
              endif
           end do              

           do i=1, numRealParms
              if (len_trim(realParmNames(i)) < 27) then
                 write(parmNameFormat,999) 27
              else
                 write(parmNameFormat,999) len_trim(realParmNames(i))
              endif
              if (realParmChanged(i)) then
                 write (log_lun, parmNameFormat // "' = ', E25.3, ' [CHANGED]')") realParmNames(i), realParmValues(i) 
              else
                 write (log_lun, parmNameFormat // "' = ', E25.3)") realParmNames(i), realParmValues(i) 
              endif
           end do              

           do i=1, numStrParms
              if (len_trim(strParmNames(i)) < 27) then
                 write(parmNameFormat,999) 27
              else
                 write(parmNameFormat,999) len_trim(strParmNames(i))
              endif
              if (strParmChanged(i)) then
                 write (log_lun, parmNameFormat // "' = ', A30, ' [CHANGED]')") strParmNames(i), strParmValues(i) 
              else
                 write (log_lun, parmNameFormat // "' = ', A30)") strParmNames(i), strParmValues(i) 
              endif
           end do              

           do i=1, numLogParms
              if (len_trim(logParmNames(i)) < 27) then
                 write(parmNameFormat,999) 27
              else
                 write(parmNameFormat,999) len_trim(logParmNames(i))
              endif
              if (logParmChanged(i)) then
                 write (log_lun, parmNameFormat // "' = ', L2, ' [CHANGED]')") logParmNames(i), logParmValues(i) 
              else
                 write (log_lun, parmNameFormat // "' = ', L2)") logParmNames(i), logParmValues(i) 
              endif
           end do              

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
           

        end if


        call RuntimeParameters_getNumIgn(numIgnoredParms)
        if (numIgnoredParms > 0) then
           call Logfile_break(" ")
           call Logfile_break("=")
           call Logfile_break(" ")

           write (log_lun,*) 'WARNING: Ignored Parameters :'
           write (log_lun,*) 'These parameters were found in the flash.par file, but they were'
           write (log_lun,*) 'not declared in any Config file for the simulation!'


           call Logfile_break(" ")

           call RuntimeParameters_stampIgnored()

        end if

        call Logfile_break(" ")
        call Logfile_break("=")
        call Logfile_break(" ")

        write (log_lun,*) 'Known units of measurement:'
        call Logfile_break(" ")

        !write physical constants units
        call PhysicalConstants_listUnits(log_lun)

        call Logfile_break(" ")
        write (log_lun,*) 'Known physical constants:'
        call Logfile_break(" ")

        !write physical constants
        call PhysicalConstants_list(log_lun)
        call Logfile_break("=")


        call Logfile_break(" ")
#ifdef FLASH_MULTISPECIES
        write (log_lun,*) 'Multifluid database contents:'
        call Logfile_break(" ")

        !write multispecies
        call Multispecies_list(log_lun)
#else
        call Logfile_stampMessage('Multifluid database: not configured in')
#endif

        call Logfile_break(" ")
        call Logfile_break("=")


        close (log_lun)
        log_fileOpen = .false.

     else
        
        write (*,*) 'Logfile_create:  could not create log file'
        call Driver_abortFlash("Logfile_create : could not create logfile")

     endif

  endif

  return
end subroutine Logfile_create
