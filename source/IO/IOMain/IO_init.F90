!!****if* source/IO/IOMain/IO_init
!!
!! NAME
!!  IO_init
!!
!! SYNOPSIS
!!
!!  IO_init()
!!
!! DESCRIPTION
!!
!!  Perform IO initialization, which includes:
!!  getting the runtime parameters and if it is a restart
!!  calling IO_readCheckpoint.
!!
!!
!!  The IO unit uses a number of runtime parameters to determine
!!  if and when various types of output files need to be written.
!!  The IO unit writes checkpoint(restart) files, plotfiles for
!!  visualization, particle plotfiles, and .dat files which hold
!!  the diagnostic data like total energy, pressure, mass etc.
!!
!!  To determine exactly which runtime parameters control these
!!  files, please check the Config file in IO/IOMain or the 
!!  setup_params file in the object directory.
!!
!!
!! 
!! ARGUMENTS
!!
!!
!!
!!*** 
subroutine IO_init()

  use IO_data
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getComm,&
       Driver_getMype, Driver_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_getNumReal, RuntimeParameters_getNumInt, &
    RuntimeParameters_getNumStr, RuntimeParameters_getNumLog, RuntimeParameters_set
  use Simulation_interface, ONLY : Simulation_mapStrToInt,&
                                   Simulation_mapIntToStr
  use IO_interface, ONLY : IO_readCheckpoint,IO_getPrevScalar
  use Cosmology_interface, ONLY : Cosmology_getRedshift
  use io_ptInterface, ONLY : io_ptInit
  use Logfile_interface, ONLY:  Logfile_stamp
  use Grid_interface, only: Grid_formatNonRep

  implicit none
#include "Flash.h"
#include "constants.h"
 include "Flash_mpi.h"


  integer :: i, j, blocksPerFile, color, key, ierr, posBlank, collectiveHDF5
  character (len=MAX_STRING_LENGTH) :: plotVarNum

!!  save attribute is necessary for the XLF compiler when the variables are
!!  "OUT" argument for an overloaded function
  integer, save :: step
  real,save :: simTime,currentRedshift 

  integer :: error !catch errors in IO_getPrevScalar calls
  character(len=4) :: isErr
  logical :: havePlotVars
  
  integer, parameter :: nonrep_locunk1(0:NONREP_COUNT) = NONREP_LOCUNK1
  integer :: nonrep_maxlocs(0:NONREP_COUNT) = NONREP_MAXLOCS ! this should be a parameter but that was cuasing gfortran to choke
  character(len=*), parameter :: nonrep_rpcount_flat = NONREP_RPCOUNT_FLAT
  integer, parameter :: nonrep_rpcount_start(1:NONREP_COUNT+1) = NONREP_RPCOUNT_START
  integer :: nonrep, nonrep_globs(0:NONREP_COUNT)
  
  io_forcedFilename = "forced"

  io_lastWallClockCheckpoint = MPI_Wtime()
  io_CPUSeconds = 0.0

  call Driver_getMype(GLOBAL_COMM, io_globalMe)
  call Driver_getMype(MESH_COMM, io_meshMe)
  call Driver_getMype(MESH_ACROSS_COMM, io_acrossMe)
  
  call Driver_getNumProcs(GLOBAL_COMM, io_globalNumProcs)
  call Driver_getNumProcs(MESH_COMM, io_meshNumProcs)  
  call Driver_getNumProcs(MESH_ACROSS_COMM, io_acrossNumProcs)

  call Driver_getComm(GLOBAL_COMM,io_globalComm)

  call RuntimeParameters_get('alwaysRestrictCheckpoint', io_alwaysRestrictCheckpoint)
  call RuntimeParameters_get('alwaysComputeUserVars', io_alwaysComputeUserVars)
  call RuntimeParameters_get('ignoreForcedPlot', io_ignoreForcedPlot)
  call RuntimeParameters_get('plotFileNumber', io_plotFileNumber)
  call RuntimeParameters_get('checkpointFileNumber', io_checkpointFileNumber)
  call RuntimeParameters_get('forcedPlotFileNumber', io_forcedPlotFileNumber)

  call RuntimeParameters_get('plotFileIntervalTime',    io_plotFileIntervalTime)
  call RuntimeParameters_get('plotFileIntervalStep',    io_plotFileIntervalStep)
  call RuntimeParameters_get('plotFileIntervalZ', io_plotFileIntervalZ)
  call RuntimeParameters_get('checkpointFileIntervalTime',   io_checkpointFileIntervalTime)
  call RuntimeParameters_get('checkpointFileIntervalStep',   io_checkpointFileIntervalStep)
  call RuntimeParameters_get('checkpointFileIntervalZ', io_checkpointFileIntervalZ)
  

  call RuntimeParameters_get('tinitial', io_tinitial)
  
  call RuntimeParameters_get('restart',  io_restart)
  call RuntimeParameters_get('rolling_checkpoint', io_rollingCheckpoint)
  call RuntimeParameters_get('wall_clock_checkpoint', io_wallClockCheckpoint)
  call RuntimeParameters_get('plotfileMetadataDP', io_plotfileMetadataDP)
  call RuntimeParameters_get('plotfileGridQuantityDP', io_plotfileGridQuantityDP)
  call RuntimeParameters_get('memory_stat_freq', io_memoryStatFreq)
  call RuntimeParameters_get('rss_limit', io_maxRSS)
  !------------------------------------------------------------------------------
  ! Dump out memory usage statistics as soon as we know if we are monitoring them
  !------------------------------------------------------------------------------
  if (io_memoryStatFreq > 0) call io_memoryReport()

  call RuntimeParameters_get('wr_integrals_freq', io_integralFreq)


  ! get the runtime parameters
  call RuntimeParameters_get('stats_file', io_statsFileName)  
  call RuntimeParameters_get('basenm', io_baseName)
  !! if basenm is changed from default and stats_file isn't, then make it conform
  if ((io_baseName .NE. 'flash_') .AND.                    &
 &      (io_statsFileName .EQ. "flash.dat")) then
     posBlank = index(io_baseName,' ')
     io_statsFileName = io_baseName(:posBlank-2) // '.dat'  ! Hopefully remove trailing _
     call RuntimeParameters_set('stats_file',io_statsFileName)
  endif
  call RuntimeParameters_get('io_writeMscalarIntegrals', io_writeMscalarIntegrals)
 
  call RuntimeParameters_get("output_directory", io_outputDir)
  call RuntimeParameters_get('geometry', io_geometry)
  call RuntimeParameters_get('chkGuardCellsInput', io_chkGuardCellsInput)
  call RuntimeParameters_get('chkGuardCellsOutput', io_chkGuardCellsOutput)
!!#ifndef FIXEBLOCKSIZE
!!$  call RuntimeParameters_get('iguard', io_iguard)
!!$  call RuntimeParameters_get('jguard', io_jguard)
!!$  call RuntimeParameters_get('kguard', io_kguard)
   io_iguard = NGUARD
   io_jguard = NGUARD
   io_kguard = NGUARD
  !accounting for various dimension problems
  io_jguard = io_jguard*K2D
  io_kguard = io_kguard*K3D
  

!!#endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! generate the lists of local and global unks
  io_unkNonRep(:) = 0
  io_unkNonRepIdx(:) = 0
  io_unkActive(:) = .true.
  io_unkToGlobal(:) = 0
  nonrep_globs(:) = 0
  do nonrep=1, NONREP_COUNT
     call RuntimeParameters_get( &
        nonrep_rpcount_flat(nonrep_rpcount_start(nonrep):nonrep_rpcount_start(nonrep+1)-1), &
        nonrep_globs(nonrep))
     do i=1, nonrep_maxlocs(nonrep)
        j = nonrep_locunk1(nonrep) - 1 + i
        io_unkNonRep(j) = nonrep
        io_unkNonRepIdx(j) = NONREP_LOC2GLOB(i, io_acrossMe, io_acrossNumProcs)
        io_unkActive(j) = io_unkNonRepIdx(j) <= nonrep_globs(nonrep)
        if(.not. io_unkActive(j)) io_unkNonRepIdx(j) = 0
     end do
  end do
  ! the size of the global unks = (global nonrep unks) + (local unks) - (local nonrep unks)
  allocate(io_unklabelsGlobal(sum(nonrep_globs) + UNK_VARS_END-UNK_VARS_BEGIN+1 - sum(nonrep_maxlocs)))
  j = 1 ! keeps track of our position in io_unklabelsGlobal
  ! start with the regular (nonrep) unks
  do i=UNK_VARS_BEGIN, UNK_VARS_END
     call Simulation_mapIntToStr(i, io_unklabels(i), MAPBLOCK_UNK)
     if(io_unkNonRep(i) == 0) then ! if this unk doesnt belong to a nonrep array, then store it in the global list
        io_unklabelsGlobal(j) = io_unklabels(i)
        io_unkToGlobal(i) = j
        j = j + 1
     end if
  end do
  ! now add the nonrep unks
  do nonrep=1, NONREP_COUNT
     do i=1, nonrep_globs(nonrep) ! generate the name for each element of this array
        call Grid_formatNonRep(nonrep, i, io_unklabelsGlobal(j))
        if(NONREP_MESHOFGLOB(i,io_acrossNumProcs) == io_acrossMe) then
           io_unkToGlobal(nonrep_locunk1(nonrep)-1 + NONREP_GLOB2LOC(i,io_acrossMe,io_acrossNumProcs)) = j
        end if
        j = j + 1
     end do
  end do
  if(j-1 .ne. ubound(io_unklabelsGlobal,1)) call Driver_abortFlash('BUG AT: ' // FILE_AT_LINE)

#if defined(IO_HDF5_PARALLEL) || (defined(FLASH_IO_PNETCDF) && defined(FLASH_IO_EXPERIMENTAL))
#else
  if(io_acrossNumProcs > 1) then
     call Driver_abortFlash('[' // FILE_AT_LINE // '] ERROR: ' // &
        'mesh replication is only supported by parallel HDF5 ' //&
        'and derived datatype pnetCDF I/O units.')
  end if
  if(any(io_unkToGlobal == 0)) then

     call Logfile_stamp('ERROR: all unks allocated for a NONREP variable array must be in use.')
     call Logfile_stamp('This restriction does not apply to the PM, UG, and NOFBS implementations')
     call Logfile_stamp('of IO/IOMain/hdf5/parallel.')
     call Logfile_stamp('')
     call Logfile_stamp('Possible Solutions:')
     call Logfile_stamp('1) Try using parallel IO (+parallelIO in setup line)')
     call Logfile_stamp('2) Remove NONREP config file commands')
     call Logfile_stamp('3) Make sure that number of unks in use matches the global NONREP array size.')
     call Logfile_stamp('   For example, when using, MGD, make sure that rt_mgdNumGroups is equal to ')
     call Logfile_stamp('   meshCopyCount*{number of local MGDR NONREPs} This number comes from the MGDR ')
     call Logfile_stamp('   NONREP command in Config File')

     call Driver_abortFlash('[' // FILE_AT_LINE // '] ERROR: Unused UNK ' // & 
          'due to NONREP Config directive and IO unit. LOOK AT THE LOG FILE.')
  end if
#endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! get all the IO plot vars.  These are the variables that will be written to plotfiles

  !check if we have a phantom variable
  if(NUNK_VARS == 1) then
    call Simulation_mapIntToStr(1,isErr, MAPBLOCK_UNK)
    havePlotVars = .not.(isErr .EQ. 'dumm' .or. isErr .EQ. 'err ')
  else
    havePlotVars = .TRUE.
  end if

  if(havePlotVars) then
     do i = 1,io_maxPlotVars
        write(plotVarNum, *) i
        plotVarNum = ADJUSTL(plotVarNum)
        plotVarNum = 'plot_var_' // plotVarNum
        call RuntimeParameters_get(plotVarNum, io_plotVarStr(i))
        call makeLowercase(io_plotVarStr(i))
     end do !i = 1,io_maxPlotVar
  end if

  !translate the plotvars from strings to integers
  !and calculate the number of them
  j = 0
  do i=1, io_maxPlotVars
     if(any(io_plotVarStr(i) == io_unklabelsGlobal)) then
        j = j + 1
        if (j < i) io_plotVarStr(j) =  io_plotVarStr(i)
        call Simulation_mapStrToInt(trim(io_plotVarStr(j)), io_plotVar(j), MAPBLOCK_UNK)
     else                       ! Not a recognized name? Then do not count it.
        if (io_plotVarStr(i) .NE. "none") &
             call Logfile_stamp(io_plotVarStr(i),"[IO_init] ignoring unrecognized plot variable")
     end if
  end do  
  io_nPlotVars = j


  if(NFACE_VARS .GT. 0) then
    do i= 1, NFACE_VARS
      WRITE(io_faceXVarLabels(i),'(''fcx'',I1)') i
      WRITE(io_faceYVarLabels(i),'(''fcy'',I1)') i
      WRITE(io_faceZVarLabels(i),'(''fcz'',I1)') i
    end do
  endif


  ! get all the IO scratch plot grid vars.  These are the scratch grid vars to
  ! be written to the checkpoint and plotfiles 
  call RuntimeParameters_get("plot_grid_var_1", io_plotGridVarStr(1))
  call RuntimeParameters_get('plot_grid_var_2', io_plotGridVarStr(2))
  call RuntimeParameters_get('plot_grid_var_3', io_plotGridVarStr(3))
  call RuntimeParameters_get('plot_grid_var_4', io_plotGridVarStr(4))
  call RuntimeParameters_get('plot_grid_var_5', io_plotGridVarStr(5))
  call RuntimeParameters_get('plot_grid_var_6', io_plotGridVarStr(6))
  call RuntimeParameters_get('plot_grid_var_7', io_plotGridVarStr(7))
  call RuntimeParameters_get('plot_grid_var_8', io_plotGridVarStr(8))
  call RuntimeParameters_get('plot_grid_var_9', io_plotGridVarStr(9))
  call RuntimeParameters_get('plot_grid_var_10', io_plotGridVarStr(10))
  call RuntimeParameters_get('plot_grid_var_11', io_plotGridVarStr(11))
  call RuntimeParameters_get('plot_grid_var_12', io_plotGridVarStr(12))


  !translate the scratch grid vars from strings to integers
  !and calculate the number of them

  io_nPlotGridVars = 0

  do i=SCRATCH_GRID_VARS_BEGIN,SCRATCH_GRID_VARS_END
     call Simulation_mapIntToStr(i, io_scratchGridVarlabels(i),MAPBLOCK_SCRATCH)
  end do

  do i=1, io_maxPlotGridVars
     call Simulation_mapStrToInt(trim(io_plotGridVarStr(i)), io_plotGridVar(i),MAPBLOCK_SCRATCH)
     if(io_plotGridVar(i) /= NONEXISTENT) then
        print *, "got the plotvarstr ", io_plotGridVarStr(i) 
        io_nPlotGridVars = io_nPlotGridVars + 1
     end if
  end do




  !get the runtime parameter to indicate if bytepacking is turned on or
  !off in plotfiles
  call RuntimeParameters_get('bytePack', io_bytePack)


  !get the runtime parameter for the number of files to split
  !and hdf5 file into
  call RuntimeParameters_get('outputSplitNum', io_outputSplitNum)
 
  if((mod(io_meshNumProcs, io_outputSplitNum) /= 0)) then
     call Driver_abortFlash('mod(io_meshNumProcs, io_outputSplitNum) /= 0!')
  end if

  !set each files specific split number

  !prepare the io communicator if necessary
  call Driver_getComm(MESH_COMM, io_meshComm)
  io_comm = io_globalComm !io_meshComm
  call MPI_Comm_rank(io_meshComm, io_meshMe, ierr)
  call MPI_Comm_size(io_meshComm, io_meshNumProcs, ierr)
  !io_splitFileNum = io_globalMe / io_meshNumProcs
  if(io_outputSplitNum /= 1) then
     !This line can only be good for UG. Paramesh figures this differently
     blocksPerFile = io_meshNumProcs / io_outputSplitNum
     color = io_meshMe /blocksPerfile
     key = mod(io_meshMe, blocksPerFile)
     call MPI_Comm_split(io_comm,color,key,io_comm,ierr)
     io_splitFileNum = color
  endif
  !io_outputSplitNum=io_outputSplitNum*io_globalNumProcs/io_meshNumProcs
  !print*,'io_outputSplitNum',io_outPutSplitNum, io_meshNumProcs, io_comm

  call RuntimeParameters_get("fileFormatVersion", io_fileFormatVersion)

  !useCollectiveHDF5=.true. switches on collective I/O optimizations.
  !The call to io_h5set_xfer_mode must happen before IO_readCheckpoint.
  call RuntimeParameters_get("useCollectiveHDF5", io_useCollectiveHDF5)

#ifdef IO_HDF5_PARALLEL
#if defined(FLASH_GRID_UG) && \
  !defined(FIXEDBLOCKSIZE) && \
  !defined(FLASH_IO_EXPERIMENTAL)

  !We need to switch off collective I/O optimizations for nofbs UG.
  !This is because only the master process writes metadata to file.
  !This leads to a deadlock in the io_h5writeXXX functions as they
  !do not enable a process to contribute zero data.  The extra code
  !that allows a process to contribute zero data exists in io_h5_xfer.c.
  !This single function replaces the dozen or so io_h5writeXXX
  !functions and is used by experimental I/O implementation.  As the
  !name suggests a lot of this work is still slightly experimental, so
  !it is not used by default.  To switch on collective I/O optimizations
  !for nofbs UG use experimental I/O with the following setup options:
  !./setup ... -nofbs +ug -unit=IO/IOMain/hdf5/parallel/PM_argonne
  if (io_useCollectiveHDF5) then
     io_useCollectiveHDF5 = .false.
     call Logfile_stamp( &
          "*** Switching off collective I/O for NOFBS grid ***", &
          "[IO_init]")
  end if
#endif

  !This function only makes sense in parallel HDF5 mode.
  collectiveHDF5 = 0
  if(io_useCollectiveHDF5) collectiveHDF5 = 1
  call io_h5set_xfer_mode(collectiveHDF5)
#endif


#ifdef FLASH_IO_EXPERIMENTAL
  !Similarly the MPI types must be set up before we read a checkpoint.
  !Only needed when we use type based (i.e. experimental) I/O.
  call io_typeInit(io_globalMe, io_globalNumProcs)
#endif


  call RuntimeParameters_getNumReal(io_numRealParms)
  call RuntimeParameters_getNumInt(io_numIntParms)
  call RuntimeParameters_getNumStr(io_numStrParms)
  call RuntimeParameters_getNumLog(io_numLogParms)

  if(io_restart) then
     !if we are restarting from a checkpoint file
     !initialize num parameters
     io_numRealParmsPrev = io_numRealParms + io_prevParmPad
     io_numIntParmsPrev = io_numIntParms + io_prevParmPad
     io_numStrParmsPrev = io_numStrParms + io_prevParmPad
     io_numLogParmsPrev = io_numLogParms + io_prevParmPad

     io_numRealScalars = io_maxParms
     io_numIntScalars = io_maxParms
     io_numStrScalars = io_maxParms
     io_numLogScalars = io_maxParms

     call RuntimeParameters_get('checkpointFileNumber', io_checkpointFileNumber)
     call IO_readCheckpoint()
     
  end if

  if(io_restart) then
     call IO_getPrevScalar("nstep", step)
     call IO_getPrevScalar("time", simTime)
     call IO_getPrevScalar("redshift", currentRedshift, error)
     !we can survive this error.
     if(error /= NORMAL) then
        if(error == NOTFOUND) then
           currentRedshift = 0.0
        else
           call Driver_abortFlash("ERROR: Error in lookup of 'redshift' scalar.")
        end if
     end if
     ! Values from the checkpoint file for next checkpoint and plot times
     ! will be used only if they are compatible with the current
     ! checkpointFileIntervalTime and plotfileIntervalTime, repectively
     ! (which may have changed from those of the original run) -
     ! see code below. - KW

     call IO_getPrevScalar("nextCheckpointTime", io_nextCheckpointTime)
     call IO_getPrevScalar("nextCheckpointZ", io_nextCheckpointZ, error)
     !we can survive this error.
     if(error /= NORMAL) then
        if(error == NOTFOUND) then
           io_nextCheckpointZ = HUGE(1.)
        else
           call Driver_abortFlash("ERROR: Error in lookup of nextCheckpointZ scalar.")
        end if
     end if
     
     call IO_getPrevScalar("nextPlotfileTime", io_nextPlotfileTime)
     call IO_getPrevScalar("nextPlotFileZ", io_nextPlotFileZ, error)
     !we can survive this one, too.
     if(error /= NORMAL) then
        if(error == NOTFOUND) then
           io_nextPlotFileZ = HUGE(1.)
        else
           call Driver_abortFlash("ERROR: Error in lookup of nextPlotFileZ scalar.")
        end if
     end if
  else
     call RuntimeParameters_get("nbegin", step)
     call RuntimeParameters_get("tinitial",simTime)
     
     io_nextCheckpointTime = 0.0 !any number to make sure logical expression below is valid, overwritten below
     io_nextPlotfileTime = 0.0   !any number to make sure logical expression below is valid, overwritten below
     io_nextPlotFileZ = 0.0     !prevent uninitialized data errors when runtime checking for such is enabled
     io_splitNumBlks = 0
     io_splitParts = 0
  end if
  
  !initialize the next step at which point to 
  !write a checkpoint file 
  io_nextCheckPointStep = step + io_checkpointFileIntervalStep
 
  !calculate next simulation time in which to 
  !write a checkpoint file
  if (.not. io_restart .or. (io_nextCheckpointTime <= simTime) .or. &
       (io_nextCheckpointTime > simTime + io_checkpointFileIntervalTime)) then
     if ( io_checkpointFileIntervalTime > 0.e0 ) then
        io_nextCheckPointTime = simTime + io_checkpointFileIntervalTime
     else
        io_nextCheckpointTime = 0.e0
     end if
  end if
  
 
  
  if(.not. io_restart) call Cosmology_getRedshift(currentRedshift)
  if(.not. io_restart )then !.or. (io_nextCheckpointZ >= currentRedshift) .or. &
       !(io_nextCheckpointZ < currentRedshift +io_checkpointFileIntervalZ)) then
     if(io_checkpointFileIntervalZ < HUGE(1.)) then
        io_nextCheckpointZ = currentRedshift - io_checkpointFileIntervalZ
     else
        io_nextCheckpointZ = HUGE(1.)
     end if
  end if

  !initialize next step at which to write a plotfile
  io_nextPlotFileStep = step + io_plotfileIntervalStep

  !calculate next simulation time in which to 
  !write a plotfile file
  if (.not. io_restart .or. (io_nextPlotFileTime <= simTime) .or. &
       (io_nextPlotFileTime > simTime + io_plotFileIntervalTime)) then
     if ( io_plotFileIntervalTime > 0.e0 ) then
        io_nextPlotFileTime = simTime + io_plotFileIntervalTime

     else
        io_nextPlotFileTime = 0.e0
     end if
  end if


  call RuntimeParameters_get('typeMatchedXfer', io_type_matched_xfer)
  call RuntimeParameters_get ("reduceGcellFills", io_reduceGcellFills)
  call RuntimeParameters_get ("summaryOutputOnly", io_summaryOutputOnly)

  !initialize runtime parameters specific to IOParticles subunit
  !if particles are not included in a simulation this is a stub
  call io_ptInit()


  !------------------------------------------------------------------------------
  ! Dump out memory usage statistics again if we are monitoring them
  !------------------------------------------------------------------------------
  if (io_memoryStatFreq > 0) call io_memoryReport()

  !Set the initial state of io_outputInStack.  This variable is used by 
  !io_restrictBeforeWrite.F90.  It is .false. for most of the simulation, 
  !and is only briefly .true. when an IO_output routine is in the call stack.
  io_outputInStack = .false.

  if (io_summaryOutputOnly) then
     !Ensure that the variables nextCheckpointTime, nextCheckpointZ, 
     !nextPlotfileTime and nextPlotFileZ retain their original values so that
     !they will contain expected values if we write an emergency checkpoint
     !file, e.g. with .dump_restart.

     !IO_output.F90 plot file conditions are shown in parentheses:
     !------------------------------------------------------------
     !(nstep == io_nextPlotFileStep)
     io_nextPlotFileStep = -1 !Normally set from io_plotFileIntervalStep
     !(io_plotFileIntervalTime > 0.e0 .AND. simTime >= io_nextPlotFileTime)
     io_plotFileIntervalTime = -1.0 !No need to modify io_nextPlotFileTime
     !(io_plotFileIntervalZ < HUGE(1.) .AND. currentRedshift <= io_nextPlotFileZ)
     io_plotFileIntervalZ = HUGE(1.)

     !IO_output.F90 checkpoint file conditions are shown in parentheses:
     !------------------------------------------------------------------
     !(nstep == io_nextCheckpointStep)
     io_nextCheckpointStep = -1 !Normally set from io_checkpointFileIntervalStep
     !(io_checkpointFileIntervalTime > 0.e0  .and. simTime >= io_nextCheckpointTime)
     io_checkpointFileIntervalTime = -1.0 !No need to modify io_nextCheckpointTime
     !(io_checkpointFileIntervalZ < HUGE(1.) .and. currentRedshift <= io_nextCheckpointZ)
     io_checkpointFileIntervalZ = HUGE(1.)
     !(io_wallClockCheckpoint > 0.e0)
     io_wallClockCheckpoint = -1.0
     !(io_maxRSS > 0.0)
     io_maxRSS = -1.0
  end if

end subroutine IO_init
