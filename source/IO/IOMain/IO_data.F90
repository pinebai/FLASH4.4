!!****if* source/IO/IOMain/IO_data
!!
!! NAME
!!  IO_data
!!
!! SYNOPSIS
!!
!!  use IO_data
!!
!! DESCRIPTION 
!!  
!!  Holds all the IO data that is needed by the IO unit     
!!
!! ARGUMENTS
!!
!!  none    
!!
!!
!!***

#include "constants.h"
#include "Flash.h"

module IO_data

  use nameValueLL_data

#ifdef USE_IO_C_INTERFACE
  use iso_c_binding
#else
  !Older compilers do not support Fortran 2003.
  integer, parameter :: c_char = KIND('A')
#endif

  integer, save :: io_checkpointFileNumber
  integer, save :: io_plotFileNumber
  integer, save :: io_forcedPlotFileNumber

  integer, save :: io_checkpointFileIntervalStep
  integer, save :: io_plotFileIntervalStep

  real, save :: io_checkpointFileIntervalTime
  real, save :: io_plotFileIntervalTime

  real, save :: io_checkpointFileIntervalZ
  real, save :: io_plotFileIntervalZ
  
  integer, save :: io_nextCheckpointStep
  integer, save :: io_nextPlotFileStep

  real, save :: io_nextCheckpointTime
  real, save :: io_nextPlotFileTime

  real, save :: io_nextCheckpointZ
  real, save :: io_nextPlotFileZ

  real, save :: io_wallClockCheckpoint
  real, save :: io_tinitial

  !split checkpoint files in io_outputSplitNum files
  integer, save :: io_outputSplitNum
  
  !number given to split file
  integer, save :: io_splitFileNum

  !checkpoint file io handle
  integer, save :: io_chkptFileID

  !communicator for io if necessary!
  integer, save :: io_comm = 0
  integer, save :: io_meshComm, io_meshMe, io_meshNumProcs
  integer, save :: io_acrossComm, io_acrossMe, io_acrossNumProcs
  integer, save :: io_globalComm=0, io_globalMe=0, io_globalNumProcs

  !used for split IO 
  integer, save :: io_splitNumBlks !what is the number of blocks in the file?
  integer, save :: io_splitParts !how man particles in a split file.

  integer, save :: io_rollingCheckpoint
  integer, save :: io_memoryStatFreq, io_integralFreq
  logical, save :: io_writeMscalarIntegrals = .FALSE.
  character(len=MAX_STRING_LENGTH), save :: io_geometry
  real, save :: io_redshift = 1.0
  real, save :: io_CPUSeconds, io_lastCPUSeconds, io_lastWallClockCheckpoint
  logical, save :: io_restart, io_bytePack, io_justCheckpointed
  
  !For clean exit when memory usage gets scary
  real, save :: io_maxRSS, io_measRSS

  integer,parameter :: io_maxParms = 500
  integer, parameter :: io_prevParmPad = 200

  character (len=MAX_STRING_LENGTH) :: io_parmStrName = "runtime parameters"
  character (len=MAX_STRING_LENGTH) :: io_scalarStrName = "scalars"


  type (context_type), save :: io_scalar  


  !!These are for the checkpoint file
  character (len=MAX_STRING_LENGTH), save :: io_baseName, io_outputDir

  !!This is to denote a "forced" plotfile
  character (len=MAX_STRING_LENGTH), save :: io_forcedFilename

  !!Store the integral quantities file name (the .dat file)
  character (len=MAX_STRING_LENGTH), save :: io_statsFileName


  !!These are for the linked list output
  integer, save :: io_numRealParms, io_numIntParms, io_numStrParms, io_numLogParms
  integer, save :: io_numRealParmsPrev, io_numIntParmsPrev
  integer, save :: io_numStrParmsPrev, io_numLogParmsPrev
  integer :: io_numRealScalars, io_numIntScalars, io_numStrScalars, io_numLogScalars

  integer, allocatable, target, save :: io_intParmValues(:), io_intScalarValues(:), io_intParmValuesPrev(:)
  real, allocatable, target, save :: io_realParmValues(:), io_realScalarValues(:) , io_realParmValuesPrev(:)
  logical, allocatable, target, save :: io_logParmValues(:), io_logScalarValues(:), io_logParmValuesPrev(:)
  character (len=MAX_STRING_LENGTH), allocatable, target, save :: &
    io_strParmValues(:), io_strScalarValues(:), io_strParmValuesPrev(:)
  integer, allocatable, target, save :: io_logToIntParmValues(:), io_logToIntScalarValues(:), io_logToIntParmValuesPrev(:)
  
  character (kind=c_char,len=MAX_STRING_LENGTH), allocatable, target, save :: &
    io_intParmNames(:), io_intScalarNames(:), io_intParmNamesPrev(:), &
    io_realParmNames(:), io_realScalarNames(:), io_realParmNamesPrev(:), &
    io_strParmNames(:), io_strScalarNames(:), io_strParmNamesPrev(:), &
    io_logParmNames(:), io_logScalarNames(:), io_logParmNamesPrev(:)

  !This is hard coded now, max number of vars to output to the plotfile
  integer, parameter :: io_maxPlotVars = MAX_PLOT_VARS
  integer, parameter :: io_maxPlotGridVars = 12 !max scratch grid vars to output
  integer, parameter :: io_maxPlotFaceVars = 12 !max face vars to output - NOT YET IMPLEMENTED


  integer, save :: io_plotVar(io_maxPlotVars), io_nPlotVars
  integer, save :: io_plotGridVar(io_maxPlotGridVars), io_nPlotGridVars

  logical, save :: io_unkActive(UNK_VARS_BEGIN:UNK_VARS_END)
  integer, save :: io_unkToGlobal(UNK_VARS_BEGIN:UNK_VARS_END) ! stores for each unk the corresponding index into io_unklabelsGlobal, or 0 if unk is inactive
  integer, save :: io_unkNonRep(UNK_VARS_BEGIN:UNK_VARS_END) ! the nonrep array each unk belongs to, or 0 if not in a nonrep
  integer, save :: io_unkNonRepIdx(UNK_VARS_BEGIN:UNK_VARS_END) ! the index of this unk in its nonrep array, or 0 if not in a nonrep
  
  ! create a temporary array to hold the 4 character variable names
  character (len=4), save :: io_unklabels(UNK_VARS_BEGIN:UNK_VARS_END) ! the global name of each unk.  corresponds witih the indices of unk.  since a nonrep array wont necessarily use all unks in each process, some of these entries may be "err"
  character (len=4), allocatable, save :: io_unklabelsGlobal(:) ! list of global unks, in no particular order.  this will be the same on all processes.
  character (len=4), save :: io_scratchGridVarlabels(SCRATCH_GRID_VARS_BEGIN:SCRATCH_GRID_VARS_END)
  character (len=4), save :: io_faceXVarLabels(NFACE_VARS)
  character (len=4), save :: io_faceYVarLabels(NFACE_VARS)
  character (len=4), save :: io_faceZVarLabels(NFACE_VARS)
  character (len=4), save :: io_plotVarStr(io_maxPlotVars)
  character (len=4), save :: io_plotGridVarStr(io_maxPlotGridVars)
  character (len=4), save :: io_plotFaceVarStr(io_maxPlotFaceVars) ! NOT YET IMPLEMENTED


  logical, save :: io_chkGuardCellsInput, io_chkGuardCellsOutput !currently only implemented in hdf5 parallel for paramesh

  
  character (len=MAX_STRING_LENGTH), save                  :: io_flashRelease = ' '
  character (len=MAX_STRING_LENGTH), save  :: io_buildDate, io_buildDir, io_buildMachine, io_fileCreationTime
  character (len=MAX_STRING_LENGTH), save  :: io_setupTimeStamp, io_buildTimeStamp
  character (len=400), save                :: io_setupCall, io_cflags, io_fflags

  !single precision or double precision

  logical, save :: io_doublePrecision

#ifdef FIXEDBLOCKSIZE    
       
   logical, save :: io_fixedBlockSize = .TRUE.   
       
#ifdef FL_NON_PERMANENT_GUARDCELLS

   integer, save :: io_ilo = GRID_ILO - NGUARD
   integer, save :: io_ihi = GRID_IHI - NGUARD
   integer, save :: io_jlo = GRID_JLO - NGUARD * K2D
   integer, save :: io_jhi = GRID_JHI - NGUARD * K2D
   integer, save :: io_klo = GRID_KLO - NGUARD * K3D
   integer, save :: io_khi = GRID_KHI - NGUARD * K3D

#else

   integer, save :: io_ilo = GRID_ILO
   integer, save :: io_ihi = GRID_IHI
   integer, save :: io_jlo = GRID_JLO
   integer, save :: io_jhi = GRID_JHI
   integer, save :: io_klo = GRID_KLO
   integer, save :: io_khi = GRID_KHI

#endif
 
#else
   logical, save :: io_fixedBlockSize = .FALSE.        
   integer, save :: io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi       
#endif

   integer, save :: io_iguard, io_jguard, io_kguard    
   integer, save :: io_iloGC, io_ihiGC, io_jloGC, io_jhiGC, io_kloGC, io_khiGC


   logical, save :: io_alwaysRestrictCheckpoint
   logical, save :: io_alwaysComputeUserVars
   logical, save :: io_ignoreForcedPlot

   !only used with HDF5, but needs to be made available
   logical, save :: io_useCollectiveHDF5
   
   !plotfile accuracy controls: enabled, these allow us to output plotfile data
   !in double precision both default to false.
   logical, save :: io_plotfileMetadataDP
   logical, save :: io_plotfileGridQuantityDP

   !Used in data restriction decision.
   logical, save :: io_outputInStack

   integer, save :: io_fileFormatVersion

   !This parameter ensures that all floating point data transfers (reads and
   !writes) involve the same memory and disk type, e.g. single precision in
   !both memory and file.  This prevents HDF5 library reverting to independent
   !I/O as would happen if the data transfer involves data that is double
   !precision in memory and single precision in file.
   logical, save :: io_type_matched_xfer

   ! This flag is set to true when a plot is written. It is set to
   ! false after one cycle
   logical, save :: io_wrotePlot = .false.

   ! The name of the plot file that was just written
   character(len=MAX_STRING_LENGTH), save :: io_oldPlotFileName

   integer, save :: io_rayFileId
   integer, save :: io_protonFileId

   logical, save :: io_reduceGcellFills
   logical, save :: io_summaryOutputOnly


   type tree_data_t
      real, dimension(:,:,:), pointer :: bnd_box
      real, dimension(:,:), pointer :: coord
      real, dimension(:,:), pointer :: bsize
      integer, dimension(:,:), pointer :: bflags
      integer, dimension(:,:), pointer :: gid
      integer, dimension(:), pointer :: nodetype
      integer, dimension(:), pointer :: lrefine
      integer, dimension(:), pointer :: which_child
      integer, dimension(:), pointer :: procnumber
      integer, dimension(:,:,:,:,:), pointer :: gsurr_blks
   end type tree_data_t

end module IO_data
