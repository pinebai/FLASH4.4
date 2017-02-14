!!****if* source/flashUtilities/Pipeline/PipelineMain/Pipeline_localCreate
!!
!! NAME
!!  
!!  Pipeline_localCreate
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localCreate (integer,           intent (in)           :: itemSize,
!!                             integer,           intent (in)           :: maxItems,
!!                             integer,           intent (in)           :: channelSize,
!!                             character (len=*), intent (in), optional :: logName)
!!
!! DESCRIPTION
!!
!!  Creates a local pipeline section on a processor inside a FLASH application. A (global)
!!  pipeline is defined as a set of processors connected by a closed system of channels
!!  through which items will 'flow' from one processor to another. Feeding of the pipeline
!!  with items is done by a special feeding API, which is called from the unit where the
!!  pipeline is being used for. As a default, the pipeline structure is defined from the
!!  current structure of the grid, but can also be defined in an abstract way, as long as
!!  the channels between each of the processors are set up in a consistent (closed) way.
!!  This means that no open channels can exist, for which sending takes place but no
!!  corresponding receive is ever posted, a situation which would lead to deadlocks in
!!  communication.
!!
!!  An optional 'logName' can be passed, indicating that the user wants to know details
!!  details about the workings of the pipeline for this particular application. The log
!!  file created is specific for each processor, i.e. each processor writes only to its
!!  own log file with name: <basenm>Pipeline<logName>_<procID>.log , where <basenm> is the
!!  base name of the current application and <procID> is the processor ID of the current
!!  processor in 000x format.
!!
!!  The only arguments the user has to provide is: 1) the size (number of elements) of each
!!  item, 2) the maximum number of such items and 3) the size of each pipeline channel, i.e.
!!  the number of items each channel can maximally hold, before sending them to the
!!  connecting processor. Details about the pipeline structure to be set up is deferred to
!!  an internal routine 'pl_localPipelineSetup', thus giving the user the possibility
!!  for setting up his personal pipeline (i.e. not based on the default -> grid structure)
!!  in case he needs it.
!!
!! ARGUMENTS
!!
!!  itemSize    : size of each item (number of elements per item)
!!  maxItems    : maximum number of items
!!  channelSize : channel size (number of items each channel can hold at a time)
!!  logName     : optional log file name
!!
!! NOTES
!!
!!  1) This operation is local. No other processors are involved.
!!
!!  2) Channels are only defined as processor ID's that are not equal to the current processor
!!     ID. If no such channels are defined, the number of channels is zero.
!!
!!  3) Pipelines can only be used sequentially (one after the other). No simultaneous pipelines
!!     can be in operation right now.
!!
!!***

subroutine Pipeline_localCreate (itemSize, maxItems, channelSize, logName)

  use Pipeline_data,     ONLY : pl_baseName,           &
                                pl_channelSize,        &
                                pl_comm,               &
                                pl_commGlobal,         &
                                pl_commMesh,           &
                                pl_doLog,              &
                                pl_itemBuf,            &
                                pl_itemCount,          &
                                pl_itemSize,           &
                                pl_logName,            &
                                pl_logUnit,            &
                                pl_maxItems,           &
                                pl_numChannels,        &
                                pl_pipelineCreated,    &
                                pl_procStatusGlobal,   &
                                pl_procStatusLocal,    &
                                pl_rank,               &
                                pl_recvBuf,            &
                                pl_recvCount,          &
                                pl_recvIndex,          &
                                pl_recvRequest,        &
                                pl_recvStatus,         &
                                pl_sendBuf,            &
                                pl_sendCount,          &
                                pl_sendIndex,          &
                                pl_sendRequest,        &
                                pl_sendStatus,         &
                                pl_size

  use Driver_interface,  ONLY : Driver_abortFlash,       &
                                Driver_checkMPIErrorCode

  use pl_interface,      ONLY : pl_localPipelineSetup

  implicit none

#include "Flash.h"
#include "Pipeline.h"
#include "constants.h"
 include "Flash_mpi.h"

  character (len=*), intent (in), optional :: logName
  integer,           intent (in)           :: itemSize
  integer,           intent (in)           :: maxItems
  integer,           intent (in)           :: channelSize

  character (len=4) :: charRank
  integer           :: error
  integer           :: ut_getFreeFileUnit
!
!
!     ...Check first, if any previous pipeline has been created. If the case, we
!        need to abort, since we cannot have two or more pipelines created at the
!        same time.
!
!    
  if (pl_pipelineCreated) then
      call Driver_abortFlash ('[Pipeline_localCreate] ERROR: Only 1 pipeline at a time possible!')
  end if
!
!
!     ...Copy the arguments to the internal variables and prepare the log file (if needed).
!        Set also the internal communicator.
!
!
  pl_itemSize    = itemSize
  pl_maxItems    = maxItems
  pl_channelSize = channelSize
  pl_doLog       = .false.
  pl_comm        = pl_commGlobal      ! based on global communicator for now

  call MPI_Comm_rank            (pl_comm, pl_rank, error)
  call Driver_checkMPIErrorCode (error)
  call MPI_Comm_size            (pl_comm, pl_size, error)
  call Driver_checkMPIErrorCode (error)

  if (present (logName)) then
      write (charRank,'(i4.4)') pl_rank
      pl_logName = trim (pl_baseName) // "Pipeline" // trim (adjustl (logName)) // "_" // charRank // ".log"
      pl_logUnit = ut_getFreeFileUnit ()
      pl_doLog = (pl_logUnit >= 0)
      if (pl_doLog) open (unit = pl_logUnit, file = pl_logName)
  end if
!
!
!     ...Set up the structure of the pipeline. After calling this routine, the number of
!        channels and their processor lists are ready. Allocate the necessary internal arrays.
!
!
  call pl_localPipelineSetup ()

  if (pl_numChannels > 0) then

      allocate (pl_sendRequest     (1:pl_numChannels))
      allocate (pl_sendIndex       (1:pl_numChannels))
      allocate (pl_sendCount       (1:pl_numChannels))
      allocate (pl_recvRequest     (1:pl_numChannels))
      allocate (pl_recvIndex       (1:pl_numChannels))
      allocate (pl_recvCount       (1:pl_numChannels))
      allocate (pl_itemBuf         (1:pl_itemSize , 1:pl_maxItems))
      allocate (pl_sendStatus      (1:MPI_STATUS_SIZE , 1:pl_numChannels))
      allocate (pl_recvStatus      (1:MPI_STATUS_SIZE , 1:pl_numChannels))
      allocate (pl_sendBuf         (1:pl_itemSize , 1:pl_channelSize , 1:pl_numChannels))
      allocate (pl_recvBuf         (1:pl_itemSize , 1:pl_channelSize , 1:pl_numChannels))

      pl_sendRequest (:) = MPI_REQUEST_NULL
      pl_recvRequest (:) = MPI_REQUEST_NULL

  end if
!
!
!     ...Set up the array which will monitor global status of the pipeline at any stage.
!        The global status array will contain communication info (how many sends/receive pairs
!        are still active on each processor) as well as receive/item buffer info (empty or
!        non-empty). A send/receive pair is considered inactive, if a receive operation has
!        completed (the receive buffer can be used). The global status array has three entries
!        for each processor: PL_STATUS_COMM (communication status), PL_STATUS_RECV (receive
!        buffer status) and PL_STATUS_ITEM (item buffer status). They are explained in detail
!        in the include file Pipeline.h. The local status array is used to accumulate local
!        processor info, which will then be reduced to global info.
!
!        Note that even if all pipeline communication has ceased and all receive/item buffers
!        have been processed, the pipeline can still be used by adding new items from outside
!        to the pipeline. It is the application users responsibility to decide when to shut
!        (deactivate + destroy) the pipeline down.
!
!
  allocate (pl_procStatusGlobal (0:pl_size-1, 1:PL_STATUS_SIZE))
  allocate (pl_procStatusLocal  (0:pl_size-1, 1:PL_STATUS_SIZE))

  pl_procStatusGlobal (:,PL_STATUS_COMM) = 0
  pl_procStatusGlobal (:,PL_STATUS_RECV) = 0
  pl_procStatusGlobal (:,PL_STATUS_ITEM) = 0

  pl_procStatusLocal  (:,PL_STATUS_COMM) = 0
  pl_procStatusLocal  (:,PL_STATUS_RECV) = 0
  pl_procStatusLocal  (:,PL_STATUS_ITEM) = 0
!
!
!     ...The pipeline is considered created now, even if no channels are present.
!
!    
  pl_pipelineCreated = .true.

  if (pl_doLog) then
      write (pl_logUnit,'(a,i6)') ' Pipeline created on proc ID ', pl_rank
  end if
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_localCreate
