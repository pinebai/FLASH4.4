!!****if* source/multiprocessorTools/Pipeline/PipelineMain/Pipeline_data
!!
!! NAME
!!
!!  Pipeline_data
!!
!! SYNOPSIS
!!
!!  use Pipeline_data
!!  
!! DESCRIPTION
!!
!!  Data module for a Pipeline
!!  --------------------------
!!   
!!   Legend: (P) means data that is set as (P)arameters
!!           (G) means data that is (G)et from other units (driver, physical constants, etc ...)
!!           (R) means data that is supplied as input through (R)untime parameters
!!           (I) means data that is supplied as input at pipeline initialization (in argument list)
!!           (C) means data that is (C)alculated internally by the pipeline code
!!
!!   pl_baseName             (R) : The base name for the application
!!   pl_channelSize          (I) : number of items each pipeline channel can hold
!!   pl_comm                 (G) : the communicator within which the pipeline is defined
!!   pl_commGlobal           (G) : the global communicator
!!   pl_commMesh             (G) : the mesh communicator
!!   pl_doLog                (C) : is set true, if a pipeline log file is requested
!!   pl_itemBuf              (C) : total received items buffer on the current processor
!!   pl_itemCount            (C) : counts the total items received by the current processor at each stage
!!   pl_itemSize             (I) : the size (number of reals) defining each item
!!   pl_logName            (I,C) : The log name for pipeline info (each processor has its own log file)
!!   pl_logUnit              (G) : the unit number for the pipeline log file
!!   pl_masterRank           (C) : the rank that is considered the master of the pipeline
!!   pl_maxChannels          (C) : the maximum number of channels (at master rank only)
!!   pl_maxItems             (I) : the maximum number of items that the pipeline can handle
!!   pl_numChannels          (C) : the number of channels of the pipeline for the current rank 
!!   pl_pipelineActive       (C) : is the pipeline still considered active?
!!   pl_pipelineCreated      (C) : has the pipeline been created (all arrays allocated)?
!!   pl_pipelineItemCount    (C) : the number of items currently in the pipeline
!!   pl_procItemBalance      (C) : the number of items added/removed (+/-1) from pipeline for each processor
!!   pl_procList             (C) : list of processor channels that will communicate with the current rank
!!   pl_rank                 (C) : rank of current processor
!!   pl_recvBuf              (C) : items receive buffer (itemSize:itemsCount:channel)
!!   pl_recvCount            (C) : number of items received by each channel
!!   pl_recvIndex            (C) : will contain channel indices for a particular mpi receive test operation
!!   pl_recvRequest          (C) : receive request handles for each channel
!!   pl_recvStatus           (C) : the receive status of each channel
!!   pl_sendBuf              (C) : items send buffer (itemSize:itemsCount:channel)
!!   pl_sendCount            (C) : number of items to be sent by each channel
!!   pl_sendIndex            (C) : will contain channel indices for a particular mpi send test operation
!!   pl_sendRequest          (C) : send request handles for each channel
!!   pl_sendStatus           (C) : the sending status of each channel
!!   pl_size                 (C) : size (number of processors) of the pipeline communicator
!!   pl_tag                  (P) : tag identifier for sending/receiving items
!!
!!***

Module Pipeline_data

  implicit none

#include "constants.h"

  character (len = MAX_STRING_LENGTH), save :: pl_baseName
  character (len = MAX_STRING_LENGTH), save :: pl_logName

  logical, save :: pl_doLog
  logical, save :: pl_pipelineActive    = .false.
  logical, save :: pl_pipelineCreated   = .false.

  integer, save :: pl_channelSize
  integer, save :: pl_comm
  integer, save :: pl_commGlobal
  integer, save :: pl_commMesh
  integer, save :: pl_itemCount
  integer, save :: pl_itemSize
  integer, save :: pl_logUnit
  integer, save :: pl_masterRank
  integer, save :: pl_maxChannels     ! at master rank only
  integer, save :: pl_maxItems
  integer, save :: pl_numChannels
  integer, save :: pl_pipelineItemCount
  integer, save :: pl_procItemBalance
  integer, save :: pl_rank
  integer, save :: pl_size

  integer, parameter :: pl_tag = 1234

  integer, save, allocatable :: pl_procList            (:)
  integer, save, allocatable :: pl_recvCount           (:)
  integer, save, allocatable :: pl_recvIndex           (:)
  integer, save, allocatable :: pl_recvRequest         (:)
  integer, save, allocatable :: pl_sendCount           (:)
  integer, save, allocatable :: pl_sendIndex           (:)
  integer, save, allocatable :: pl_sendRequest         (:)

  real,    save, allocatable :: pl_itemBuf             (:,:)
  integer, save, allocatable :: pl_pipelineMap         (:,:)      ! at master rank only
  integer, save, allocatable :: pl_recvStatus          (:,:)
  integer, save, allocatable :: pl_sendStatus          (:,:)

  real,    save, allocatable :: pl_recvBuf             (:,:,:)    ! (itemSize:itemsCount:channel)
  real,    save, allocatable :: pl_sendBuf             (:,:,:)    ! (itemSize:itemsCount:channel)

end Module Pipeline_data
