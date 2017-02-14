!!****if* source/multiprocessorTools/Pipeline/PipelineMain/pl_masterPipelineMap
!!
!! NAME
!!  
!!  pl_masterPipelineMap
!!
!! SYNOPSIS
!! 
!!  call pl_masterPipelineMap ()
!!
!! DESCRIPTION
!!
!!  This routine sets up the pipeline map at the master rank. The pipeline map consists
!!  of a two-dimensional integer array, which contains info about how many channels (#)
!!  each rank has and also all the processor ID's of these channels:
!!
!!               Rank| # ||    channel proc ID's
!!               -----------------------------------
!!                 0 | 4 || 11 |  5 |  6 |  7 |
!!                 1 | 1 || 13 |
!!                 2 | 3 ||  1 |  8 |  5 |  4 |
!!                 3 | 5 || 17 | 15 | 16 |  6 |  4 |
!!                  .............................
!!
!!  For convenience, the pipeline map is declared such that the ranks are located in the
!!  columns and the channels in the rows of the pipeline map array.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  After exiting this routine, the maximum number of channels per processor has been
!!  set on all processors and the pipeline map has been allocated and created on the master.
!!
!!***

subroutine pl_masterPipelineMap ()

  use Pipeline_data,     ONLY : pl_comm,        &
                                pl_maxChannels, &
                                pl_masterRank,  &
                                pl_numChannels, &
                                pl_pipelineMap, &
                                pl_procList,    &
                                pl_rank,        &
                                pl_size

  use pl_interface,      ONLY : pl_printPipelineMap

  implicit none

  include "Flash_mpi.h"

  integer :: error

  integer, parameter :: notPresent = -1   ! must be negative

  integer, allocatable :: procPipelineMap (:)
!
!
!     ...Determine first on all processors the maximum number of channels per processor.
!
!
  call MPI_Allreduce (pl_numChannels, &
                      pl_maxChannels, &
                      1,              &
                      FLASH_INTEGER,  &
                      MPI_Max,        &
                      pl_comm,        &
                      error           )
!
!
!     ...Allocate and create necessary arrays.
!
!
  allocate (procPipelineMap (1:1+pl_maxChannels))

  procPipelineMap (1) = pl_numChannels

  if (pl_numChannels > 0) then
      procPipelineMap (2:1+pl_numChannels) = pl_procList (1:pl_numChannels)
  end if

  if (pl_maxChannels > pl_numChannels) then
      procPipelineMap (2+pl_numChannels:1+pl_maxChannels) = notPresent
  end if
!
!
!     ...We allocate the pipeline map on all processors, even if only the master
!        processor is going to use it. If this is not done, older versions of intel
!        compilers in debug mode will complain with an error message 'attempt to
!        fetch from allocatable variable when it is not allocated' and the code
!        crashes.
!
!
!  if (pl_rank == pl_masterRank) then
      allocate (pl_pipelineMap (1:1+pl_maxChannels, 0:pl_size-1))
!  end if
!
!
!     ...Create the pipeline map on the master.
!
!
  call MPI_Gather (procPipelineMap (1),        &
                   1+pl_maxChannels,           &
                   FLASH_INTEGER,              &
                   pl_pipelineMap (1,pl_rank), &
                   1+pl_maxChannels,           &
                   FLASH_INTEGER,              &
                   pl_masterRank,              &
                   pl_comm,                    &
                   error                       )

  deallocate (procPipelineMap)
!
!
!     ...Print the pipeline map to the log file (if requested).
!
!
  call pl_printPipelineMap ()
!
!
!    ...Ready!
!
!
  return
end subroutine pl_masterPipelineMap
