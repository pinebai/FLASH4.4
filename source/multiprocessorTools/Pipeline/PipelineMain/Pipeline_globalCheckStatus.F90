!!****if* source/multiprocessorTools/Pipeline/PipelineMain/Pipeline_globalCheckStatus
!!
!! NAME
!!  
!!  Pipeline_globalCheckStatus
!!
!! SYNOPSIS
!! 
!!  call Pipeline_globalCheckStatus (logical, intent (out) :: empty)
!!
!! DESCRIPTION
!!
!!  Returns the global status of the pipeline. This is a blocking operation that
!!  currently determines, if the pipeline is considered empty, i.e. there are no
!!  more items in the pipeline.
!!
!!  Procedure:
!!
!!  The item balance counts (number of items added vs. number of items removed from
!!  pipeline) on each processor are all summation reduced and the resulting total
!!  pipeline item count redistributed on all processors. Since the number of total
!!  items in the pipeline cannot be < 0, a negative number indicates a bug someplace.
!!  Note that the individual processor item balance counts can be < 0. This happens
!!  if on a processor no items were fed into the pipeline, but many items were retrieved
!!  on that same processor.
!!
!! ARGUMENTS
!!
!!  empty : is set true, if no more items are in pipeline
!!
!! NOTES
!!
!!  This is a collective operation, hence all processors must call it.
!!
!!***

subroutine Pipeline_globalCheckStatus (empty)

  use Pipeline_data, ONLY : pl_comm,              &
                            pl_doLog,             &
                            pl_logUnit,           &
                            pl_pipelineItemCount, &
                            pl_procItemBalance,   &
                            pl_size

  implicit none

  include "Flash_mpi.h"

  logical, intent (out) :: empty

  integer :: error
!
!
!     ...Check total item count in pipeline. If > 0, pipeline is not finished.
!
!
  call MPI_Allreduce (pl_procItemBalance,       &
                      pl_pipelineItemCount,     &
                      1,                        &
                      FLASH_INTEGER,            &
                      MPI_Sum,                  &
                      pl_comm,                  &
                      error                     )

  if (pl_doLog) then
      write (pl_logUnit,'(a,i10)') ' Processor Item Balance =        ', pl_procItemBalance
      write (pl_logUnit,'(a,i10)') ' Pipeline  Item Count   =        ', pl_pipelineItemCount
  end if

  empty = pl_pipelineItemCount == 0
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_globalCheckStatus
