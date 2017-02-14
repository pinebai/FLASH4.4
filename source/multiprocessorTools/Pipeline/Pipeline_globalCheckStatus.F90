!!****f* source/multiprocessorTools/Pipeline/Pipeline_globalCheckStatus
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

  implicit none

  logical, intent (out) :: empty

  empty = .false.

  return
end subroutine Pipeline_globalCheckStatus
