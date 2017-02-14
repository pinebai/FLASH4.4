!!****if* source/multiprocessorTools/Pipeline/PipelineMain/Pipeline_finalize
!!
!! NAME
!!  
!!  Pipeline_finalize
!!
!! SYNOPSIS
!! 
!!  call Pipeline_finalize ()
!!
!! DESCRIPTION
!!
!!  Finalizes the Pipeline unit. It deactivates and destroys any pipeline
!!  that is currently in use. The deactivation is done in a synchronous way
!!  to allow all processors to finalize the Pipeline unit at the same time.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  none.
!!
!!***

subroutine Pipeline_finalize ()

  use Pipeline_interface, ONLY : Pipeline_localDeactivate, &
                                 Pipeline_localDestroy

  implicit none
!
!
!     ...Deactivate and destroy any active local pipeline.
!
!
  call Pipeline_localDeactivate (doAsyncReturn = .false.)
  call Pipeline_localDestroy    ()
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_finalize
