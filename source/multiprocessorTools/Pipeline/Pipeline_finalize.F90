!!****f* source/multiprocessorTools/Pipeline/Pipeline_finalize
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

  implicit none

  return
end subroutine Pipeline_finalize
