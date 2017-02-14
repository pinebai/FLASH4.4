!!****f* source/multiprocessorTools/Pipeline/Pipeline_localDestroy
!!
!! NAME
!!  
!!  Pipeline_localDestroy
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localDestroy ()
!!
!! DESCRIPTION
!!
!!  Destroys the local pipeline structure. Destruction means that all arrays
!!  (the plumbing of the pipeline) are deallocated. Also the local log file
!!  will be closed (if present). 
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  This routine should onle be called after the local pipeline has been
!!  deactivated.
!!
!!***

subroutine Pipeline_localDestroy ()

  implicit none

  return
end subroutine Pipeline_localDestroy
