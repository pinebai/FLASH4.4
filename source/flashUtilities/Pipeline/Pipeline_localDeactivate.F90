!!****f* source/flashUtilities/Pipeline/Pipeline_localDeactivate
!!
!! NAME
!!  
!!  Pipeline_localDeactivate
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localDeactivate (logical, optional, intent (in) :: doAsyncReturn)
!!
!! DESCRIPTION
!!
!!  Deactivates the local pipeline environment. Deactivation means smooth shutdown
!!  of the local pipeline communitation environment. There are two ways the local
!!  pipeline communicator can shut down: 1) all processors return at will or 2) all
!!  processors return at the same time.
!!
!! ARGUMENTS
!!
!!  doAsyncReturn : if set true, each processor returns at will
!!
!! NOTES
!!
!!  This routine should be called at the end of a pipeline application. Deactivation
!!  is essential before destroying the pipeline. One should never destruct a pipeline
!!  before deactivating it.
!!
!!***

subroutine Pipeline_localDeactivate (doAsyncReturn)

  implicit none

  logical, optional, intent (in) :: doAsyncReturn

  return
end subroutine Pipeline_localDeactivate
