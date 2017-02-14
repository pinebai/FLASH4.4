!!****f* source/multiprocessorTools/Pipeline/Pipeline_localProgress
!!
!! NAME
!!  
!!  Pipeline_localProgress
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localProgress ()
!!
!! DESCRIPTION
!!
!!  This routine can be considered the motor of the pipeline. It makes sure the items
!!  are transported through the channels on the local processor. As such the routine
!!  should be called often to ensure best possible communication progress between the
!!  channels. Progress on a processor can stall, if the item buffer is full.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  All MPI ranks must call this subroutine to avoid possible deadlock. This subroutine
!!  must remain non-blocking.
!!
!!***

subroutine Pipeline_localProgress ()

  implicit none

  return
end subroutine Pipeline_localProgress
