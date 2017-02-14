!!****if* source/multiprocessorTools/Pipeline/localAPI/pl_printPipelineMap
!!
!! NAME
!!
!!  pl_printPipelineMap
!!
!! SYNOPSIS
!!
!!  call pl_printPipelineMap ()
!!
!! DESCRIPTION
!!
!!  Prints the pipeline map to the master's log file (if log files are requested).
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  Only the master rank can do the printing as he is the only one who has the
!!  pipeline map.
!!
!!***

subroutine pl_printPipelineMap ()

  implicit none

  return
end subroutine pl_printPipelineMap
