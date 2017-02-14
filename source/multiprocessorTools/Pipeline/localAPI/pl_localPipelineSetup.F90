!!****if* source/multiprocessorTools/Pipeline/localAPI/pl_localPipelineSetup
!!
!! NAME
!!  
!!  pl_localPipelineSetup
!!
!! SYNOPSIS
!! 
!!  call pl_localPipelineSetup ()
!!
!! DESCRIPTION
!!
!!  This routine sets up the local pipeline structure. It determines the number of
!!  channels and the channel processor list on the current processor. The default
!!  structure is based on the current grid structure and this is what is done below.
!!  If the user whishes a pipeline based on a different structure, (s)he needs to
!!  export this routine to her/his personal application (simmulation) unit.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  The pipeline structure is set once: 1) pl_numChannels and 2) pl_procList are defined.
!!  pl_procList is an array and is allocated here.
!!
!!***

subroutine pl_localPipelineSetup ()

  implicit none

  return
end subroutine pl_localPipelineSetup
