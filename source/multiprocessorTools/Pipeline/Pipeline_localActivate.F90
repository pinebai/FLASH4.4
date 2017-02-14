!!****f* source/multiprocessorTools/Pipeline/Pipeline_localActivate
!!
!! NAME
!!  
!!  Pipeline_localActivate
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localActivate ()
!!
!! DESCRIPTION
!!
!!  Activates locally the pipeline created. Activation consists in posting speculative
!!  receive messages, awaiting for items to be received on the local processor.
!!  If no channels were created, no communication will happen on the local processor,
!!  but the processor is still considered to be active inside the pipeline.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  It is always a good idea to check beforehand, if the created global pipeline
!!  structure is ok by calling the routine 'Pipeline_globalCheckStructure'.
!!  Otherwise sends might be posted on channels which have no receive channel
!!  counterpart on another processor.
!!
!!***

subroutine Pipeline_localActivate ()

  implicit none

  return
end subroutine Pipeline_localActivate
