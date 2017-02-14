!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Async/ed_commProgressPipeline
!!
!!  NAME     
!!   ed_commProgressPipeline
!!
!!  SYNOPSIS
!!   ed_commProgressPipeline(logical, intent(IN) :: doFlush)
!!
!!  DESCRIPTION 
!!
!!  ARGUMENTS
!!
!!  NOTES
!!
!!  Global progress is only possible if *all* MPI ranks call this subroutine.
!!
!!  SIDE EFFECTS
!!
!!***

subroutine ed_commProgressPipeline(doFlush)
  use ed_commData, ONLY : ed_commPipelineState, ed_commCounterState, &
       PIPELINE_ACTIVE, PIPELINE_CLOSING, PIPELINE_INACTIVE, &
       COUNTER_INACTIVE
  use ed_commInterface, ONLY : ed_commCheckRayExchange
  use UTPipeline, ONLY : UTPipeline_progressComm, UTPipeline_closeSendChannels, &
       UTPipeline_isDone, UTPipeline_finalize
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  logical, intent(IN) :: doFlush
  logical :: isClosing, isPipelineDone
  
  select case (ed_commPipelineState)
  case (PIPELINE_ACTIVE)
     call UTPipeline_progressComm(doFlush)
     if (COUNTER_INACTIVE == ed_commCounterState) then
        !Attempt a clean shut-down of the pipeline if the counter target is
        !met.  It may take more than one attempt due to the non-blocking
        !design of the underlying pipeline module.
        call UTPipeline_closeSendChannels(isClosing)
        if (isClosing) ed_commPipelineState = PIPELINE_CLOSING
     end if

  case (PIPELINE_CLOSING)
     call UTPipeline_isDone(isPipelineDone)
     if (isPipelineDone) then
        ed_commPipelineState = PIPELINE_INACTIVE
        call ed_commCheckRayExchange()

        call UTPipeline_finalize() !tmp
     end if

  case (PIPELINE_INACTIVE)
     continue !no-op

  case default
     call Driver_abortFlash("Unknown pipeline state")
  end select

end subroutine ed_commProgressPipeline
