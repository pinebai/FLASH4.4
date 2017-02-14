!!****if* source/multiprocessorTools/Pipeline/PipelineMain/Pipeline_localActivate
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

  use Pipeline_data,     ONLY : pl_doLog,           &
                                pl_itemCount,       &
                                pl_logUnit,         &
                                pl_numChannels,     &
                                pl_pipelineActive,  &
                                pl_pipelineCreated, &
                                pl_procItemBalance, &
                                pl_rank,            &
                                pl_recvCount,       &
                                pl_sendCount

  use Driver_interface,  ONLY : Driver_abortFlash

  use pl_interface,      ONLY : pl_localPostRecvMsg

  implicit none

#include "Flash.h"
#include "constants.h"
 include "Flash_mpi.h"

  integer :: channel
!
!
!     ...Check first, if any pipeline has been created. If not the case, we abort.
!
!    
  if (.not. pl_pipelineCreated) then
      call Driver_abortFlash ('[Pipeline_localActivate] ERROR: No pipeline created!')
  end if
!
!
!     ...Check next, if any previous pipeline is still active. If the case, we need
!        to abort, since we cannot have two or more pipelines active at the same time.
!        Otherwise set the pipeline inactive for the moment. The pipeline is considered
!        active, once the first initial communication is activated, i.e. the initial
!        receives have been posted.
!
!    
  if (pl_pipelineActive) then
      call Driver_abortFlash ('[Pipeline_localActivate] ERROR: Only 1 active pipeline possible!')
  end if
!
!
!     ...Activate (locally) the pipeline:
!
!         1) Sets all sending/receiving counts to 0
!         2) Posts all receiving messages on all channels
!            (waiting for the sending messages to arrive)
!
!    
  pl_itemCount = 0
  pl_procItemBalance = 0

  if (pl_numChannels > 0) then

      pl_sendCount (1:pl_numChannels) = 0
      pl_recvCount (1:pl_numChannels) = 0

      do channel = 1, pl_numChannels
         call pl_localPostRecvMsg (channel)
      end do

  end if
!
!
!     ...The pipeline is considered active now.
!
!    
  pl_pipelineActive = .true.

  if (pl_doLog) then
      write (pl_logUnit,'(a,i6)') ' Pipeline is active       on proc ID ', pl_rank
  end if
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_localActivate
