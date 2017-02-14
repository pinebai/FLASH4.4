!!****if* source/multiprocessorTools/Pipeline/PipelineMain/Pipeline_localDestroy
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

  use Pipeline_data, ONLY : pl_doLog,              &
                            pl_itemBuf,            &
                            pl_logUnit,            &
                            pl_pipelineActive,     &
                            pl_pipelineCreated,    &
                            pl_pipelineMap,        &
                            pl_procList,           &
                            pl_rank,               &
                            pl_recvBuf,            &
                            pl_recvCount,          &
                            pl_recvIndex,          &
                            pl_recvRequest,        &
                            pl_recvStatus,         &
                            pl_sendBuf,            &
                            pl_sendCount,          &
                            pl_sendIndex,          &
                            pl_sendRequest,        &
                            pl_sendStatus

  implicit none
!
!
!     ...If pipeline is still active, we cannot destroy it.
!
!
  if (pl_pipelineActive) then
      call Driver_abortFlash ('[Pipeline_localDestroy] ERROR: Cannot destroy an active pipeline!')
  end if
!
!
!     ...Deallocate any pipeline arrays (if any).
!
!
  if (allocated (pl_itemBuf)           ) deallocate (pl_itemBuf)
  if (allocated (pl_sendBuf)           ) deallocate (pl_sendBuf)
  if (allocated (pl_sendStatus)        ) deallocate (pl_sendStatus)
  if (allocated (pl_sendRequest)       ) deallocate (pl_sendRequest)
  if (allocated (pl_sendIndex)         ) deallocate (pl_sendIndex)
  if (allocated (pl_sendCount)         ) deallocate (pl_sendCount)
  if (allocated (pl_recvBuf)           ) deallocate (pl_recvBuf)
  if (allocated (pl_recvStatus)        ) deallocate (pl_recvStatus)
  if (allocated (pl_recvRequest)       ) deallocate (pl_recvRequest)
  if (allocated (pl_recvIndex)         ) deallocate (pl_recvIndex)
  if (allocated (pl_recvCount)         ) deallocate (pl_recvCount)
  if (allocated (pl_procList)          ) deallocate (pl_procList)
  if (allocated (pl_pipelineMap)       ) deallocate (pl_pipelineMap)
!
!
!     ...Close the log file (if present) after recording the termination message.
!
!
  if (pl_doLog) then
      write (pl_logUnit,'(a,i6)') ' Pipeline destroyed       on proc ID ', pl_rank
      close (pl_logUnit)
  end if
!
!
!     ...The pipeline is considered now destroyed.
!
!    
  pl_pipelineCreated = .false.
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_localDestroy
