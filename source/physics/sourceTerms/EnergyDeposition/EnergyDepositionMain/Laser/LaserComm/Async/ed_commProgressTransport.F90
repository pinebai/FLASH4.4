!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Async/ed_commProgressTransport
!!
!!  NAME     
!!   ed_commProgressTransport
!!
!!  SYNOPSIS
!!   call ed_commProgressTransport(logical, optional, intent(IN) :: optionalForceProgress)
!!
!!  DESCRIPTION
!!   This subroutine progresses communication by testing completion of
!!   non-blocking communication calls in the shared counter and
!!   pipeline exchange modules.  It must be called frequently to
!!   ensure best possible communication progress:
!!
!!    1). All MPI ranks test receive buffers for newly received rays.
!!        This is the only place where MPI ranks with ed_rayCount = 0
!!        progress global communication!
!!    2). MPI ranks with rays in send buffers are forced to flush
!!        the fullest send buffer.  This is absolutely necessary
!!        when processing the last few rays.
!!
!!  ARGUMENTS
!!   optionalForceProgress: Forces communication progress.
!!
!!  NOTES
!!   All MPI ranks must call this subroutine to avoid possible
!!   deadlock.
!!
!!   This subroutine must remain non-blocking.  The only way for the
!!   pipeline to shutdown is for there to be no more rays in internal
!!   buffers inc. send promises.  This may mean ed_commUpdateRays still
!!   has useful work to do.
!!
!!  SIDE EFFECTS
!!   ed_commRaysUntilMsgTest
!!
!!***

subroutine ed_commProgressTransport(optionalForceProgress)
  use ed_commData, ONLY : ed_commRaysUntilMsgTest, ed_commRaysBetweenMsgTest
  !$ use ed_commData, ONLY : ed_commLock
  use ed_commInterface, ONLY : ed_commProgressPipeline, &
       ed_commProgressCounter
  implicit none
  logical, optional, intent(IN) :: optionalForceProgress
  logical :: doProgress, doFlush, forceProgress

  if (present(optionalForceProgress)) then
     forceProgress = optionalForceProgress
  else
     forceProgress = .false.
  end if

  if (forceProgress) then
     !Periodic flushing prevents deadlock when processing the last few rays!
     doProgress = .true.
     doFlush = .true.
  else
     !$omp critical (ShouldCommHappen)
     ed_commRaysUntilMsgTest = ed_commRaysUntilMsgTest - 1
     if (ed_commRaysUntilMsgTest == 0) then
        ed_commRaysUntilMsgTest = ed_commRaysBetweenMsgTest
        doProgress = .true.
     else
        doProgress = .false.
     end if
     !$omp end critical (ShouldCommHappen)
     doFlush = .false.
  end if

  if (doProgress) then
     !$ call omp_set_lock(ed_commLock)
     call ed_commProgressCounter()
     call ed_commProgressPipeline(doFlush)
     !$ call omp_unset_lock(ed_commLock)
  end if

end subroutine ed_commProgressTransport
