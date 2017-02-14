!!****if* source/Simulation/SimulationMain/unitTest/Pipeline/Linear/pl_localPipelineSetup
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
!!  Overriding version of the original. This routine sets up the local linear pipeline
!!  structure. It determines the number of channels and the channel processor list on the
!!  current processor.
!!
!!  The linear setup is simple. For n processors we have for their ranks:
!!
!!                      0 -> 1 -> 2 -> 3 -> 4 -> 5 -> ... -> n-1
!!
!!  Procedure:
!!
!!  Given a rank p, we need to establish all channels belonging to that particular
!!  rank. We have:
!!
!!                       channel 1  ->  p - 1  (if p > 0)
!!                       channel 2  ->  p + 1  (if p < n-1)
!!
!!  The sending channel is always the channel with the highest rank.
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

  use Pipeline_data,     ONLY : pl_doLog,       &
                                pl_logUnit,     &
                                pl_numChannels, &
                                pl_procList,    &
                                pl_rank,        &
                                pl_size

  use Driver_interface,  ONLY : Driver_abortFlash

  use Simulation_data,   ONLY : sim_depthPascalTriangle, &
                                sim_globalMe,            &
                                sim_globalNumProcs,      &
                                sim_numSendingChannels,  &
                                sim_sendingChannels

  implicit none

  integer :: n,p
  integer :: proc

  integer :: channel (1:2)
!
!
!     ...Set up the simulation pipeline structure.
!
!
  if (sim_globalMe /= pl_rank) then
      call Driver_abortFlash ('[pl_localPipelineSetup] ERROR: Simulation/pipeline rank mismatch!')
  end if

  if (sim_globalNumProcs /= pl_size) then
      call Driver_abortFlash ('[pl_localPipelineSetup] ERROR: Simulation/pipeline size mismatch!')
  end if

  pl_numChannels = 0
  sim_numSendingChannels = 0

  p = pl_rank

  if (p > 0) then
      pl_numChannels = pl_numChannels + 1
      channel (pl_numChannels) = p - 1
  end if

  if (p < pl_size - 1) then
      pl_numChannels = pl_numChannels + 1
      sim_numSendingChannels = sim_numSendingChannels + 1
      channel (pl_numChannels) = p + 1
      sim_sendingChannels (sim_numSendingChannels) = p + 1
  end if

  allocate (pl_procList (1:pl_numChannels))

  pl_procList (1:pl_numChannels) = channel (1:pl_numChannels)
!
!
!     ...Print the local pipeline structure to the log file (if requested).
!
!
  if (pl_doLog) then

      write (pl_logUnit,'(a)'   ) ' ----------------------------------------'
      write (pl_logUnit,'(a,i6)') '          Processor ID = ',pl_rank 
      write (pl_logUnit,'(a,i4)') '    Number of channels = ',pl_numChannels 

      if (pl_numChannels > 0) then
          do n = 1,pl_numChannels
             proc = pl_procList (n)
             write (pl_logUnit,'(a,i4,a,i6)'   ) '        Channel ',n,' -> Processor ',proc 
          end do
      end if

      write (pl_logUnit,'(a)'   ) ' ----------------------------------------'

  end if
!
!
!    ...Ready!
!
!
  return
end subroutine pl_localPipelineSetup
