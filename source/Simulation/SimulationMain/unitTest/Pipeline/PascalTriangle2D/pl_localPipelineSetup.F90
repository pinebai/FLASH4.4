!!****if* source/Simulation/SimulationMain/unitTest/Pipeline/PascalTriangle2D/pl_localPipelineSetup
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
!!  Overriding version of the original. This routine sets up the local pipeline structure.
!!  It determines the number of channels and the channel processor list on the current
!!  processor.
!!
!!  For a Pascal triangle of depth d, the pipeline will be set up the following way (Pascal
!!  triangle arrangement of ranks):
!!
!!                                         0              d = 0
!!                                        / \
!!                                       1   2            d = 1
!!                                      / \ / \
!!                                     3   4   5          d = 2
!!                                    / \ / \ / \
!!                                   6   7   8   9        d = 3
!!                                  / \ / \ / \ / \
!!                                10  11  12  13  14      d = 4
!!
!!  Procedure:
!!
!!  Given a rank p, we need to establish 2 things: 1) to what depth level in the Pascal
!!  triangle does p belong to and 2) determine all channels belonging to that particular
!!  rank. The depth level d is found from the following equation:
!!
!!                       d = int ( (-3+sqrt(8p+1))/2 ) + 1
!!
!!  The possible 4 channels for a particular p are:
!!
!!                       channel 1  ->  p - d - 1  (only if it exists at d - 1 depth)
!!                       channel 2  ->  p - d      (only if it exists at d - 1 depth)
!!                       channel 3  ->  p + d + 1  (only if it exists at d + 1 depth)
!!                       channel 4  ->  p + d + 2  (only if it exists at d + 1 depth)
!!
!!  At a particular depth d, the minimum/maximum allowed ranks are:
!!
!!                                p(min) = d * (d + 1) / 2
!!                                p(max) = p(min) + d
!!
!!  This allows us to check the existence of the channels as follows:
!!
!!                 channel 1  ->  p - d - 1  (only if >= d*(d-1)/2)
!!                 channel 2  ->  p - d      (only if <= (d-1)*(d+2)/2)
!!                 channel 3  ->  p + d + 1  (only if >= (d+1)*(d+2)/2 and d+1 <= total depth)
!!                 channel 4  ->  p + d + 2  (only if <= (d+1)*(d+4)/2 and d+1 <= total depth)
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

  integer :: d,n,p
  integer :: dmin, dmax
  integer :: pmaxDown, pmaxUp
  integer :: pminDown, pminUp
  integer :: proc

  integer :: channel (1:4)
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
  d = int ((-3.0 + sqrt (8.0 * p + 1.0)) / 2) + 1

  dmin = 0
  dmax = sim_depthPascalTriangle

  if (d <= dmax) then

      pminDown = (d - 1) *  d      / 2
      pminUp   = (d + 1) * (d + 2) / 2
      pmaxDown = pminDown + d - 1
      pmaxUp   = pminUp   + d + 1

      proc = p - d - 1
      if (proc >= pminDown .and. d - 1 >= dmin) then
          pl_numChannels = pl_numChannels + 1
          channel (pl_numChannels) = proc
      end if

      proc = p - d
      if (proc <= pmaxDown .and. d - 1 >= dmin) then
          pl_numChannels = pl_numChannels + 1
          channel (pl_numChannels) = proc
      end if

      proc = p + d + 1
      if (proc >= pminUp .and. d + 1 <= dmax) then
          pl_numChannels = pl_numChannels + 1
          sim_numSendingChannels = sim_numSendingChannels + 1
          channel (pl_numChannels) = proc
          sim_sendingChannels (sim_numSendingChannels) = proc
      end if

      proc = p + d + 2
      if (proc <= pmaxUp .and. d + 1 <= dmax) then
          pl_numChannels = pl_numChannels + 1
          sim_numSendingChannels = sim_numSendingChannels + 1
          channel (pl_numChannels) = proc
          sim_sendingChannels (sim_numSendingChannels) = proc
      end if

      allocate (pl_procList (1:pl_numChannels))

      pl_procList (1:pl_numChannels) = channel (1:pl_numChannels)

  end if
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
