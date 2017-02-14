!!****if* source/Simulation/SimulationMain/unitTest/Pipeline/Cubic/pl_localPipelineSetup
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
!!  For a cubic lattice with N number of nodes per side, we will label their vertices such
!!  that the numbering is consecutive along each line L for each surface S. The following
!!  shows all 4 surfaces S = 1,2,3,4 for a cube with N = 4:
!!
!!  Surfaces ->      S = 1                 S = 2                 S = 3                 S = 4  exit
!!  Lines                                                                                      |
!!  L = 4       13--14--15--16        29--30--31--32        45--46--47--48        61--62--63--64
!!              /   /   /   /         /   /   /   /         /   /   /   /         /   /   /   /
!!  L = 3      9--10--11--12        25--26--27--28        41--42--43--44        57--58--59--60
!!            /   /   /   /         /   /   /   /         /   /   /   /         /   /   /   /
!!  L = 2    5---6---7---8        21--22--23--24        37--38--39--40        53--54--55--56
!!          /   /   /   /         /   /   /   /         /   /   /   /         /   /   /   /
!!  L = 1  1---2---3---4        17--18--19--20        33--34--35--36        49--50--51--52
!!         |
!!       entry
!!
!!  Each node with index p (and thus of rank p - 1) has potentially 6 channels:
!!
!!                                       p+N^2   p+N
!!                                          |   /
!!                                          |  /
!!                                          | /
!!                                          |/
!!                              p-1 --------p-------- p+1
!!                                         /|
!!                                        / |
!!                                       /  |
!!                                      /   |
!!                                   p-N   p-N^2
!!
!!  of which those with index > p are sending channles (p+1, p+N, p+N^2). Some nodes however
!!  do have less than 6 channels like for example node 1, which has only 3 sending channels
!!  2, 5, 17 (1+1, 1+4, 1+16) but no receiving channels. For each given p we need to determine
!!  the appropriate set of channels.
!!
!!  Given N, the node index p belongs to the following surface S and line L number:
!!
!!                                                       Example: N = 4 and p = 43
!!              S = int [(p-1)/N^2] + 1                     S = int[42/16]+1 = 3
!!              L = int [(p-1-N^2{S-1})/N] + 1              L = int[(42-16*2)/4]+1 = 3
!!
!!  The minimum and maximum node indices on each surface and line are given by the following
!!  expressions (included are those for the entire cube):
!!
!!      pmin       = 1                              minimum of entire cube
!!      pmax       = N^3                            maximum of entire cube
!!      pmin (S)   = 1 + N^2 * (S - 1)              minimum on entire surface S
!!      pmax (S)   = N^2 * S                        maximum on entire surface S
!!      pmin (S,L) = pmin (S) + N * (L - 1)         minimum on surface S and line L
!!      pmax (S,L) = pmin (S) + N * L - 1           maximum on surface S and line L
!!
!!  The existence of all 6 channels on a particular node with index p is governed by
!!  the following 6 inequalities:
!!
!!             Channel   p + N^2   exists, if it is   <= pmax
!!                       p - N^2                      >= pmin
!!                       p + N                        <= pmax (S)
!!                       p - N                        >= pmin (S)
!!                       p + 1                        <= pmax (S,L)
!!                       p - 1                        >= pmin (S,L)
!!
!!  Of course, no channels exist if p > N^3. As an example we take node index p = 32 for N = 4.
!!  We have S = 2 and L = 4 and therefore:
!!
!!                  pmin, pmax, pmin (S), pmax (S), pmin (S,L), pmax (S,L) = 
!!                     1,   64,       17,       32,         29,         32
!!
!!  The possible channels to be compared are
!!
!!                 p-N^2,p+N^2,      p-N,      p+N,        p-1,        p+1
!!                    16,   48,       28,       36,         31,         33
!!
!!  From the 6 inequalities, we see that node index p = 32 has only 4 channels: 16,28,31,48 of
!!  which only channel 48 is a sending channel.
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

  use Simulation_data,   ONLY : sim_cubeSide,            &
                                sim_globalMe,            &
                                sim_globalNumProcs,      &
                                sim_numSendingChannels,  &
                                sim_sendingChannels

  implicit none

  integer :: L,n,n2,n3,p,S
  integer :: pmin, pmax, pminS, pmaxS, pminSL, pmaxSL
  integer :: proc

  integer :: channel (1:6)
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

  n  = sim_cubeSide
  p  = pl_rank + 1
  n2 = n * n
  n3 = n * n2
!
!
!     ...Consider only those nodes which make up the cube. Excessive number of nodes
!        results in nodes having no channels.
!
!
  if (p <= n3) then

      S = (p - 1) / n2 + 1
      L = (p - 1 - n2 * (S - 1)) / n + 1

      pmin   = 1
      pmax   = n3
      pminS  = 1 + n2 * (S - 1)
      pmaxS  = n2 * S
      pminSL = pminS + n * (L - 1)
      pmaxSL = pminS + n * L - 1

      if (p >= pminSL + 1) then
          pl_numChannels = pl_numChannels + 1
          channel (pl_numChannels) = p - 2                             ! convert to rank
      end if

      if (p <= pmaxSL - 1) then
          pl_numChannels = pl_numChannels + 1
          sim_numSendingChannels = sim_numSendingChannels + 1
          channel (pl_numChannels) = p                                 ! convert to rank
          sim_sendingChannels (sim_numSendingChannels) = p             ! convert to rank
      end if

      if (p >= pminS + n) then
          pl_numChannels = pl_numChannels + 1
          channel (pl_numChannels) = p - n - 1                         ! convert to rank
      end if

      if (p <= pmaxS - n) then
          pl_numChannels = pl_numChannels + 1
          sim_numSendingChannels = sim_numSendingChannels + 1
          channel (pl_numChannels) = p + n - 1                         ! convert to rank
          sim_sendingChannels (sim_numSendingChannels) = p + n - 1     ! convert to rank
      end if

      if (p >= pmin + n2) then
          pl_numChannels = pl_numChannels + 1
          channel (pl_numChannels) = p - n2 - 1                        ! convert to rank
      end if

      if (p <= pmax - n2) then
          pl_numChannels = pl_numChannels + 1
          sim_numSendingChannels = sim_numSendingChannels + 1
          channel (pl_numChannels) = p + n2 - 1                        ! convert to rank
          sim_sendingChannels (sim_numSendingChannels) = p + n2 - 1    ! convert to rank
      end if

  end if

  if (pl_numChannels > 0) then
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
