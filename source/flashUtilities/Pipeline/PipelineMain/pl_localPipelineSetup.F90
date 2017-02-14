!!****if* source/flashUtilities/Pipeline/PipelineMain/pl_localPipelineSetup
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

  use Pipeline_data,     ONLY : pl_doLog,       &
                                pl_logUnit,     &
                                pl_numChannels, &
                                pl_procList,    &
                                pl_rank

  use Grid_interface,    ONLY : Grid_getNeighProcList

  implicit none

  logical, parameter :: includeMyProc = .false.
  integer, pointer   :: tmpNeighborsProcList (:)

  integer :: channel
  integer :: numNeighbors
  integer :: proc
!
!
!     ...The default is to set up the pipeline structure from the current structure
!        of the grid. Obtain the proc IDs of all surrounding neighbor blocks in
!        order to initialize the communication pipeline.  Note that the Grid subroutine
!        will return a null pointer if all neighboring blocks exist on my MPI rank and
!        includeMyProc is set to false. If all neighboring blocks exist on my MPI rank
!        then we set a placeholder 1-element array.
!
!
  call Grid_getNeighProcList (includeMyProc, tmpNeighborsProcList, numNeighbors)

  if (numNeighbors > 0) then
      if (.not. associated (tmpNeighborsProcList)) then
           call Driver_abortFlash ("[pl_localPipelineSetup] Neighbor processor list not associated.?")
      end if

      pl_numChannels = numNeighbors
      allocate (pl_procList (1:pl_numChannels))
      pl_procList (:) = tmpNeighborsProcList (:)
  else
      pl_numChannels = 0
      allocate (pl_procList (1:1))
      pl_procList (1) = -1
  end if

  if (associated (tmpNeighborsProcList)) deallocate (tmpNeighborsProcList)
  nullify (tmpNeighborsProcList)
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
          do channel = 1,pl_numChannels
             proc = pl_procList (channel)
             write (pl_logUnit,'(a,i4,a,i6)'   ) '        Channel ',channel,' -> Processor ',proc 
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
