!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Async/ed_commCheckTermination
!!
!!  NAME     
!!   ed_commCheckTermination
!!
!!  SYNOPSIS
!!   ed_commCheckTermination()
!!
!!  DESCRIPTION 
!!
!!  ARGUMENTS
!!
!!  NOTES
!!
!!  SIDE EFFECTS
!!
!!***

subroutine ed_commCheckTermination(finalCount)
  use ed_commData, ONLY : ed_commLog, ed_commLogUnit, ed_commGlobalRays
  use Logfile_interface, ONLY : Logfile_open, Logfile_close
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer, intent(IN) :: finalCount
  character(len=*), parameter :: badTerminationMsg = &
       'Initial active ray count /= Final inactive ray count.'
  logical :: isBadTermination

  isBadTermination = (finalCount /= ed_commGlobalRays)
  if (isBadTermination) then
     if (.not.ed_commLog) call Logfile_open(ed_commLogUnit,.true.)
     write(ed_commLogUnit,'(a, 2(/,a,i12))') &
          'Bad inactive ray count.', &
          'Initial active ray count ', ed_commGlobalRays, &
          'Final inactive ray count ', finalCount
     call Logfile_close(.true.)
     call Driver_abortFlash(badTerminationMsg)
  end if
end subroutine ed_commCheckTermination
