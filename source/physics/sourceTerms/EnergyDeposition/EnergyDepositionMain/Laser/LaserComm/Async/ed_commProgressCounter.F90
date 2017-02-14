!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Async/ed_commProgressCounter
!!
!!  NAME     
!!   ed_commProgressCounter
!!
!!  SYNOPSIS
!!   ed_commProgressCounter()
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

subroutine ed_commProgressCounter()
  use ed_commInterface, ONLY : ed_commCheckTermination
  use ed_commData, ONLY : ed_commCounterState, COUNTER_ACTIVE, COUNTER_INACTIVE
  use UTCounter_sharedCounter, ONLY : UTCounter_progressCounter
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer :: finalCount
  logical :: isCounterTargetMet  

  select case (ed_commCounterState)
  case (COUNTER_ACTIVE)
     call UTCounter_progressCounter(isCounterTargetMet, finalCount)
     if (isCounterTargetMet) then
        ed_commCounterState = COUNTER_INACTIVE
        call ed_commCheckTermination(finalCount)
     end if

  case (COUNTER_INACTIVE)
     continue !no-op

  case default
     call Driver_abortFlash("Unknown counter state")
  end select

end subroutine ed_commProgressCounter
