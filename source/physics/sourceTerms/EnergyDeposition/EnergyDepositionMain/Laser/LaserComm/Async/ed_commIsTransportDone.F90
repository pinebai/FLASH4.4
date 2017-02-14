!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Async/ed_commIsTransportDone
!!
!!  NAME     
!!   ed_commIsTransportDone
!!
!!  SYNOPSIS
!!   ed_commIsTransportDone(logical, intent(OUT) :: isTransportDone)
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

subroutine ed_commIsTransportDone(isTransportDone)
  use ed_commData, ONLY : ed_commPipelineState, ed_commCounterState, &
       PIPELINE_INACTIVE, COUNTER_INACTIVE
  implicit none
  logical, intent(OUT) :: isTransportDone

  isTransportDone = &
       ( COUNTER_INACTIVE == ed_commCounterState .and. &
       PIPELINE_INACTIVE == ed_commPipelineState )
end subroutine ed_commIsTransportDone
