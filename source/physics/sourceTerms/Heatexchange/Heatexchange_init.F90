!!****f* source/physics/sourceTerms/Heatexchange/Heatexchange_init
!!
!! NAME
!!
!!  Heatexchange_init
!!
!! SYNOPSIS
!!
!!  call Heatexchange_init(logical(IN)  :: restart)
!!
!! DESCRIPTION
!!
!!   Initializes the Heatexchange unit.
!!
!! ARGUMENTS
!!
!!   restart : indicates whether the run is restarting from a checkpoint
!!
!! PARAMETERS
!!
!!  useHeatexchange -- Boolean, True.  Turns on Heatexchange unit
!!
!!  Other runtime parameters are specific to certain implementations.
!!
!!***

subroutine Heatexchange_init(restart)
  implicit none
  
  logical, intent(IN) :: restart


end subroutine Heatexchange_init
