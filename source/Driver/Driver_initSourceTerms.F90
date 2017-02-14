!!****f* source/Driver/Driver_initSourceTerms
!!
!! NAME
!!   
!!  Driver_initSourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_initSourceTerms(logical(in) :: restart)
!!
!! DESCRIPTION
!!
!!   Initializes all source terms Units by 
!!   calling their respective initialization routines,
!!   viz. Stir_init, Burn_init, Heat_init, Cool_init, etc.
!!  
!!
!! ARGUMENTS
!!   restart - indicates if run is starting from scratch (.false.)
!!             or restarting from checkpoint (.true.)
!!
!!***

subroutine Driver_initSourceTerms( restart)

  implicit none

  logical, intent(in) :: restart

end subroutine Driver_initSourceTerms
