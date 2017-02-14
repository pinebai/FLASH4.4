!!****f* source/Driver/Driver_finalizeSourceTerms
!!
!! NAME
!!   
!!  Driver_finalizeSourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_finalizeSourceTerms(logical(in) :: restart)
!!
!! DESCRIPTION
!!
!!   Finalizes all source terms Units by 
!!   calling their respective termination  routines,
!!   viz. Stir_finalize, Burn_finalize, Heat_finalize, Cool_finalize, etc.
!!  
!!
!! ARGUMENTS
!!   restart - indicates if run is starting from scratch (.false.)
!!             or restarting from checkpoint (.true.)
!!
!!***

subroutine Driver_finalizeSourceTerms( restart)

  implicit none

  logical, intent(in) :: restart

end subroutine Driver_finalizeSourceTerms
