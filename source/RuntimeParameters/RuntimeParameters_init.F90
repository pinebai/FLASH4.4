!!****f* source/RuntimeParameters/RuntimeParameters_init
!!
!! NAME
!!  RuntimeParameters_init
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_init(logical(out) :: restart)
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the data need in the Runtime Parameters
!!  Unit.  Reads the parameter file, usually flash.par,
!!  and broadcasts parameters to the other processors.
!!    
!!
!! ARGUMENTS
!!
!! restart - true if run is restarted from checkpoint, false if starting
!!           from scratch
!!           
!!
!!
!!***



subroutine RuntimeParameters_init( restart)

implicit none

  logical, intent(out) :: restart

  restart = .false.

end subroutine RuntimeParameters_init
