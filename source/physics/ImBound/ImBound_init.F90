!!****f* source/physics/ImBound/ImBound_init
!!
!! NAME
!!
!!  ImBound_init
!!
!!
!! SYNOPSIS
!!
!!  call ImBound_init(LOGICAL(IN) :: restart)
!!  
!! ARGUMENTS
!!
!!  restart - restart flag.
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which typically take values from runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine ImBound_init(restart)

  implicit none

  logical, INTENT(IN) :: restart

  return

end subroutine ImBound_init
