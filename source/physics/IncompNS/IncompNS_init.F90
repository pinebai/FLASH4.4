!!****f* source/physics/IncompNS/IncompNS_init
!!
!! NAME
!!
!!  IncompNS_init
!!
!!
!! SYNOPSIS
!!
!!  call IncompNS_init(restart)
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine IncompNS_init(restart)

  implicit none
  logical, intent(IN) :: restart

end subroutine IncompNS_init

