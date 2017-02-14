!!****f* source/physics/Cosmology/Cosmology_init
!!
!! NAME
!!
!!  Cosmology_init
!!
!!
!! SYNOPSIS
!!
!!  Cosmology_init(logical(IN) :: restart)
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!! ARGUMENTS
!!
!!  
!!  restart -- note if we are restarting. 
!!
!! NOTE
!!
!!  From a "cold start" this unit has to determine the cosmological scaling 
!!  and will solve the Friedmann equations to do so.  On a restart, the
!!  correct values will be read from the checkpoint file
!!
!!***

subroutine Cosmology_init(restart)

  implicit none
  logical, intent(IN) :: restart
  
end subroutine Cosmology_init
