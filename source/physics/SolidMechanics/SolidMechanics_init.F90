!!****f* source/physics/SolidMechanics/SolidMechanics_init
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine SolidMechanics_init(restart)

  implicit none

  logical, INTENT(IN) :: restart

  return

end subroutine SolidMechanics_init
