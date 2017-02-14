!!****f* source/physics/Gravity/Gravity_init
!!
!! NAME
!!
!!  Gravity_init
!!  
!! SYNOPSIS
!!
!!  Gravity_init()
!!
!! DESCRIPTION
!!
!!  Initialize unit scope variables in the Gravity unit, which are typically the 
!!  runtime parameters.  This routine must be called once by Driver_initFlash.F90. 
!!  Calling multiple times will not cause any harm but is unnecessary.
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!! 
!!  useGravity  BOOLEAN true  Controls turning on/off the compiled gravity unit
!!
!! NOTES
!!   
!!  Each implementation of Gravity has its own runtime parameters.  Be sure to check
!!  the documentation or Config files to see them.
!!
!!***

subroutine Gravity_init ()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
   

  logical, save :: testUseGravity

!==============================================================================

  !! It is a failure to invoke the stub when useGravity is set TRUE.

  call RuntimeParameters_get ("useGravity", testUseGravity)
  if (testUseGravity) then
     call Driver_abortFlash("Gravity unit seems not to be compiled in, and the Gravity_init stub does not &
          &allow the value of useGravity to be TRUE.")
  end if

end subroutine Gravity_init
