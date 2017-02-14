!!****f* source/physics/sourceTerms/Deleptonize/Deleptonize_init
!!
!! NAME
!!  
!!  Deleptonize_init
!!
!!
!! SYNOPSIS
!! 
!!  call Deleptonize_init()
!!
!!  
!! DESCRIPTION
!!
!!  Perform various initializations (apart from the problem-dependent ones)
!!  for the heat module.
!!
!!
!! ARGUMENTS
!!
!!  
!!
!!***
subroutine Deleptonize_init()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  

  logical :: testUseDeleptonize

!==============================================================================

  !! It is a failure to invoke the stub when useDeleptonize is set TRUE.

  call RuntimeParameters_get ("useDeleptonize", testUseDeleptonize)
  if (testUseDeleptonize) then
     call Driver_abortFlash("Deleptonize unit seems not to be compiled in, and the Deleptonize_init stub does not &
          &allow the value of useDeleptonize to be TRUE.")
  end if

end subroutine Deleptonize_init
