!!****f* source/physics/TreeRay/TreeRay_init
!!
!! NAME
!!
!!  TreeRay_init
!!  
!! SYNOPSIS
!!
!!  TreeRay_init()
!!
!! DESCRIPTION
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!! 
!!  useTreeRay  BOOLEAN true  Controls turning on/off the compiled gravity unit
!!
!! NOTES
!!   
!!  Each implementation of TreeRay has its own runtime parameters.  Be sure to check
!!  the documentation or Config files to see them.
!!
!!***

subroutine TreeRay_init ()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
   

  logical, save :: testUseTreeRay

!==============================================================================

  !! It is a failure to invoke the stub when useTreeRay is set TRUE.

  call RuntimeParameters_get ("useTreeRay", testUseTreeRay)
  if (testUseTreeRay) then
     call Driver_abortFlash("TreeRay unit seems not to be compiled in, and the TreeRay_init stub does not &
          &allow the value of useTreeRay to be TRUE.")
  end if

end subroutine TreeRay_init
