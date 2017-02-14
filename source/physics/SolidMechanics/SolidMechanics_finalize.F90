!!****f* source/physics/SolidMechanics/SolidMechanics_finalize
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!  
!!
!! DESCRIPTION
!! 
!!  Finalize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_finalizeFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine SolidMechanics_finalize()

  implicit none

end subroutine SolidMechanics_finalize

