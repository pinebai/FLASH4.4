!!****f* source/physics/IncompNS/IncompNS_finalize
!!
!! NAME
!!
!!  IncompNS_finalize
!!
!!
!! SYNOPSIS
!!
!!  IncompNS_finalize()
!!  
!!
!! DESCRIPTION
!! 
!!  Finalize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_finalizeFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine IncompNS_finalize()

  implicit none

end subroutine IncompNS_finalize

