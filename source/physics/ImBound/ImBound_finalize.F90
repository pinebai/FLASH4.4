!!****f* source/physics/ImBound/ImBound_finalize
!!
!! NAME
!!
!!  Imbound_finalize
!!
!!
!! SYNOPSIS
!!
!!  ImBound_finalize()
!!  
!!
!! DESCRIPTION
!! 
!!  Finalize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_finalizeFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine ImBound_finalize()

  implicit none

end subroutine ImBound_finalize

