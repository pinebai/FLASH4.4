!!****f* source/physics/RadTrans/RadTrans_init
!!
!!  NAME 
!!
!!  RadTrans_init
!!
!!  SYNOPSIS
!!
!!  call RadTrans_init()
!!
!!  DESCRIPTION 
!!    Initialize radiative transfer unit
!!
!! ARGUMENTS
!!
!!
!!***
subroutine RadTrans_init()

  ! Stub implementation

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  logical, save :: testUseRadTrans

  !! It is a failure to invoke the stub when useRadTrans is set TRUE.

  call RuntimeParameters_get ("useRadTrans", testUseRadTrans)
  if (testUseRadTrans) then
     call Driver_abortFlash("RadTrans unit seems not to be compiled in, and the RadTrans_init stub does not &
          &allow the value of useRadTrans to be TRUE.")
  end if

  return
end subroutine RadTrans_init
