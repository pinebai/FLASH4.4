!!****if* source/numericalTools/Roots/RootsMain/Roots_init
!!
!! NAME
!!  
!!  Roots_init
!!
!! SYNOPSIS
!! 
!!  call Roots_init ()
!!
!! DESCRIPTION
!!
!!  Initializes the Roots unit.
!!
!! ARGUMENTS
!!
!!***

subroutine Roots_init ()

  use Roots_data

  implicit none
!
!
!     ...Set all machine/precision dependent parameters that control
!        quality/accuracy of the roots.
!
!
  rt_macheps = epsilon (1.0)
  rt_LPN     = huge    (1.0)
  rt_sqrtLPN = sqrt (rt_LPN)
!
!
!    ...Ready!
!
!
  return
end subroutine Roots_init
