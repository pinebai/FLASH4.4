!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhBTInit
!!
!! NAME
!!
!!  gr_bhBTInit
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhBTInit()
!!
!! DESCRIPTION
!!
!!  This subroutine is called at the beginning of the tree build
!!  and it calls subroutine of tree physical modules. They are 
!!  typically used to prepare field arrays for the tree build.
!!
!! ARGUMENTS
!!
!!
!!***


subroutine gr_bhBTInit()

  use Gravity_interface, ONLY : Gravity_bhBTInit
  use TreeRay_interface, ONLY : TreeRay_bhBTInit
  implicit none

  ! call _bhBTInit subroutines of physical modules
  call Gravity_bhBTInit()
  call TreeRay_bhBTInit()

  return
end subroutine gr_bhBTInit


