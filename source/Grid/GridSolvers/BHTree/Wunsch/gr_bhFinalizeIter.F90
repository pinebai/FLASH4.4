!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhFinalizeIter
!!
!! NAME
!!
!!  gr_bhFinalizeIter
!!
!! 
!! SYNOPSIS
!!
!!  call gr_bhFinalizeIter()
!!
!!
!! DESCRIPTION
!!
!!  Called by Grid_solvePoisson after all iterations of the tree solver
!!  have finished. Calls corresponding subroutines (*_bhFinalizeIter) of
!!  physical modules (Gravity and TreeRay). 
!!
!!***

subroutine gr_bhFinalizeIter()
  use Gravity_interface, ONLY : Gravity_bhFinalizeIter
  use TreeRay_interface, ONLY : TreeRay_bhFinalizeIter

  implicit none

  call Gravity_bhFinalizeIter()
  call TreeRay_bhFinalizeIter()

  return
end subroutine gr_bhFinalizeIter
