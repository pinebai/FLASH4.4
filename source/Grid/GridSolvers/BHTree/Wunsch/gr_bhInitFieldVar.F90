!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhInitFieldVar
!!
!! NAME
!!
!!  gr_bhInitFieldVar
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhInitFieldVar(
!!                 integer(in) :: gpotVar
!!  )
!!
!! DESCRIPTION
!!
!!  Called before the tree walk. Calls corresponding subroutines of physical
!!  modules to initialize field variables needed for calculation of the tree
!!  code. Passes variable gpotVar (where the potential is stored obtained 
!!  from Grid_solvePoisson - ipotvar.
!!
!! ARGUMENTS
!!
!!  gpotVar : index of the grid variable where the gravitational potential is
!!            stored. Passed from Grid_solvePoisson.
!!
!!
!!***

subroutine gr_bhInitFieldVar(gpotVar)
  use Gravity_interface, ONLY : Gravity_bhInitFieldVar
  use TreeRay_interface, ONLY : TreeRay_bhInitFieldVar
  implicit none
  integer, intent(IN) :: gpotVar

  ! call _bhInitFieldVar subroutines of physical modules
  call Gravity_bhInitFieldVar(gpotVar)
  call TreeRay_bhInitFieldVar(gpotVar)

  return
end subroutine gr_bhInitFieldVar
