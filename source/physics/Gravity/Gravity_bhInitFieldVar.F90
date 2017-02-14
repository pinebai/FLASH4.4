!!****f* source/physics/Gravity/Gravity_bhInitFieldVar
!!
!! NAME
!!
!!  Gravity_bhInitFieldVar
!!
!!
!! SYNOPSIS
!!
!!  call Gravity_bhInitFieldVar(
!!                 integer(in) :: gpotVar
!!  )
!!
!! DESCRIPTION
!!
!!  Called before the tree walk. Initializes fiels variables befor the potential
!!  is calculated. Sets index of the grid variable where the gravitational 
!!  potential is stored (passed from Grid_solvePoisson - ipotvar. Calculate
!!  gravitational acceleration at each grid cell (inverted maximum allowed error
!!  in acceleration stored in field variable ACEI) - needed for normalization of 
!!  the error.
!!
!!
!! ARGUMENTS
!!
!!  gpotVar : index of the grid variable where the gravitational potential is
!!            stored. Passed from Grid_solvePoisson.
!!
!!
!!***

subroutine Gravity_bhInitFieldVar(gpotVar)
  implicit none
#include "constants.h"
#include "Flash.h"
  integer, intent(IN) :: gpotVar
  
  return
end subroutine Gravity_bhInitFieldVar



