!!****f* source/physics/SolidMechanics/SolidMechanics
!!
!!
!! NAME
!!
!! 
!!
!!
!! SYNOPSIS
!!
!!  
!!
!!
!! DESCRIPTION
!!
!!
!!
!!***

subroutine SolidMechanics(selector_flag, restart, convflag_all)
  implicit none
#include "SolidMechanics.h"
  integer, intent(in) :: selector_flag
  logical, intent(in),optional :: restart
  integer, optional, intent(inout) :: convflag_all 

  if (present(convflag_all)) convflag_all = SM_CONVERGED

  return
end subroutine SolidMechanics
