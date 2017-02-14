!!****if* source/physics/sourceTerms/Turb/TurbMain/Turb_getFilterScale
!!
!! NAME
!!
!!  Turb_getFilterScale
!!
!! SYNOPSIS
!!
!!  call Turb_getFilterScale(real(in) :: dx,
!!                           real(out) :: de)
!!
!! DESCRIPTION
!!
!! Aaron Jackson 2010
!! returns filter scale for turbulence operator.
!! The current implementation does not use dx, but rather the
!! minimum dx in the simulation. This assumes that the highest
!! refinement is used to resolve the flame.
!! The filter scale associated with the turbulence operator is
!! returned in de.
!!
!! ARGUMENTS
!!
!!   dx : minimum dx of the simulation
!!
!!   de : filter scale 
!!
!!
!!
!!***


#include "constants.h"
#include "Flash.h"
subroutine Turb_getFilterScale(dx, de)

  use Turb_data, only : turb_stepSize
  use Grid_interface, only : Grid_getMinCellSize

  implicit none

  real, intent(in) :: dx
  real, intent(out) :: de

  real :: dx_min

  ! For this implementation, we only care about the resolution used
  ! to resolve the flame which is the minimum
  call Grid_getMinCellSize(dx_min)

  de = 4.0e0 * turb_stepSize * dx_min

  return
end subroutine Turb_getFilterScale
