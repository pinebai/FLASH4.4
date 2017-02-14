!!****if* source/physics/sourceTerms/Flame/FlameEffects/EIP/Flame_computeAbarZbar
!!
!! NAME
!!
!!  Flame_computeAbarZbar
!!
!! SYNOPSIS
!!
!!  call Flame_computeAbarZbar(real, dimension(:,:)(in) :: solnscalars,
!!                             real, dimension(:)(inout) :: abardata,
!!                             real, dimension(:)(inout) :: zbardata)
!!
!! DESCRIPTION
!! 
!! Return abar and zbar
!!
!! ARGUMENTS
!!
!!   solnscalars : solution scalars
!!
!!   abardata : abar 
!!
!!   zbardata : zbar 
!!
!!
!!
!!***

#include "Flash.h"
subroutine Flame_computeAbarZbar(solnScalars, abarData, zbarData)

  use fl_effData, only : fl_eff_ye_u, fl_eff_ye_b, fl_eff_sumy_u, fl_eff_sumy_b

  implicit none

  real, dimension(:,:), intent(in)  :: solnScalars
  real, dimension(:), intent(inout) :: abarData, zbarData

  ! index of flam variable in scalars array
  integer :: flami = FLAM_MSCALAR-SPECIES_BEGIN+1

  real :: flam, ye, yi
  integer i

  do i = 1, ubound(solnScalars, 2)
     flam = solnScalars(flami,i)
     abarData(i) = 1.0 / ( (1.0-flam)*fl_eff_sumy_u + flam*fl_eff_sumy_b )
     zbarData(i) = abarData(i)* ( (1.0-flam)*fl_eff_ye_u + flam*fl_eff_ye_b )
  enddo

end subroutine Flame_computeAbarZbar
