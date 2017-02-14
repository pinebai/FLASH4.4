!!****if* source/Simulation/SimulationMain/RTFlame/Simulation_guardCellMaskHook
!!
!! NAME
!!
!!  Simulation_guardCellMaskHook
!!
!! SYNOPSIS
!!
!!  call Simulation_guardCellMaskHook(logical(INOUT)  :: ccmask,
!!                                    logical(IN)  :: needeos)
!!
!! DESCRIPTION
!!
!!  A hook that lets a simulation modify the mask to use for guard cell filling.
!!
!!  Indirectly called from gr_makeMaskConsistent, which may get called from
!!  Grid_fillGuardCells (depending on the arguments with which Grid_fillGuardCells
!!  is called).
!!
!! ARGUMENTS
!!
!!   ccmask : the mask
!!
!!   needeos : switch for the need of Eos
!!
!!
!!***

#include "Flash.h"

subroutine Simulation_guardCellMaskHook(ccMask, needEos)
  implicit none
  logical,intent(INOUT) :: ccMask(*)
  logical,intent(IN)    :: needEos

  !!  Additional logic necessary due to the hydrostatic boundary conditions.
  ! For constant isothermal BC we need temperature and material info (ye, sumy) in interior.
  ! This temperature setting does not trigger the "needEos" flag.
  if ( ccMask(DENS_VAR) .or. ccMask(EINT_VAR) .or. ccMask(TEMP_VAR) &
       .or. ccMask(PRES_VAR) .or. ccMask(ENER_VAR) &
       .or. ccMask(GAMC_VAR) .or. ccMask(GAME_VAR) ) then
     ccMask(DENS_VAR) = .true.
     ccMask(TEMP_VAR) = .true.
     ! material information (ye, yi) stored in flame variable
     ccMask(FLAM_MSCALAR) = .true.
  endif

end subroutine Simulation_guardCellMaskHook

