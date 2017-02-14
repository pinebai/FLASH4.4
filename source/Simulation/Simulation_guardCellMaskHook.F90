!!****f* source/Simulation/Simulation_guardCellMaskHook
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


subroutine Simulation_guardCellMaskHook(ccMask, needEos)
  implicit none
  logical,intent(INOUT) :: ccMask(*)
  logical,intent(IN)    :: needEos

  ! Stub does nothing.
end subroutine Simulation_guardCellMaskHook

