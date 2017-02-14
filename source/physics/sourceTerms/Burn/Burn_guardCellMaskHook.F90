!!****f* source/physics/sourceTerms/Burn/Burn_guardCellMaskHook
!!
!! NAME
!!
!!  Burn_guardCellMaskHook
!!
!! SYNOPSIS
!!
!!  call Burn_guardCellMaskHook(logical,intent(INOUT)  :: ccMask,
!!                              logical,intent(IN)  :: needEos)
!!
!! DESCRIPTION
!! GC mask hook
!!
!! ARGUMENTS
!!
!!   ccmask : cc mask
!!
!!   needeos : check for eos
!!
!!
!!
!!***

subroutine Burn_guardCellMaskHook(ccMask, needEos)
  implicit none
  logical,intent(INOUT) :: ccMask(*)
  logical,intent(IN)    :: needEos

end subroutine Burn_guardCellMaskHook

