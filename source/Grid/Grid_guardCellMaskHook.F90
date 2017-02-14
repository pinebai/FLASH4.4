!!****f* source/Grid/Grid_guardCellMaskHook
!!
!! NAME
!!
!!  Grid_guardCellMaskHook
!!
!! SYNOPSIS
!!
!!  call Grid_guardCellMaskHook(logical,intent(INOUT)  :: ccmask,
!!                              logical,intent(IN)  :: needeos)
!!
!! DESCRIPTION
!!  
!! Stub
!!
!! ARGUMENTS
!!
!!   ccmask : the mask
!!
!!   needeos : switch for the need of Eos
!!
!!
!!
!!***

subroutine Grid_guardCellMaskHook(ccMask, needEos)
  implicit none
  logical,intent(INOUT) :: ccMask(*)
  logical,intent(IN)    :: needEos

end subroutine Grid_guardCellMaskHook
