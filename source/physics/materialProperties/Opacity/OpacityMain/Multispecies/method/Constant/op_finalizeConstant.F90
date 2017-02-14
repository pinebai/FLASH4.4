!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Constant/op_finalizeConstant
!!
!! NAME
!!
!!  op_finalizeConstant
!!
!! SYNOPSIS
!!
!!  call op_finalizeConstant ()
!!
!! DESCRIPTION
!!
!!  Finalizes the constant section of the opacity unit.
!!
!! ARGUMENTS
!!
!!***
subroutine op_finalizeConstant ()

  use op_constantData, ONLY : op_absorptionConstant,  &
                              op_emissionConstant,    &
                              op_transportConstant

  implicit none
!
!
!   ...Deallocate the arrays.
!
!
  if (allocated (op_absorptionConstant)) deallocate (op_absorptionConstant)
  if (allocated (op_emissionConstant)  ) deallocate (op_emissionConstant)
  if (allocated (op_transportConstant) ) deallocate (op_transportConstant)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_finalizeConstant
