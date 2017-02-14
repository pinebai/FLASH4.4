!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Constcm2g/op_finalizeConstcm2g
!!
!! NAME
!!
!!  op_finalizeConstcm2g
!!
!! SYNOPSIS
!!
!!  call op_finalizeConstcm2g ()
!!
!! DESCRIPTION
!!
!!  Finalizes the constant-cm2g section of the opacity unit.
!!
!! ARGUMENTS
!!
!!***
subroutine op_finalizeConstcm2g ()

  use op_constcm2gData, ONLY : op_absorptionConstcm2g,  &
                               op_emissionConstcm2g,    &
                               op_transportConstcm2g

  implicit none
!
!
!   ...Deallocate the arrays.
!
!
  if (allocated (op_absorptionConstcm2g)) deallocate (op_absorptionConstcm2g)
  if (allocated (op_emissionConstcm2g)  ) deallocate (op_emissionConstcm2g)
  if (allocated (op_transportConstcm2g) ) deallocate (op_transportConstcm2g)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_finalizeConstcm2g
