!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Integrate/op_finalizeIntegrate
!!
!! NAME
!!
!!  op_finalizeIntegrate
!!
!! SYNOPSIS
!!
!!  call op_finalizeIntegrate ()
!!
!! DESCRIPTION
!!
!!  Finalizes the integration section of the opacity unit.
!!
!! ARGUMENTS
!!
!!***
subroutine op_finalizeIntegrate ()

  use op_integrateData, ONLY : op_AuxPolynomialA,      &
                               op_AuxPolynomialB,      &
                               op_Moments,             &
                               op_JmatrixDiagonals,    &
                               op_JmatrixOffdiagonals, &
                               op_OrthoPolynomialA,    &
                               op_OrthoPolynomialB,    &
                               op_work1,               &
                               op_work2,               &
                               op_work3,               &
                               op_RombergIntegral,     &
                               op_RombergRow

  implicit none
!
!
!   ...Deallocate the quadradture arrays (if any).
!
!
  if (allocated (op_AuxPolynomialA)     ) deallocate (op_AuxPolynomialA)
  if (allocated (op_AuxPolynomialB)     ) deallocate (op_AuxPolynomialB)
  if (allocated (op_Moments)            ) deallocate (op_Moments)
  if (allocated (op_JmatrixDiagonals)   ) deallocate (op_JmatrixDiagonals)
  if (allocated (op_JmatrixOffdiagonals)) deallocate (op_JmatrixOffdiagonals)
  if (allocated (op_OrthoPolynomialA)   ) deallocate (op_OrthoPolynomialA)
  if (allocated (op_OrthoPolynomialB)   ) deallocate (op_OrthoPolynomialB)
  if (allocated (op_work1)              ) deallocate (op_work1)
  if (allocated (op_work2)              ) deallocate (op_work2)
  if (allocated (op_work3)              ) deallocate (op_work3)
!
!
!   ...Deallocate the romberg arrays (if any).
!
!
  if (allocated (op_RombergIntegral)    ) deallocate (op_RombergIntegral)
  if (allocated (op_RombergRow)         ) deallocate (op_RombergRow)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_finalizeIntegrate
