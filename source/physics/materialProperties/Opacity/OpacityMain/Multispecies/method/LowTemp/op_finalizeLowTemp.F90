!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_finalizeLowTemp
!!
!! NAME
!!
!!  op_finalizeLowTemp
!!
!! SYNOPSIS
!!
!!  call op_finalizeLowTemp ()
!!
!! DESCRIPTION
!!
!!  Finalizes the low temperature section of the opacity unit.
!!
!! ARGUMENTS
!!
!!***
subroutine op_finalizeLowTemp ()

  use op_lowTempData,  ONLY : op_A1group,               &
                              op_A2group,               &
                              op_A3group,               &
                              op_A4group,               &
                              op_intLimits,             &
                              op_Aij4,                  &
                              op_Jmax,                  &
                              op_PEenergyRange,         &
                              op_elementAij4,           &
                              op_elementJmax,           &
                              op_elementPEenergyRange,  &
                              op_tableLowTemp,          &
                              op_PlanckLowTempTables,   &
                              op_RosselandLowTempTables
  implicit none
!
!
!   ...Deallocate the arrays.
!
!
  if (allocated (op_A1group)               ) deallocate (op_A1group)
  if (allocated (op_A2group)               ) deallocate (op_A2group)
  if (allocated (op_A3group)               ) deallocate (op_A3group)
  if (allocated (op_A4group)               ) deallocate (op_A4group)
  if (allocated (op_intLimits)             ) deallocate (op_intLimits)

  if (allocated (op_Aij4)                  ) deallocate (op_Aij4)
  if (allocated (op_Jmax)                  ) deallocate (op_Jmax)
  if (allocated (op_PEenergyRange)         ) deallocate (op_PEenergyRange)

  if (allocated (op_elementAij4)           ) deallocate (op_elementAij4)
  if (allocated (op_elementJmax)           ) deallocate (op_elementJmax)
  if (allocated (op_elementPEenergyRange)  ) deallocate (op_elementPEenergyRange)

  if (allocated (op_tableLowTemp)          ) deallocate (op_tableLowTemp)
  if (allocated (op_PlanckLowTempTables)   ) deallocate (op_PlanckLowTempTables)
  if (allocated (op_RosselandLowTempTables)) deallocate (op_RosselandLowTempTables)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_finalizeLowTemp
