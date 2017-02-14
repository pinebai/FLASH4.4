!!****if* source/physics/materialProperties/Opacity/OpacityMain/OPAL/op_finalizeTabulated
!!
!! NAME
!!
!!  op_finalizeTabulated
!!
!! SYNOPSIS
!!
!!  call op_finalizeTabulated ()
!!
!! DESCRIPTION
!!
!!  Finalizes the tabulated section of the opacity unit.
!!
!! ARGUMENTS
!!
!!***
subroutine op_finalizeTabulated ()

  use Opacity_data,                ONLY : op_opalNumHydrogenAbundances

  use op_opalData, ONLY : OP_OPAL_LOWT,OP_OPAL_HIGHT,       &
                          op_tableKind,                 &
                               op_tableNameLowT,                 &
                               op_tableNameHighT,                 &
                               op_opalAllTab

  implicit none

  integer :: m, tempRange
!
!
!   ...Deallocate the arrays.
!
!
  if (allocated (op_tableKind)                ) deallocate (op_tableKind)
  if (allocated (op_tableNameLowT)                ) deallocate (op_tableNameLowT)
  if (allocated (op_tableNameHighT)                ) deallocate (op_tableNameHighT)

  do tempRange = OP_OPAL_LOWT,OP_OPAL_HIGHT
     do m = 1,op_opalNumHydrogenAbundances
        ! Currently, the only way these pointers get different from null is by being allocated. - KW
        if (associated(op_opalAllTab(m,tempRange)%tg%table)) then
              if (associated(op_opalAllTab(m,tempRange)%tg%table%table)) &
                   deallocate(op_opalAllTab(m,tempRange)%tg%table%table)
           deallocate(op_opalAllTab(m,tempRange)%tg%table)
        end if
        if (associated(op_opalAllTab(m,tempRange)%tg%mgTable)) deallocate(op_opalAllTab(m,tempRange)%tg%mgTable)
        if (associated(op_opalAllTab(m,tempRange)%tg%td%Temperatures)) deallocate(op_opalAllTab(m,tempRange)%tg%td%Temperatures)
        if (associated(op_opalAllTab(m,tempRange)%tg%td%Densities)) deallocate(op_opalAllTab(m,tempRange)%tg%td%Densities)
     end do
  end do
  deallocate(op_opalAllTab)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_finalizeTabulated
