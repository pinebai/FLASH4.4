!!****if* source/physics/materialProperties/Opacity/OpacityMain/OPAL/op_readOpalTables
!!
!! NAME
!!
!!  op_readOpalTables
!!
!! SYNOPSIS
!!
!!  call op_readOpalTables (character (in) :: tableNameLowT (len=80),
!!                      character (in) :: tableNameHighT (len=80),
!!                      logical   (in) :: needLowTable,
!!                      logical   (in) :: needHighTable,
!!                      logical   (in) :: needROTable)
!!
!! DESCRIPTION
!!
!!  Reads in the necessary data for processing and interpolating tabulated
!!  opacities from specific tables. The routine calls appropriate subroutines
!!  according to the kind of opacity tables specified. Currently only the
!!  IONMIX opacities can be read in.
!!
!! ARGUMENTS
!!
!!  tableNameLowT  : the name of tabulated Opacity file where low temp data is going to be read
!!  tableNameHighT : the name of tabulated Opacity file where high temp data is going to be read
!!  needLowTable : if yes, Planck Absorption Opacities are needed from the Opacity table
!!  needHighTable : if yes, Planck   Emission Opacities are needed from the Opacity table
!!  needROTable : if yes,         Rosseland Opacities are needed from the Opacity table
!!
!!***
subroutine op_readOpalTables (tableNameLowT,   &
                          tableNameHighT,   &
                          needLowTable, &
                          needHighTable, &
                          needROTable,      &
                              td, &
                              tbLowT,tbHighT  )

  use Driver_interface,   ONLY : Driver_abortFlash
  use op_opalInterface,   ONLY : op_readOpalTable
  use op_opalData,        ONLY : OP_OPAL_LOWT,         &
                                 OP_OPAL_HIGHT,        &
                                 opT_tableGroupDescT,  &
                                 opT_oneVarTablePT

  implicit none

  character (len=80), intent (in) :: tableNameLowT
  character (len=80), intent (in) :: tableNameHighT
  logical,            intent (in) :: needLowTable
  logical,            intent (in) :: needHighTable
  logical,            intent (in) :: needROTable
  type(opT_tableGroupDescT),intent(inout) :: td(OP_OPAL_LOWT:)
  type(opT_oneVarTablePT),pointer :: tbLowT,tbHighT
!
!
!   ...Call the appropriate routine.
!
!
#ifdef DEBUG_OPACITY
  print*,'op_readOpalTables NEEDS:',needROTable, needLowTable, needHighTable
#endif

  if (needROTable) then
    if (needLowTable ) then
      call op_readOpalTable (tableNameLowT,  td(OP_OPAL_LOWT) , tbLowT)
    end if

    if (needHighTable ) then
      call op_readOpalTable (tableNameHighT, td(OP_OPAL_HIGHT), tbHighT)
    end if

  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_readOpalTables
