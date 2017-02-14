!!****if* source/physics/materialProperties/Opacity/localAPI/op_opalInterface
!!
!! NAME
!!
!!  op_opalInterface
!!
!! SYNOPSIS
!!
!!  use op_opalInterface
!!
!! DESCRIPTION
!!
!!  Defines interfaces for some local routines for OPAL tabulated opacities.
!!  
!!***

module op_opalInterface
  
  interface
     subroutine op_readOpalTables (tableNameLowT,   &
          tableNameHighT,   &
          needLowTable, &
          needHighTable, &
          needROTable,      &
          td, &
          tbLowT,tbHighT  )
       use op_opalData, ONLY : opT_tableGroupDescT,  &
                               opT_oneVarTablePT,    &
                               OP_OPAL_LOWT
       
       implicit none
       character (len=80), intent (in) :: tableNameLowT
       character (len=80), intent (in) :: tableNameHighT
       logical,            intent (in) :: needLowTable
       logical,            intent (in) :: needHighTable
       logical,            intent (in) :: needROTable
       type(opT_tableGroupDescT),intent(inout) :: td(OP_OPAL_LOWT:)
       type(opT_oneVarTablePT),pointer :: tbLowT,tbHighT
       end
  end interface

  interface
!!$     subroutine op_readOpalTable (tableName,   &
!!$                                      needPATable, &
!!$                                      needPETable, &
!!$                                      needROTable, &
!!$                                      indexPA,     &
!!$                                      indexPE,     &
!!$                                      indexRO      )
!!$
!!$       character (len=80), intent (in) :: tableName
!!$       logical,            intent (in) :: needPATable
!!$       logical,            intent (in) :: needPETable
!!$       logical,            intent (in) :: needROTable
!!$       integer,            intent (in) :: indexPA
!!$       integer,            intent (in) :: indexPE
!!$       integer,            intent (in) :: indexRO
!!$     end subroutine op_readOpalTable
     subroutine op_readOpalTable (tableName,   &
                                      td, &
                                      tb      )
       use op_opalData, ONLY : opT_tableGroupDescT,  &
                               opT_oneVarTablePT
       implicit none
       character (len=80), intent (in) :: tableName
       type(opT_tableGroupDescT),intent(inout) :: td
       type(opT_oneVarTablePT),pointer :: tb
     end subroutine op_readOpalTable
  end interface

end module op_opalInterface

