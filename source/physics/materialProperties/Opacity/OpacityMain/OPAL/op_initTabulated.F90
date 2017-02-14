!!****if* source/physics/materialProperties/Opacity/OpacityMain/OPAL/op_initTabulated
!!
!! NAME
!!
!!  op_initTabulated
!!
!! SYNOPSIS
!!
!!  call op_initTabulated ()
!!
!! DESCRIPTION
!!
!!  Inititalizes the tabulated section of the opacity unit.
!!
!! ARGUMENTS
!!
!!***
subroutine op_initTabulated ()

  use Driver_interface,            ONLY : Driver_abortFlash

  use Opacity_data,                ONLY : op_opalNumHydrogenAbundances,              &
                                          op_nEnergyGroups,             &
                                          op_absorptionKind,            &
                                          op_emissionKind,              &
                                          op_transportKind,             &
                                          op_writeOpacityInfo

  use op_opalData,            ONLY :  OP_OPAL_LOWT,OP_OPAL_HIGHT,       &
                                          op_initializedTabulated,      &
                                          op_useLogTables,              &
                                          op_maxTablesLowT,               &
                                          op_maxTablesHighT,               &
                                          op_maxTablesRO,               &
                                          op_tableKind,                 &
                                          op_tableNameLowT,                 &
                                          op_tableNameHighT,                 &
                                          op_opalTableAbundMax,         &
                                          op_opalAllTab

  use op_interface,                ONLY : op_browseTables,              &
                                          op_writeTables
  use op_opalInterface,            ONLY : op_readOpalTables

  use Simulation_interface,        ONLY : Simulation_mapIntToStr

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt

  implicit none

#include "Flash.h"
#include "Opacity.h"
#include "OpacityOPAL.h"
#include "constants.h"

  real, parameter :: ten   = 10.0

  character (len=20) :: spec_str
  character (len=MAX_STRING_LENGTH) :: str
  character (len=MAX_STRING_LENGTH) :: rtpar

  logical :: needLowTable
  logical :: needHighTable
  logical :: needROTable
  logical :: needTable

  integer :: i
  integer :: m, tempRange

  integer :: absorptionKind
  integer :: emissionKind
  integer :: transportKind
  integer :: nstepsDensityLowT
  integer :: nstepsDensityHighT
  integer :: nstepsDensityRO
  integer :: nstepsTemperatureLowT
  integer :: nstepsTemperatureHighT
  integer :: nstepsTemperatureRO
  integer :: status
!
!
!    ...Safety net. Runtime parameters.
!
!
  if (op_opalNumHydrogenAbundances < 1) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: No species present!')
  end if

  call RuntimeParameters_get ("opacity_useLogTables",   op_useLogTables)

!
!
!    ...Allocate those arrays that depend on the # of species.
!
!
  allocate (op_tableKind (1:op_opalNumHydrogenAbundances), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_tableKind allocation failed')
  end if

  allocate (op_tableNameLowT (1:op_opalNumHydrogenAbundances), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_tableNameLowT allocation failed')
  end if
  allocate (op_tableNameHighT (1:op_opalNumHydrogenAbundances), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_tableNameHighT allocation failed')
  end if


  allocate (op_opalTableAbundMax (1:op_opalNumHydrogenAbundances), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_opalNumHydrogenAbundances allocation failed')
  end if


  ! LOOP # 0 - Allocation & zeroing of pointers

  allocate(op_opalAllTab(op_opalNumHydrogenAbundances,OP_OPAL_LOWT:OP_OPAL_HIGHT))
  do tempRange = OP_OPAL_LOWT,OP_OPAL_HIGHT
     do m = 1,op_opalNumHydrogenAbundances
        !nullify pointers for tables of all types
        nullify(op_opalAllTab(m,tempRange)%tg%table)
        nullify(op_opalAllTab(m,tempRange)%tg%mgTable)
        op_opalAllTab(m,tempRange)%tg%numTables = 0 ! to be modified below
        op_opalAllTab(m,tempRange)%tg%td%ntemp = 0 ! to be modified below
        op_opalAllTab(m,tempRange)%tg%td%ndens = 0 ! to be modified below
        op_opalAllTab(m,tempRange)%tg%td%nmg = 0
        nullify(op_opalAllTab(m,tempRange)%tg%td%Temperatures)
        nullify(op_opalAllTab(m,tempRange)%tg%td%Densities)
        if (.TRUE.) then
           op_opalAllTab(m,tempRange)%tg%td%isLog = op_useLogTables
        else
           op_opalAllTab(m,tempRange)%tg%td%isLog = .FALSE.
        end if
     end do
  end do

!
!    ...Initialize:   counter arrays  -> zero
!                     index arrays    -> zero
!                     max Temp arrays -> - 1.
!                     min Temp arrays -> - 1.
!
!
  op_maxTablesLowT            = 0
  op_maxTablesHighT            = 0
  op_maxTablesRO            = 0

#ifdef DEBUG_OPACITY
  print*,'op_initTabulated: starting.'
#endif
!
!
!    ...Read in the needed info from the runtime parameters and set
!       the handle arrays.
!
!
  ! LOOP # 1 

  do m = 1,op_opalNumHydrogenAbundances

     call concatStringWithInt("op_opalTableAbundMax_", m, rtpar)
     call RuntimeParameters_get(rtpar, op_opalTableAbundMax(m))

     absorptionKind = OP_CONSTCM2G

     emissionKind = OP_CONSTCM2G

     transportKind = OP_UNDEFINED

     op_tableKind(m) = "OPAL"

     call concatStringWithInt("op_opalTableLowT_", m, rtpar)
     call RuntimeParameters_get(rtpar, op_tableNameLowT(m))

     call concatStringWithInt("op_opalTableHighT_", m, rtpar)
     call RuntimeParameters_get(rtpar, op_tableNameHighT(m))


     op_absorptionKind (m) = absorptionKind
     if (absorptionKind == OP_OPAL_LOWT) op_absorptionKind (m) = OP_OPAL_LOWT
     if (absorptionKind == OP_OPAL_HIGHT) op_absorptionKind (m) = OP_OPAL_HIGHT
     if (absorptionKind == OP_TABULAR_RO) op_absorptionKind (m) = OP_TABULAR_RO

     op_emissionKind   (m) = emissionKind
     if (  emissionKind == OP_OPAL_LOWT) op_emissionKind   (m) = OP_OPAL_LOWT
     if (  emissionKind == OP_OPAL_HIGHT) op_emissionKind   (m) = OP_OPAL_HIGHT
     if (  emissionKind == OP_TABULAR_RO) op_emissionKind   (m) = OP_TABULAR_RO

     if ( transportKind == OP_OPAL_LOWT) op_transportKind  (m) = OP_OPAL_LOWT
     if ( transportKind == OP_OPAL_HIGHT) op_transportKind  (m) = OP_OPAL_HIGHT
     if ( transportKind == OP_TABULAR_RO) op_transportKind  (m) = OP_TABULAR_RO

  end do
!
!
  ! LOOP # 2 - Browse & Allocate
  !
!    ...Determine first the overall maximal dimensions needed for allocating
!       the tables.
!
!
  do m = 1,op_opalNumHydrogenAbundances

     needLowTable  =      (op_absorptionKind (m) == OP_OPAL_LOWT) &
                    .or. (op_emissionKind   (m) == OP_OPAL_LOWT) &
                    .or. (op_transportKind  (m) == OP_OPAL_LOWT) &
                    .or. (op_tableKind(m) == "OPAL" .AND. op_tableNameLowT(m) .NE. "-none-")

     needHighTable  =      (op_absorptionKind (m) == OP_OPAL_HIGHT) &
                    .or. (op_emissionKind   (m) == OP_OPAL_HIGHT) &
                    .or. (op_transportKind  (m) == OP_OPAL_HIGHT) &
                    .or. (op_tableKind(m) == "OPAL" .AND. op_tableNameHighT(m) .NE. "-none-")

     needROTable  =      (op_absorptionKind (m) == OP_TABULAR_RO) &
                    .or. (op_emissionKind   (m) == OP_TABULAR_RO) &
                    .or. (op_transportKind  (m) == OP_TABULAR_RO) &
                    .or.  needLowTable .or. needHighTable

     needTable    =       needLowTable &
                    .or.  needHighTable &
                    .or.  needROTable

     if (needLowTable) op_maxTablesLowT = op_maxTablesLowT + 1
     if (needHighTable) op_maxTablesHighT = op_maxTablesHighT + 1
     if (needROTable) op_maxTablesRO = op_maxTablesRO + 1

#ifdef DEBUG_OPACITY
     print*,'op_initTabulated: hydro abundance #, needTable =',m, needTable, needLowTable, needHighTable, needROTable
#endif
     if (needTable) then

         call op_browseTables (op_tableNameLowT (m),           &
                               op_tableNameHighT (m),           &
                               needLowTable,                      &
                               needHighTable,                      &
                               needROTable,                      &
                                            nstepsDensityLowT,     &
                                            nstepsDensityHighT,     &
                                            nstepsDensityRO,     &
                                            nstepsTemperatureLowT, &
                                            nstepsTemperatureHighT, &
                                            nstepsTemperatureRO  )

         if (needLowTable) then
            op_opalAllTab(m,OP_OPAL_LOWT)%tg%td%ntemp = nstepsTemperatureLowT
            op_opalAllTab(m,OP_OPAL_LOWT)%tg%td%ndens = nstepsDensityLowT
            nullify(op_opalAllTab(m,OP_OPAL_LOWT)%tg%td%Temperatures)
            nullify(op_opalAllTab(m,OP_OPAL_LOWT)%tg%td%Densities)
            !The following allocations will be done in op_readOpalTable:
!!$            allocate(op_opalAllTab(m,OP_OPAL_LOWT)%tg%td%Temperatures(nstepsTemperatureLowT))
!!$            allocate(op_opalAllTab(m,OP_OPAL_LOWT)%tg%td%Densities(nstepsDensityLowT))
            allocate(op_opalAllTab(m,OP_OPAL_LOWT)%tg%table)
            op_opalAllTab(m,OP_OPAL_LOWT)%tg%numTables = 1
               nullify(op_opalAllTab(m,OP_OPAL_LOWT)%tg%table%table)
               op_opalAllTab(m,OP_OPAL_LOWT)%tg%td%isLog = .TRUE.
         End if
         if (needHighTable) then
            op_opalAllTab(m,OP_OPAL_HIGHT)%tg%td%ntemp = nstepsTemperatureHighT
            op_opalAllTab(m,OP_OPAL_HIGHT)%tg%td%ndens = nstepsDensityHighT
            nullify(op_opalAllTab(m,OP_OPAL_HIGHT)%tg%td%Temperatures)
            nullify(op_opalAllTab(m,OP_OPAL_HIGHT)%tg%td%Densities)
            !The following allocations will be done in op_readOpalTable:
!!$            allocate(op_opalAllTab(m,OP_OPAL_HIGHT)%tg%td%Temperatures(nstepsTemperatureHighT))
!!$            allocate(op_opalAllTab(m,OP_OPAL_HIGHT)%tg%td%Densities(nstepsDensityHighT))
            allocate(op_opalAllTab(m,OP_OPAL_HIGHT)%tg%table)
            op_opalAllTab(m,OP_OPAL_HIGHT)%tg%numTables = 1
               nullify(op_opalAllTab(m,OP_OPAL_HIGHT)%tg%table%table)
               op_opalAllTab(m,OP_OPAL_HIGHT)%tg%td%isLog = .TRUE.
         End if
     end if

  end do
!
!
!    ...Read the tables and store the needed opacity values and all associated data
!       into the arrays.
!
!

  do m = 1,op_opalNumHydrogenAbundances

     needLowTable  =      (op_absorptionKind (m) == OP_OPAL_LOWT) &
                    .or. (op_emissionKind   (m) == OP_OPAL_LOWT) &
                    .or. (op_transportKind  (m) == OP_OPAL_LOWT) &
                    .or. (op_tableKind(m) == "OPAL" .AND. op_tableNameLowT(m) .NE. "-none-")

     needHighTable  =      (op_absorptionKind (m) == OP_OPAL_HIGHT) &
                    .or. (op_emissionKind   (m) == OP_OPAL_HIGHT) &
                    .or. (op_transportKind  (m) == OP_OPAL_HIGHT) &
                    .or. (op_tableKind(m) == "OPAL" .AND. op_tableNameHighT(m) .NE. "-none-")

     needROTable  =      (op_absorptionKind (m) == OP_TABULAR_RO) &
                    .or. (op_emissionKind   (m) == OP_TABULAR_RO) &
                    .or. (op_transportKind  (m) == OP_TABULAR_RO) &
		    .or. needLowTable .or. needHighTable

     needTable    =       needLowTable &
                    .or.  needHighTable &
                    .or.  needROTable


#ifdef DEBUG_OPACITY
     print*,'op_initTabulated: R-hydro abundance #, needTable =',m, needTable, needLowTable, needHighTable, needROTable
#endif
     if (needTable) then

         call op_readOpalTables (op_tableNameLowT (m),           &
                             op_tableNameHighT (m),           &
                             needLowTable,            &
                             needHighTable,            &
                             needROTable,              &
                             op_opalAllTab(m,:)%tg%td, &
                             op_opalAllTab(m,OP_OPAL_LOWT )%tg%table, &
                             op_opalAllTab(m,OP_OPAL_HIGHT)%tg%table)

     end if

  end do


!
!
!    ...Write out the Opacity tables and associated data (if requested).
!
!
  if (op_writeOpacityInfo) then
      call op_writeTables ()
  end if
!
!
!    ...Set initialization status.
!
!
  op_initializedTabulated = .true.
!
!
!   ...Ready! 
!
!
  return
end subroutine op_initTabulated
