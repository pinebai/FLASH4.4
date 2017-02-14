!!****if* source/physics/materialProperties/Opacity/OpacityMain/OPAL/op_readTables
!!
!! NAME
!!
!!  op_readTables
!!
!! SYNOPSIS
!!
!!  call op_readTables (character (in) :: tableKind (len=80),
!!                      character (in) :: tableName (len=80),
!!                      logical   (in) :: needLowTable,
!!                      logical   (in) :: needHighTable,
!!                      logical   (in) :: needROTable,
!!                      integer   (in) :: indexLOWT,
!!                      integer   (in) :: indexHIGHT,
!!                      integer   (in) :: indexRO)
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
!!  tableKind   : the kind of tabulated Opacity file where data is going to be read
!!  tableName   : the name of tabulated Opacity file where data is going to be read
!!  needLowTable : if yes, Planck Absorption Opacities are needed from the Opacity table
!!  needHighTable : if yes, Planck   Emission Opacities are needed from the Opacity table
!!  needROTable : if yes,         Rosseland Opacities are needed from the Opacity table
!!  indexLOWT     : table counting index where Planck Absorption Opacities will be placed
!!  indexHIGHT     : table counting index where Planck   Emission Opacities will be placed
!!  indexRO     : table counting index where Planck  Transport Opacities will be placed
!!
!!***
subroutine op_readTables (tableKind,   &
                          tableName,   &
                          needLowTable, &
                          needHighTable, &
                          needROTable, &
                          indexLOWT,     &
                          indexHIGHT,     &
                          indexRO      )

  use Driver_interface,   ONLY : Driver_abortFlash
  use op_interface,       ONLY : op_readIonmixTables,  &
                                 op_readIonmix4Tables
  use op_opalInterface,   ONLY : op_readOpalTable

  implicit none

  character (len=80), intent (in) :: tableKind
  character (len=80), intent (in) :: tableName
  logical,            intent (in) :: needLowTable
  logical,            intent (in) :: needHighTable
  logical,            intent (in) :: needROTable
  integer,            intent (in) :: indexLowT
  integer,            intent (in) :: indexHighT
  integer,            intent (in) :: indexRO
!
!
!   ...Call the appropriate routine.
!
!
  if (tableKind == "IONMIX" .or. tableKind == "ionmix") then

      call op_readIonmixTables (tableName,   &
                                needLowTable, &
                                needHighTable, &
                                needROTable, &
                                indexLOWT,     &
                                indexHIGHT,     &
                                indexRO      )

  elseif (tableKind == "IONMIX4" .or. tableKind == "ionmix4") then

      call op_readIonmix4Tables (tableName,   &
                                 needLowTable, &
                                 needHighTable, &
                                 needROTable, &
                                 indexLOWT,     &
                                 indexHIGHT,     &
                                 indexRO      )

  elseif (tableKind == "OPAL" .or. tableKind == "opal") then

      call Driver_abortFlash ('[op_readTables] ERROR: Reading OPAL pacities not implemented here, use op_readOpalTables.')

  elseif (tableKind == "PROPACEOS" .or. tableKind == "propaceos") then

      call op_readPropaceosTables (tableName,   &
                                 needLowTable, &
                                 needHighTable, &
                                 needROTable, &
                                 indexLOWT,     &
                                 indexHIGHT,     &
                                 indexRO      )

   else
      call Driver_abortFlash ('[op_readTables] ERROR: Opacity table kind not recognized')
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_readTables
