!!****if* source/physics/materialProperties/Opacity/localAPI/op_readTables
!!
!! NAME
!!
!!  op_readTables
!!
!! SYNOPSIS
!!
!!  call op_readTables (character (in) :: tableKind (len=80),
!!                      character (in) :: tableName (len=80),
!!                      logical   (in) :: needPATable,
!!                      logical   (in) :: needPETable,
!!                      logical   (in) :: needROTable,
!!                      integer   (in) :: indexPA,
!!                      integer   (in) :: indexPE,
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
!!  needPATable : if yes, Planck Absorption Opacities are needed from the Opacity table
!!  needPETable : if yes, Planck   Emission Opacities are needed from the Opacity table
!!  needROTable : if yes,         Rosseland Opacities are needed from the Opacity table
!!  indexPA     : table counting index where Planck Absorption Opacities will be placed
!!  indexPE     : table counting index where Planck   Emission Opacities will be placed
!!  indexRO     : table counting index where Planck  Transport Opacities will be placed
!!
!!***
subroutine op_readTables (tableKind,   &
                          tableName,   &
                          needPATable, &
                          needPETable, &
                          needROTable, &
                          indexPA,     &
                          indexPE,     &
                          indexRO      )
  implicit none

  character (len=80), intent (in) :: tableKind
  character (len=80), intent (in) :: tableName
  logical,            intent (in) :: needPATable
  logical,            intent (in) :: needPETable
  logical,            intent (in) :: needROTable
  integer,            intent (in) :: indexPA
  integer,            intent (in) :: indexPE
  integer,            intent (in) :: indexRO

  return
end subroutine op_readTables
