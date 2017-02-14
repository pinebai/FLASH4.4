!!****if* source/physics/materialProperties/Opacity/localAPI/op_browseTables
!!
!! NAME
!!
!!  op_browseTables
!!
!! SYNOPSIS
!!
!!  call op_browseTables (character (in)  :: tableKind (len=80),
!!                        character (in)  :: tableName (len=80),
!!                        logical   (in)  :: needPATable,
!!                        logical   (in)  :: needPETable,
!!                        logical   (in)  :: needROTable,
!!                        integer   (out) :: nstepsDensityPA,
!!                        integer   (out) :: nstepsDensityPE,
!!                        integer   (out) :: nstepsDensityRO,
!!                        integer   (out) :: nstepsTemperaturePA,
!!                        integer   (out) :: nstepsTemperaturePE,
!!                        integer   (out) :: nstepsTemperatureRO)
!!
!! DESCRIPTION
!!
!!  This operation browses through specific tables returning the number of
!!  steps for both the density and the temperature. This is needed for establishing
!!  the maximum table allocation sizes. Currently only the IONMIX tables can be browsed.
!!
!! ARGUMENTS
!!
!!  tableKind           : the kind of tabulated Opacity file where data is going to be read
!!  tableName           : the name of tabulated Opacity file where data is going to be read
!!  needPATable         : if yes, Planck Absorption Opacities are needed from the Opacity table
!!  needPETable         : if yes, Planck   Emission Opacities are needed from the Opacity table
!!  needROTable         : if yes,         Rosseland Opacities are needed from the Opacity table
!!  nstepsDensityPA     : the size of the Planck Absorption density grid returned
!!  nstepsDensityPE     : the size of the Planck   Emission density grid returned
!!  nstepsDensityRO     : the size of the         Rosseland density grid returned
!!  nstepsTemperaturePA : the size of the Planck Absorption temperature grid returned
!!  nstepsTemperaturePE : the size of the Planck   Emission temperature grid returned
!!  nstepsTemperatureRO : the size of the         Rosseland temperature grid returned
!!
!!***
subroutine op_browseTables (tableKind,                       &
                            tableName,                       &
                            needPATable,                     &
                            needPETable,                     &
                            needROTable,                     &
                                        nstepsDensityPA,     &
                                        nstepsDensityPE,     &
                                        nstepsDensityRO,     &
                                        nstepsTemperaturePA, &
                                        nstepsTemperaturePE, &
                                        nstepsTemperatureRO  )
  implicit none

  character (len=80), intent (in)  :: tableKind
  character (len=80), intent (in)  :: tableName
  logical,            intent (in)  :: needPATable
  logical,            intent (in)  :: needPETable
  logical,            intent (in)  :: needROTable
  integer,            intent (out) :: nstepsDensityPA
  integer,            intent (out) :: nstepsDensityPE
  integer,            intent (out) :: nstepsDensityRO
  integer,            intent (out) :: nstepsTemperaturePA
  integer,            intent (out) :: nstepsTemperaturePE
  integer,            intent (out) :: nstepsTemperatureRO

  nstepsDensityPA     = 0
  nstepsDensityPE     = 0
  nstepsDensityRO     = 0
  nstepsTemperaturePA = 0
  nstepsTemperaturePE = 0
  nstepsTemperatureRO = 0

  return
end subroutine op_browseTables
