!!****if* source/physics/materialProperties/Opacity/OpacityMain/OPAL/op_browseTables
!!
!! NAME
!!
!!  op_browseTables
!!
!! SYNOPSIS
!!
!!  call op_browseTables (character (in)  :: tableKind (len=80),
!!                        character (in)  :: tableName (len=80),
!!                        logical   (in)  :: needLowTable,
!!                        logical   (in)  :: needHighTable,
!!                        logical   (in)  :: needROTable,
!!                        integer   (out) :: nstepsDensityLowT,
!!                        integer   (out) :: nstepsDensityHighT,
!!                        integer   (out) :: nstepsDensityRO,
!!                        integer   (out) :: nstepsTemperatureLowT,
!!                        integer   (out) :: nstepsTemperatureHighT,
!!                        integer   (out) :: nstepsTemperatureRO)
!!
!! DESCRIPTION
!!
!!  This operation browses through specific tables returning the number of
!!  steps for both the density and the temperature. This is needed for establishing
!!  the maximum table allocation sizes. In this implementation only OPAL tables can be browsed.
!!
!! ARGUMENTS
!!
!!  tableKind           : the kind of tabulated Opacity file where data is going to be read
!!  tableName           : the name of tabulated Opacity file where data is going to be read
!!  needLowTable         : if yes, Low Temperature Opacities are needed from the Opacity table
!!  needHighTable         : if yes, high temperature Opacities are needed from the Opacity table
!!  needROTable         : if yes,         Rosseland Opacities are needed from the Opacity table
!!  nstepsDensityLowT     : the size of the Low Temperature density grid returned
!!  nstepsDensityHighT     : the size of the high temperature density grid returned
!!  nstepsDensityRO     : the size of the         Rosseland density grid returned
!!  nstepsTemperatureLowT : the size of the Low Temperature temperature grid returned
!!  nstepsTemperatureHighT : the size of the high temperature temperature grid returned
!!  nstepsTemperatureRO : the size of the         Rosseland temperature grid returned
!!
!!***
subroutine op_browseTables (tableNameLowT,                     &
                            tableNameHighT,                    &
                            needLowTable,                      &
                            needHighTable,                      &
                            needROTable,                      &
                                         nstepsDensityLowT,     &
                                         nstepsDensityHighT,     &
                                         nstepsDensityRO,     &
                                         nstepsTemperatureLowT, &
                                         nstepsTemperatureHighT, &
                                         nstepsTemperatureRO  )

  use Driver_interface,   ONLY : Driver_abortFlash
  use op_interface,       ONLY : op_browseOpalTable

  implicit none

  character (len=80), intent (in)  :: tableNameLowT
  character (len=80), intent (in)  :: tableNameHighT
  logical,            intent (in)  :: needLowTable
  logical,            intent (in)  :: needHighTable
  logical,            intent (in)  :: needROTable
  integer,            intent (out) :: nstepsDensityLowT
  integer,            intent (out) :: nstepsDensityHighT
  integer,            intent (out) :: nstepsDensityRO
  integer,            intent (out) :: nstepsTemperatureLowT
  integer,            intent (out) :: nstepsTemperatureHighT
  integer,            intent (out) :: nstepsTemperatureRO

  nstepsDensityLowT = 0
  nstepsDensityHighT = 0
  nstepsDensityRO = 0
  nstepsTemperatureLowT = 0
  nstepsTemperatureHighT = 0
  nstepsTemperatureRO = 0
!
!
!   ...Call the appropriate routine.
!
!
  if (needROTable) then
    if (needLowTable ) then
      call op_browseOpalTable (tableNameLowT,                        &
                                               nstepsDensityLowT,     &
                                               nstepsTemperatureLowT)
    end if

    if (needHighTable ) then
      call op_browseOpalTable (tableNameHighT,                        &
                                               nstepsDensityHighT,     &
                                               nstepsTemperatureHighT)
    end if
    nstepsDensityRO = max(nstepsDensityLowT, nstepsDensityHighT)
    nstepsTemperatureRO = max(nstepsTemperatureLowT, nstepsTemperatureHighT)
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_browseTables
