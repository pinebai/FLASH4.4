!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_browseIonmixTables
!!
!! NAME
!!
!!  op_browseIonmixTables
!!
!! SYNOPSIS
!!
!!  call op_browseIonmixTables (character (in)  :: tableName (len=80),
!!                              logical   (in)  :: needPATable,
!!                              logical   (in)  :: needPETable,
!!                              logical   (in)  :: needROTable,
!!                              integer   (out) :: nstepsDensityPA,
!!                              integer   (out) :: nstepsDensityPE,
!!                              integer   (out) :: nstepsDensityRO,
!!                              integer   (out) :: nstepsTemperaturePA,
!!                              integer   (out) :: nstepsTemperaturePE,
!!                              integer   (out) :: nstepsTemperatureRO)
!!
!! DESCRIPTION
!!
!!  This routine browses through the tabulated opacities from an IONMIX datafile output in
!!  order to extract the number of steps for both the density and the temperature grid with
!!  which the IONMIX tables were generated.
!!
!! ARGUMENTS
!!
!!  tableName           : the name of the IONMIX file
!!  needPATable         : if yes, Planck Absorption Opacities are needed from the IONMIX table
!!  needPETable         : if yes, Planck   Emission Opacities are needed from the IONMIX table
!!  needROTable         : if yes,         Rosseland Opacities are needed from the IONMIX table
!!  nstepsDensityPA     : the size of the Planck Absorption density grid returned
!!  nstepsDensityPE     : the size of the Planck   Emission density grid returned
!!  nstepsDensityRO     : the size of the         Rosseland density grid returned
!!  nstepsTemperaturePA : the size of the Planck Absorption temperature grid returned
!!  nstepsTemperaturePE : the size of the Planck   Emission temperature grid returned
!!  nstepsTemperatureRO : the size of the         Rosseland temperature grid returned
!!
!!***

#include "constants.h"

subroutine op_browseIonmixTables (tableName,                        &
                                  needPATable,                      &
                                  needPETable,                      &
                                  needROTable,                      &
                                               nstepsDensityPA,     &
                                               nstepsDensityPE,     &
                                               nstepsDensityRO,     &
                                               nstepsTemperaturePA, &
                                               nstepsTemperaturePE, &
                                               nstepsTemperatureRO  )

  use Driver_interface,  ONLY : Driver_abortFlash
  use Opacity_data,  ONLY : op_globalMe

  implicit none

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

  logical :: fileExists

  integer :: fileUnit
  integer :: nstepsDensity
  integer :: nstepsTemperature
  integer :: ut_getFreeFileUnit
!
!
!   ...Check and open the IONMIX opacity file.
!
!
  inquire (file = tableName , exist = fileExists)

  if (.not.fileExists) then
     if (op_globalMe==MASTER_PE) &
          print*,'[op_browseIonmixTables] ERROR: IONMIX file not found: ',tableName 
       call Driver_abortFlash ('[op_browseIonmixTables] ERROR: no IONMIX file found')
  end if

  fileUnit = ut_getFreeFileUnit ()
  open (unit = fileUnit , file = tableName)
!
!
!   ...Read the temperature and density grids. Abort the calculation,
!      if any of the grids is not found.
!
!
  read (fileUnit,'(2I10)') nstepsTemperature , nstepsDensity

  if (nstepsTemperature <= 0) then
      call Driver_abortFlash ('[op_browseIonmixTables] ERROR: no IONMIX temperature grid found')
  end if

  if (nstepsDensity <= 0) then
      call Driver_abortFlash ('[op_browseIonmixTables] ERROR: no IONMIX density grid found')
  end if
!
!
!   ...For the IONMIX tables the grids are the same for all three Planck Absorption,
!      Planck Emission and Rosseland Transporation tables.
!
!
  nstepsDensityPA = nstepsDensity
  nstepsDensityPE = nstepsDensity
  nstepsDensityRO = nstepsDensity

  nstepsTemperaturePA = nstepsTemperature
  nstepsTemperaturePE = nstepsTemperature
  nstepsTemperatureRO = nstepsTemperature
!
!
!   ...Close the IONMIX file.
!
!
  close (fileUnit)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_browseIonmixTables
