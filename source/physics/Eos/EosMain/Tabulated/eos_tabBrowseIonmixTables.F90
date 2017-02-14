!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabBrowseIonmixTables
!!
!! NAME
!!
!!  eos_tabBrowseIonmixTables
!!
!! SYNOPSIS
!!
!!  call eos_tabBrowseIonmixTables (character (in)  :: tableName (len=80),
!!                              logical   (in)  :: needZFTable,
!!                              logical   (in)  :: needENTable,
!!                              logical   (in)  :: needHCTable,
!!                              integer   (out) :: nstepsDensityZF,
!!                              integer   (out) :: nstepsDensityEN,
!!                              integer   (out) :: nstepsDensityHC,
!!                              integer   (out) :: nstepsTemperatureZF,
!!                              integer   (out) :: nstepsTemperatureEN,
!!                              integer   (out) :: nstepsTemperatureHC)
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
!!  needZFTable         : if yes, average ionization data are needed from the IONMIX table
!!  needENTable         : if yes, internal energy data are needed from the IONMIX table
!!  needHCTable         : if yes,        specific heat data are needed from the IONMIX table
!!  nstepsDensityZF     : the size of the average ionization density grid returned
!!  nstepsDensityEN     : the size of the internal energy density grid returned
!!  nstepsDensityHC     : the size of the        specific heat density grid returned
!!  nstepsTemperatureZF : the size of the average ionization temperature grid returned
!!  nstepsTemperatureEN : the size of the internal energy temperature grid returned
!!  nstepsTemperatureHC : the size of the        specific heat temperature grid returned
!!
!!***

#include "constants.h"

subroutine eos_tabBrowseIonmixTables (tableName,                   &
                                  needZFTable,                 &
                                  needENTable,                 &
                                  needHCTable,                 &
                                          nstepsDensityZF,     &
                                          nstepsDensityEN,     &
                                          nstepsDensityHC,     &
                                          nstepsTemperatureZF, &
                                          nstepsTemperatureEN, &
                                          nstepsTemperatureHC  )

  use Driver_interface,  ONLY : Driver_abortFlash
  use Eos_data,  ONLY : eos_meshMe

  implicit none

  logical,            intent (in)  :: needZFTable
  logical,            intent (in)  :: needENTable
  logical,            intent (in)  :: needHCTable
  integer,            intent (out) :: nstepsDensityZF
  integer,            intent (out) :: nstepsDensityEN
  integer,            intent (out) :: nstepsDensityHC
  integer,            intent (out) :: nstepsTemperatureZF
  integer,            intent (out) :: nstepsTemperatureEN
  integer,            intent (out) :: nstepsTemperatureHC
  character (len=80), intent (in)  :: tableName

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
     if (eos_meshMe==MASTER_PE) &
          print*,'[eos_tabBrowseIonmixTables] ERROR: IONMIX file not found: ',tableName
     call Driver_abortFlash ('[eos_tabBrowseIonmixTables] ERROR: no IONMIX file found')
  else if (eos_meshMe==MASTER_PE) then
     print*,'[eos_tabBrowseIonmixTables] IONMIX file found: ',tableName
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
      call Driver_abortFlash ('[eos_tabBrowseIonmixTables] ERROR: no IONMIX temperature grid found')
  end if

  if (nstepsDensity <= 0) then
      call Driver_abortFlash ('[eos_tabBrowseIonmixTables] ERROR: no IONMIX density grid found')
  end if
!
!
!   ...For the IONMIX tables the grids are the same for all three tables, namely, average ionization,
!      internal energy, and specific heat tables.
!
!
  nstepsDensityZF = nstepsDensity
  nstepsDensityEN = nstepsDensity
  nstepsDensityHC = nstepsDensity

  nstepsTemperatureZF = nstepsTemperature
  nstepsTemperatureEN = nstepsTemperature
  nstepsTemperatureHC = nstepsTemperature
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
end subroutine eos_tabBrowseIonmixTables
