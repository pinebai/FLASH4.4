!!****if* source/physics/Eos/EosMain/Tabulated/Hdf5TableRead/eos_tabBrowseOpacplotTables
!!
!! NAME
!!
!!  eos_tabBrowseOpacplotTables
!!
!! SYNOPSIS
!!
!!  call eos_tabBrowseOpacplotTables (character (in)  :: tableName (len=80),
!!                                    character (in)  :: groupName (len=80),
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
!!  This routine browses through the tabulated opacities from an OPACPLOT datafile output in
!!  order to extract the number of steps for both the density and the temperature grid with
!!  which the OPACPLOT tables were generated.
!!
!!  THis routine can be be used on OPACPLOT and IONMIX6 tables.
!!
!! ARGUMENTS
!!
!!  tableName           : the name of the OPACPLOT file
!!  groupName           : the name of the hdf5 group, often material name, in Opacplot file
!!  needZFTable         : if yes, average ionization data are needed from the OPACPLOT table
!!  needENTable         : if yes, internal energy data are needed from the OPACPLOT table
!!  needHCTable         : if yes,        specific heat data are needed from the OPACPLOT table
!!  nstepsDensityZF     : the size of the average ionization density grid returned
!!  nstepsDensityEN     : the size of the internal energy density grid returned
!!  nstepsDensityHC     : the size of the        specific heat density grid returned
!!  nstepsTemperatureZF : the size of the average ionization temperature grid returned
!!  nstepsTemperatureEN : the size of the internal energy temperature grid returned
!!  nstepsTemperatureHC : the size of the        specific heat temperature grid returned
!!
!!
!!***

#include "constants.h"

subroutine eos_tabBrowseOpacplotTables (tableName,             &
                                    groupName,         &
                                    needZFTable,                 &
                                    needENTable,                 &
                                    needHCTable,                 &
                                    needEntrTable,               &
                                          nstepsDensityZF,     &
                                          nstepsDensityEN,     &
                                          nstepsDensityHC,     &
                                          nstepsDensityEntr,   &
                                          nstepsTemperatureZF, &
                                          nstepsTemperatureEN, &
                                          nstepsTemperatureHC, &
                                          nstepsTemperatureEntr)

  use Driver_interface,  ONLY : Driver_abortFlash
  use Eos_data,  ONLY : eos_meshMe

  implicit none

  logical,            intent (in)  :: needZFTable
  logical,            intent (in)  :: needENTable
  logical,            intent (in)  :: needHCTable
  logical,            intent (in)  :: needEntrTable
  integer,            intent (out) :: nstepsDensityZF
  integer,            intent (out) :: nstepsDensityEN
  integer,            intent (out) :: nstepsDensityHC
  integer,            intent (out) :: nstepsDensityEntr
  integer,            intent (out) :: nstepsTemperatureZF
  integer,            intent (out) :: nstepsTemperatureEN
  integer,            intent (out) :: nstepsTemperatureHC
  integer,            intent (out) :: nstepsTemperatureEntr
  character (len=80), intent (in)  :: tableName
  character (len=80), intent (in)  :: groupName


  character (len=81)  :: nullTermTableName, nullTermGroupName

  logical :: fileExists

  integer :: fileUnit
  integer :: nstepsDensity
  integer :: nstepsTemperature
  integer :: ut_getFreeFileUnit

!
!   ...Check and open the OPACPLOT file.
!   This should be a Two-temperature, Custom Grid, cnrdeos File, to which
!   the following descrption applies:
!   IONMIX has been extended by Duc Cao to include three new features:
!
!   1. Two-temperature (Ion/Electron) EOS data
!   2. Pressure information
!   3. Temperature/Density points that are individually specified - so that you do not have to use points that are evenly spaced logarithmically

  inquire (file = tableName , exist = fileExists)

  if (.not.fileExists) then
     if (eos_meshMe==MASTER_PE) &
          print*,'[eos_tabBrowseOpacplotTables] ERROR: OPACPLOT file not found: ',tableName
     call Driver_abortFlash ('[eos_tabBrowseOpacplotTables] ERROR: no OPACPLOT file found')
  else if (eos_meshMe==MASTER_PE) then
     print*,'[eos_tabBrowseOpacplotTables] OPACPLOT file found: ',tableName
  end if


!   ...Read the temperature and density grids. Abort the calculation,
!      if any of the grids is not found.

  nullTermTableName = trim(tableName)//char(0)
  nullTermGroupName = trim(groupName)//char(0)
  call eos_tabBrowsehdf5(nullTermTableName, nullTermGroupName, nstepsTemperature, nstepsDensity)

  if (nstepsTemperature <= 0) then
      call Driver_abortFlash ('[eos_tabBrowseOpacplotTables] ERROR: no OPACPLOT temperature grid found')
  end if

  if (nstepsDensity <= 0) then
      call Driver_abortFlash ('[eos_tabBrowseOpacplotTables] ERROR: no OPACPLOT density grid found')
  end if

  nstepsDensityZF = nstepsDensity
  nstepsDensityEN = nstepsDensity
  nstepsDensityHC = nstepsDensity
  nstepsDensityEntr = nstepsDensity

  nstepsTemperatureZF = nstepsTemperature
  nstepsTemperatureEN = nstepsTemperature
  nstepsTemperatureHC = nstepsTemperature
  nstepsTemperatureEntr = nstepsTemperature

  return
end subroutine eos_tabBrowseOpacplotTables
