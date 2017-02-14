!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabBrowseIonmix4Tables
!!
!! NAME
!!
!!  eos_tabBrowseIonmix4Tables
!!
!! SYNOPSIS
!!
!!  call eos_tabBrowseIonmix4Tables (character (in)  :: tableName (len=80),
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
!!  This routine browses through the tabulated opacities from an IONMIX4 datafile output in
!!  order to extract the number of steps for both the density and the temperature grid with
!!  which the IONMIX4 tables were generated.
!!
!!  THis routine can be be used on IONMIX4 and IONMIX6 tables.
!!
!! ARGUMENTS
!!
!!  tableName           : the name of the IONMIX4 file
!!  needZFTable         : if yes, average ionization data are needed from the IONMIX4 table
!!  needENTable         : if yes, internal energy data are needed from the IONMIX4 table
!!  needHCTable         : if yes,        specific heat data are needed from the IONMIX4 table
!!  nstepsDensityZF     : the size of the average ionization density grid returned
!!  nstepsDensityEN     : the size of the internal energy density grid returned
!!  nstepsDensityHC     : the size of the        specific heat density grid returned
!!  nstepsTemperatureZF : the size of the average ionization temperature grid returned
!!  nstepsTemperatureEN : the size of the internal energy temperature grid returned
!!  nstepsTemperatureHC : the size of the        specific heat temperature grid returned
!!
!!
!! NOTES
!!
!!  Yup, up to now the code is indentical to that for IONMIX1 (aka "IONMIX") files.
!!  Maybe that will change.
!!
!!***

#include "constants.h"

subroutine eos_tabBrowseIonmix4Tables (tableName,                   &
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

  logical :: fileExists

  integer :: fileUnit
  integer :: nstepsDensity
  integer :: nstepsTemperature
  integer :: ut_getFreeFileUnit
!
!
!   ...Check and open the IONMIX4 file.
!   This should be a Two-temperature, Custom Grid, cnrdeos File, to which
!   the following descrption applies:
!   IONMIX has been extended by Duc Cao to include three new features:
!
!   1. Two-temperature (Ion/Electron) EOS data
!   2. Pressure information
!   3. Temperature/Density points that are individually specified - so that you do not have to use points that are evenly spaced logarithmically
!
!   These new files will have the extension .cn4 instead of .cnr.
!
!
  inquire (file = tableName , exist = fileExists)

  if (.not.fileExists) then
     if (eos_meshMe==MASTER_PE) &
          print*,'[eos_tabBrowseIonmix4Tables] ERROR: IONMIX4 file not found: ',tableName
     call Driver_abortFlash ('[eos_tabBrowseIonmix4Tables] ERROR: no IONMIX4 file found')
  else if (eos_meshMe==MASTER_PE) then
     print*,'[eos_tabBrowseIonmix4Tables] IONMIX4 file found: ',tableName
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
      call Driver_abortFlash ('[eos_tabBrowseIonmix4Tables] ERROR: no IONMIX4 temperature grid found')
  end if

  if (nstepsDensity <= 0) then
      call Driver_abortFlash ('[eos_tabBrowseIonmix4Tables] ERROR: no IONMIX4 density grid found')
  end if
!
!
!   ...For the IONMIX4 tables the grids are the same for all three tables, namely, average ionization,
!      internal energy, and specific heat tables.
!
!
  nstepsDensityZF = nstepsDensity
  nstepsDensityEN = nstepsDensity
  nstepsDensityHC = nstepsDensity
  nstepsDensityEntr = nstepsDensity

  nstepsTemperatureZF = nstepsTemperature
  nstepsTemperatureEN = nstepsTemperature
  nstepsTemperatureHC = nstepsTemperature
  nstepsTemperatureEntr = nstepsTemperature
!
!
!   ...Close the IONMIX4 file.
!
!
  close (fileUnit)
!
!
!   ...Ready! 
!
!
  return
end subroutine eos_tabBrowseIonmix4Tables
