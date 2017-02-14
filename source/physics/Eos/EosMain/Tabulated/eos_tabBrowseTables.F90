!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabBrowseTables
!!
!! NAME
!!
!!  eos_tabBrowseTables
!!
!! SYNOPSIS
!!
!!  call eos_tabBrowseTables (character (in)  :: tableKind (len=80),
!!                        character (in)  :: tableName (len=80),
!!                      character (in) :: groupName (len=80),
!!                        logical   (in)  :: needZFTable,
!!                        logical   (in)  :: needENTable,
!!                        logical   (in)  :: needHCTable,
!!                        integer   (out) :: nstepsDensityZF,
!!                        integer   (out) :: nstepsDensityEN,
!!                        integer   (out) :: nstepsDensityHC,
!!                        integer   (out) :: nstepsTemperatureZF,
!!                        integer   (out) :: nstepsTemperatureEN,
!!                        integer   (out) :: nstepsTemperatureHC)
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

subroutine eos_tabBrowseTables (tableKind,                   &
                            tableName,                   &
                            groupName,   &
                            needZFTable,                 &
                            needENTable,                 &
                            needHCTable,                 &
                            needEntrTable,               &
                                    nstepsDensityZF,     &
                                    nstepsDensityEN,     &
                                    nstepsDensityHC,     &
                                    nstepsDensityEntr,     &
                                    nstepsTemperatureZF, &
                                    nstepsTemperatureEN, &
                                    nstepsTemperatureHC, &
                                    nstepsTemperatureEntr)

  use Driver_interface,   ONLY : Driver_abortFlash
  use eos_tabInterface,   ONLY : eos_tabBrowseIonmixTables
  use Eos_data,           ONLY : eos_meshMe

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
  character (len=80), intent (in)  :: tableKind
  character (len=80), intent (in)  :: tableName
  character (len=80), intent (in)  :: groupName  !DEV: UNUSED

  nstepsDensityEntr = 0
  nstepsTemperatureEntr = 0
!
!
!   ...Call the appropriate routine.
!
!
  if (tableKind == 'IONMIX') then

      call eos_tabBrowseIonmixTables (tableName,                   &
                                  needZFTable,                 &
                                  needENTable,                 &
                                  needHCTable,                 &
                                          nstepsDensityZF,     &
                                          nstepsDensityEN,     &
                                          nstepsDensityHC,     &
                                          nstepsTemperatureZF, &
                                          nstepsTemperatureEN, &
                                          nstepsTemperatureHC  )

   else if (tableKind=='IONMIX4' .OR. tableKind=='IONMIX6') then

      call eos_tabBrowseIonmix4Tables(tableName,                   &
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

  else
     if (eos_meshMe==MASTER_PE) then
        print*, '[eos_tabBrowseTables] ERROR: Did not recognize EOS table kind: ', &
             tableKind
     end if
     call Driver_abortFlash ('[eos_tabBrowseTables] ERROR: EOS table kind not recognized')
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine eos_tabBrowseTables
