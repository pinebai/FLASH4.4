!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabWriteSpeciesTables
!!
!! NAME
!!
!!  eos_tabWriteSpeciesTables
!!
!! SYNOPSIS
!!
!!  call eos_tabWriteSpeciesTables (integer (in) :: fileUnit,
!!                              integer (in) :: species,
!!                              logical (in) :: needZFTable,
!!                              logical (in) :: needENTable,
!!                              logical (in) :: needHCTable,
!!                              logical (in) :: needPRTable)
!!
!! DESCRIPTION
!!
!!  Prints out the tabulated data for the current species. Only those tables
!!  are printed which were actually generated for the current species.
!!
!! ARGUMENTS
!!
!!  fileUnit    : unit # for the output file
!!  species     : the species index
!!  needZFTable : if yes, average ionization data were stored for the current species
!!  needENTable : if yes, internal energy data were stored for the current species
!!  needHCTable : if yes,        specific heat data were stored for the current species
!!  needPRTable : if yes,        pressure data were stored for the current species
!!
!!***
subroutine eos_tabWriteSpeciesTables (fileUnit,    &
                                  species,     &
                                  needZFTable, &
                                  needENTable, &
                                  needHCTable, &
                                  needPRTable, &
                                  needEntrTable  )
  use eos_tabInterface, ONLY: eos_tabWriteSpeciesZFTable, &
                              eos_tabWriteSpeciesENTable, &
                              eos_tabWriteSpeciesHCTable, &
                              eos_tabWriteSpeciesPRTable, &
                              eos_tabWriteSpeciesEntrTable

  implicit none

  logical, intent (in) :: needZFTable
  logical, intent (in) :: needENTable
  logical, intent (in) :: needHCTable
  logical, intent (in) :: needPRTable
  logical, intent (in) :: needEntrTable
  integer, intent (in) :: fileUnit
  integer, intent (in) :: species
!
!
!   ...Print only the relevant data from the tables.
!
!
  if (needZFTable) then

      write (fileUnit,*)
      write (fileUnit,*) ' ------- IONIZATION DATA -------'
      write (fileUnit,*)

      call eos_tabWriteSpeciesZFTable (fileUnit,species)
  end if

  if (needENTable) then

      write (fileUnit,*)
      write (fileUnit,*) ' ------- INTERNAL ENERGY DATA [erg/g] -------'
      write (fileUnit,*)

      call eos_tabWriteSpeciesENTable (fileUnit,species)
  end if

  if (needHCTable) then

      write (fileUnit,*)
      write (fileUnit,*) ' ------- HEAT CAPACITY DATA [erg/g/K] -------'
      write (fileUnit,*)

      call eos_tabWriteSpeciesHCTable (fileUnit,species)
  end if

  if (needPRTable) then

      write (fileUnit,*)
      write (fileUnit,*) ' ------- PRESSURE DATA [erg/cm^3] -------'
      write (fileUnit,*)

      call eos_tabWriteSpeciesPRTable (fileUnit,species)
  end if

  if (needEntrTable) then

      write (fileUnit,*)
      write (fileUnit,*) ' ------- ENTROPY DATA [erg/K/g ?] -------'
      write (fileUnit,*)

      call eos_tabWriteSpeciesEntrTable (fileUnit,species)
  end if
!
!
!   ...Ready!
!
!
  return
end subroutine eos_tabWriteSpeciesTables
