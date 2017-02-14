!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_writeSpeciesTables
!!
!! NAME
!!
!!  op_writeSpeciesTables
!!
!! SYNOPSIS
!!
!!  call op_writeSpeciesTables (integer (in) :: fileUnit,
!!                              integer (in) :: species,
!!                              logical (in) :: needPATable,
!!                              logical (in) :: needPETable,
!!                              logical (in) :: needROTable)
!!
!! DESCRIPTION
!!
!!  Prints out the tabulated Opacities for the current species. Only those tables
!!  are printed which were actually generated for the current species.
!!
!! ARGUMENTS
!!
!!  fileUnit    : unit # for the output file
!!  species     : the species index
!!  needPATable : if yes, Planck Absorption Opacities were stored for the current species
!!  needPETable : if yes, Planck   Emission Opacities were stored for the current species
!!  needROTable : if yes,         Rosseland Opacities were stored for the current species
!!
!!***
subroutine op_writeSpeciesTables (fileUnit,    &
                                  species,     &
                                  needPATable, &
                                  needPETable, &
                                  needROTable  )

  use op_interface, ONLY : op_writeSpeciesPATable, &
                           op_writeSpeciesPETable, &
                           op_writeSpeciesROTable

  implicit none

  logical, intent (in) :: needPATable
  logical, intent (in) :: needPETable
  logical, intent (in) :: needROTable
  integer, intent (in) :: fileUnit
  integer, intent (in) :: species
!
!
!   ...Print only the relevant opacities from the tables.
!
!
  if (needPATable) then

      write (fileUnit,*)
      write (fileUnit,*) ' ------- PLANCK ABSORPTION OPACITIES -------'
      write (fileUnit,*)

      call op_writeSpeciesPATable (fileUnit,species)
  end if

  if (needPETable) then

      write (fileUnit,*)
      write (fileUnit,*) ' ------- PLANCK EMISSION OPACITIES -------'
      write (fileUnit,*)

      call op_writeSpeciesPETable (fileUnit,species)
  end if

  if (needROTable) then

      write (fileUnit,*)
      write (fileUnit,*) ' ------- ROSSELAND OPACITIES -------'
      write (fileUnit,*)

      call op_writeSpeciesROTable (fileUnit,species)
  end if
!
!
!   ...Ready!
!
!
  return
end subroutine op_writeSpeciesTables
