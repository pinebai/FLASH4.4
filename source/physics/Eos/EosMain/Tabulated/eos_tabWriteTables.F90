!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabWriteTables
!!
!! NAME
!!
!!  eos_tabWriteTables
!!
!! SYNOPSIS
!!
!!  call eos_tabWriteTables (integer   (in) :: fileUnit,
!!                       character (in) :: header)
!!
!! DESCRIPTION
!!
!!  Prints out the tabulated Opacities for each species. Only those tables
!!  are printed which were actually generated for each species.
!!
!! ARGUMENTS
!!
!!  fileUnit = unit # for the output file
!!  header   = header of the printout
!!
!!***
subroutine eos_tabWriteTables (fileUnit,header)

  use eos_tabInterface, ONLY : eos_tabWriteSpeciesTables

  use eos_tabData, ONLY : eos_tabTotalNumSpecies,   &
                           eos_tabIonizationKind, &
                           eos_tabIntEnergyKind,   &
                           eos_tabHeatCpKind, &
                           eos_tableKind

  implicit none

#include "Flash.h"
#include "Eos.h"

  character*(*), intent (in) :: header
  integer,       intent (in) :: fileUnit

  logical :: needTable
  logical :: needZFTable
  logical :: needENTable
  logical :: needHCTable
  logical :: needPRTable
  logical :: needEntrTable
  integer :: species
!
!
!   ...Print out the title. 
!
!
  write (fileUnit,*)
  write (fileUnit,*) header
  write (fileUnit,*)
!
!
!   ...Loop over all species. 
!
!
  do species = 1,eos_tabTotalNumSpecies

     needZFTable  =  ANY((eos_tabIonizationKind (:,species) == EOS_TABULAR_Z) &
                    .or. (eos_tabIntEnergyKind   (:,species) == EOS_TABULAR_Z) &
                    .or. (eos_tabHeatCpKind  (:,species) == EOS_TABULAR_Z))

     needENTable  =  ANY((eos_tabIonizationKind (:,species) == EOS_TABULAR_E) &
                    .or. (eos_tabIntEnergyKind   (:,species) == EOS_TABULAR_E) &
                    .or. (eos_tabHeatCpKind  (:,species) == EOS_TABULAR_E))

     needHCTable  =  ANY((eos_tabIonizationKind (:,species) == EOS_TABULAR_C) &
                    .or. (eos_tabIntEnergyKind   (:,species) == EOS_TABULAR_C) &
                    .or. (eos_tabHeatCpKind  (:,species) == EOS_TABULAR_C))

     needPRTable  =  ANY((eos_tabIonizationKind (:,species) == EOS_TABULAR_P) &
                    .or. (eos_tabIntEnergyKind   (:,species) == EOS_TABULAR_P) &
                    .or. (eos_tabHeatCpKind  (:,species) == EOS_TABULAR_P))

     needEntrTable = (eos_tableKind(species)=='IONMIX6')

     needTable    =       needZFTable &
                    .or.  needENTable &
                    .or.  needHCTable &
                    .or.  needPRTable .or. needEntrTable

     if (needTable) then

         write (fileUnit,*)
         write (fileUnit,*) ' ------- Species # ',species,' -------'
         write (fileUnit,*)

         call eos_tabWriteSpeciesTables (fileUnit,    &
                                     species,     &
                                     needZFTable, &
                                     needENTable, &
                                     needHCTable, &
                                     needPRTable, &
                                     needEntrTable)
     end if

  end do
!
!
!    ...Ready!
!
!
  return
end subroutine eos_tabWriteTables
