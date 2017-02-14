!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_initTabulated
!!
!! NAME
!!
!!  op_initTabulated
!!
!! SYNOPSIS
!!
!!  call op_initTabulated ()
!!
!! DESCRIPTION
!!
!!  Inititalizes the tabulated section of the opacity unit.
!!
!! ARGUMENTS
!!
!!***
subroutine op_initTabulated ()

  use Driver_interface,            ONLY : Driver_abortFlash

  use Opacity_data,                ONLY : op_totalSpecies,              &
                                          op_nEnergyGroups,             &
                                          op_absorptionKind,            &
                                          op_emissionKind,              &
                                          op_transportKind,             &
                                          op_writeOpacityInfo

  use op_tabulatedData,            ONLY : op_initializedTabulated,      &
                                          op_useLogTables,              &
                                          op_tableEnergyTolerance,      &
                                          op_maxNstepsDensityPA,        &
                                          op_maxNstepsDensityPE,        &
                                          op_maxNstepsDensityRO,        &
                                          op_maxNstepsTemperaturePA,    &
                                          op_maxNstepsTemperaturePE,    &
                                          op_maxNstepsTemperatureRO,    &
                                          op_maxTablesPA,               &
                                          op_maxTablesPE,               &
                                          op_maxTablesRO,               &
                                          op_tableKind,                 &
                                          op_tableName,                 &
                                          op_nstepsDensityPA,           &
                                          op_nstepsDensityPE,           &
                                          op_nstepsDensityRO,           &
                                          op_nstepsTemperaturePA,       &
                                          op_nstepsTemperaturePE,       &
                                          op_nstepsTemperatureRO,       &
                                          op_species2PATableIndex,      &
                                          op_species2PETableIndex,      &
                                          op_species2ROTableIndex,      &
                                          op_speciesMaxTempPATable,     &
                                          op_speciesMaxTempPETable,     &
                                          op_speciesMaxTempROTable,     &
                                          op_speciesMinTempPATable,     &
                                          op_speciesMinTempPETable,     &
                                          op_speciesMinTempROTable,     &
                                          op_tabulatedEnergyBoundaries, &
                                          op_tableDensityPA,            &
                                          op_tableDensityPE,            &
                                          op_tableDensityRO,            &
                                          op_tableTemperaturePA,        &
                                          op_tableTemperaturePE,        &
                                          op_tableTemperatureRO,        &
                                          op_PlanckAbsorptionTables,    &
                                          op_PlanckEmissionTables,      &
                                          op_RosselandTables

  use op_numericsData,             ONLY : one,ten

  use op_interface,                ONLY : op_browseTables,              &
                                          op_readTables,                &
                                          op_writeTables

  use Simulation_interface,        ONLY : Simulation_mapIntToStr

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt

  implicit none

#include "Flash.h"
#include "Opacity.h"
#include "constants.h"

  character (len=20) :: spec_str
  character (len=MAX_STRING_LENGTH) :: str
  character (len=MAX_STRING_LENGTH) :: rtpar

  logical :: needPATable
  logical :: needPETable
  logical :: needROTable
  logical :: needTable

  integer :: absorptionKind
  integer :: emissionKind
  integer :: transportKind
  integer :: indexPA
  integer :: indexPE
  integer :: indexRO
  integer :: nstepsDensityPA
  integer :: nstepsDensityPE
  integer :: nstepsDensityRO
  integer :: nstepsTemperaturePA
  integer :: nstepsTemperaturePE
  integer :: nstepsTemperatureRO
  integer :: species
  integer :: status
!
!
!    ...Safety net. Runtime parameters.
!
!
  if (op_totalSpecies < 1) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: No species present!')
  end if

  call RuntimeParameters_get ("opacity_useLogTables",   op_useLogTables)
  call RuntimeParameters_get ("op_tableEnergyTolerance", op_tableEnergyTolerance)
!
!
!    ...Allocate those arrays that depend on the # of energy groups.
!
!
  allocate (op_tabulatedEnergyBoundaries (1:op_nEnergyGroups+1), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_tabulatedEnergyBoundaries allocation failed')
  end if
!
!
!    ...Allocate those arrays that depend on the # of species.
!
!
  allocate (op_tableKind (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_tableKind allocation failed')
  end if

  allocate (op_tableName (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_tableName allocation failed')
  end if

  allocate (op_species2PATableIndex (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_species2PATableIndex allocation failed')
  end if

  allocate (op_species2PETableIndex (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_species2PETableIndex allocation failed')
  end if

  allocate (op_species2ROTableIndex (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_species2ROTableIndex allocation failed')
  end if

  allocate (op_speciesMaxTempPATable (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_speciesMaxTempPATable allocation failed')
  end if

  allocate (op_speciesMaxTempPETable (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_speciesMaxTempPETable allocation failed')
  end if

  allocate (op_speciesMaxTempROTable (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_speciesMaxTempROTable allocation failed')
  end if

  allocate (op_speciesMinTempPATable (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_speciesMinTempPATable allocation failed')
  end if

  allocate (op_speciesMinTempPETable (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_speciesMinTempPETable allocation failed')
  end if

  allocate (op_speciesMinTempROTable (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initTabulated] ERROR: op_speciesMinTempROTable allocation failed')
  end if
!
!
!    ...Initialize:   counter arrays  -> zero
!                     index arrays    -> zero
!                     max Temp arrays -> - 1.
!                     min Temp arrays -> - 1.
!
!
  op_maxTablesPA            = 0
  op_maxTablesPE            = 0
  op_maxTablesRO            = 0
  op_maxNstepsDensityPA     = 0
  op_maxNstepsDensityPE     = 0
  op_maxNstepsDensityRO     = 0
  op_maxNstepsTemperaturePA = 0
  op_maxNstepsTemperaturePE = 0
  op_maxNstepsTemperatureRO = 0

  op_species2PATableIndex   = 0
  op_species2PETableIndex   = 0
  op_species2ROTableIndex   = 0

  op_speciesMaxTempPATable  = - one
  op_speciesMaxTempPETable  = - one
  op_speciesMaxTempROTable  = - one
  op_speciesMinTempPATable  = - one
  op_speciesMinTempPETable  = - one
  op_speciesMinTempROTable  = - one
!
!
!    ...Read in the needed info from the runtime parameters and set
!       the handle arrays.
!
!
  do species = 1,op_totalSpecies
     call Simulation_mapIntToStr(species+SPECIES_BEGIN-1,spec_str,MAPBLOCK_UNK)

     write(rtpar,'(3a)') "op_", trim(spec_str), "Absorb"
     call RuntimeParameters_get(rtpar, str)
     call RuntimeParameters_mapStrToInt(str,absorptionKind)

     write(rtpar,'(3a)') "op_", trim(spec_str), "Emiss"
     call RuntimeParameters_get(rtpar, str)
     call RuntimeParameters_mapStrToInt(str,emissionKind)

     write(rtpar,'(3a)') "op_", trim(spec_str), "Trans"
     call RuntimeParameters_get(rtpar, str)
     call RuntimeParameters_mapStrToInt(str,transportKind)

     write(rtpar,'(3a)') "op_", trim(spec_str), "FileType"
     call RuntimeParameters_get(rtpar, op_tableKind(species))

     write(rtpar,'(3a)') "op_", trim(spec_str), "FileName"
     call RuntimeParameters_get(rtpar, op_tableName(species))

     if (absorptionKind == OP_TABULAR_PA) op_absorptionKind (species) = OP_TABULAR_PA
     if (absorptionKind == OP_TABULAR_PE) op_absorptionKind (species) = OP_TABULAR_PE
     if (absorptionKind == OP_TABULAR_RO) op_absorptionKind (species) = OP_TABULAR_RO

     if (  emissionKind == OP_TABULAR_PA) op_emissionKind   (species) = OP_TABULAR_PA
     if (  emissionKind == OP_TABULAR_PE) op_emissionKind   (species) = OP_TABULAR_PE
     if (  emissionKind == OP_TABULAR_RO) op_emissionKind   (species) = OP_TABULAR_RO

     if ( transportKind == OP_TABULAR_PA) op_transportKind  (species) = OP_TABULAR_PA
     if ( transportKind == OP_TABULAR_PE) op_transportKind  (species) = OP_TABULAR_PE
     if ( transportKind == OP_TABULAR_RO) op_transportKind  (species) = OP_TABULAR_RO

  end do
!
!
!    ...Determine first the overall maximal dimensions needed for allocating
!       the tables.
!
!
  do species = 1,op_totalSpecies

     needPATable  =      (op_absorptionKind (species) == OP_TABULAR_PA) &
                    .or. (op_emissionKind   (species) == OP_TABULAR_PA) &
                    .or. (op_transportKind  (species) == OP_TABULAR_PA)

     needPETable  =      (op_absorptionKind (species) == OP_TABULAR_PE) &
                    .or. (op_emissionKind   (species) == OP_TABULAR_PE) &
                    .or. (op_transportKind  (species) == OP_TABULAR_PE)

     needROTable  =      (op_absorptionKind (species) == OP_TABULAR_RO) &
                    .or. (op_emissionKind   (species) == OP_TABULAR_RO) &
                    .or. (op_transportKind  (species) == OP_TABULAR_RO)

     needTable    =       needPATable &
                    .or.  needPETable &
                    .or.  needROTable

     if (needPATable) op_maxTablesPA = op_maxTablesPA + 1
     if (needPETable) op_maxTablesPE = op_maxTablesPE + 1
     if (needROTable) op_maxTablesRO = op_maxTablesRO + 1

     if (needTable) then

         call op_browseTables (op_tableKind (species),           &
                               op_tableName (species),           &
                               needPATable,                      &
                               needPETable,                      &
                               needROTable,                      &
                                            nstepsDensityPA,     &
                                            nstepsDensityPE,     &
                                            nstepsDensityRO,     &
                                            nstepsTemperaturePA, &
                                            nstepsTemperaturePE, &
                                            nstepsTemperatureRO  )

         op_maxNstepsDensityPA     = max (nstepsDensityPA     , op_maxNstepsDensityPA    )
         op_maxNstepsDensityPE     = max (nstepsDensityPE     , op_maxNstepsDensityPE    )
         op_maxNstepsDensityRO     = max (nstepsDensityRO     , op_maxNstepsDensityRO    )
         op_maxNstepsTemperaturePA = max (nstepsTemperaturePA , op_maxNstepsTemperaturePA)
         op_maxNstepsTemperaturePE = max (nstepsTemperaturePE , op_maxNstepsTemperaturePE)
         op_maxNstepsTemperatureRO = max (nstepsTemperatureRO , op_maxNstepsTemperatureRO)

     end if

  end do
!
!
!    ...Allocate the tables and associated data for interpolation.
!       First all tables and data for the Planck Absorption (PA) opacities.
!
!
  if (op_maxTablesPA > 0) then

      if (op_maxNstepsDensityPA == 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: no PA table density grid')
      end if

      if (op_maxNstepsTemperaturePA == 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: no PA table temperature grid')
      end if

      allocate (op_nstepsDensityPA (1:op_maxTablesPA), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_nstepsDensityPA allocation failed')
      end if

      allocate (op_nstepsTemperaturePA (1:op_maxTablesPA), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_nstepsTemperaturePA allocation failed')
      end if
 
      allocate (op_tableDensityPA (1:op_maxNstepsDensityPA,1:op_maxTablesPA), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_tableDensityPA allocation failed')
      end if

      allocate (op_tableTemperaturePA (1:op_maxNstepsTemperaturePA,1:op_maxTablesPA), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_tableTemperaturePA allocation failed')
      end if

      allocate (op_PlanckAbsorptionTables (1:op_maxNstepsTemperaturePA, &
                                           1:op_maxNstepsDensityPA,     &
                                           1:op_nEnergyGroups,          &
                                           1:op_maxTablesPA             ), stat = status)
      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_PlanckAbsorptionTables allocation failed')
      end if

  end if
!
!
!    ...Next all tables and data for the Planck Emission (PE) opacities.
!
!
  if (op_maxTablesPE > 0) then

      if (op_maxNstepsDensityPE == 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: no PE table density grid')
      end if

      if (op_maxNstepsTemperaturePE == 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: no PE table temperature grid')
      end if

      allocate (op_nstepsDensityPE (1:op_maxTablesPE), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_nstepsDensityPE allocation failed')
      end if

      allocate (op_nstepsTemperaturePE (1:op_maxTablesPE), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_nstepsTemperaturePE allocation failed')
      end if

      allocate (op_tableDensityPE (1:op_maxNstepsDensityPE,1:op_maxTablesPE), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_tableDensityPE allocation failed')
      end if

      allocate (op_tableTemperaturePE (1:op_maxNstepsTemperaturePE,1:op_maxTablesPE), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_tableTemperaturePE allocation failed')
      end if

      allocate (op_PlanckEmissionTables (1:op_maxNstepsTemperaturePE, &
                                         1:op_maxNstepsDensityPE,     &
                                         1:op_nEnergyGroups,          &
                                         1:op_maxTablesPE             ), stat = status)
      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_PlanckEmissionTables allocation failed')
      end if
  end if
!
!
!    ...Next all tables and data for the Rosseland (RO) opacities.
!
!
  if (op_maxTablesRO > 0) then

      if (op_maxNstepsDensityRO == 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: no RO table density grid')
      end if

      if (op_maxNstepsTemperatureRO == 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: no RO table temperature grid')
      end if

      allocate (op_nstepsDensityRO (1:op_maxTablesRO), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_nstepsDensityRO allocation failed')
      end if

      allocate (op_nstepsTemperatureRO (1:op_maxTablesRO), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_nstepsTemperatureRO allocation failed')
      end if

      allocate (op_tableDensityRO (1:op_maxNstepsDensityRO,1:op_maxTablesRO), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_tableDensityRO allocation failed')
      end if

      allocate (op_tableTemperatureRO (1:op_maxNstepsTemperatureRO,1:op_maxTablesRO), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_tableTemperatureRO allocation failed')
      end if

      allocate (op_RosselandTables (1:op_maxNstepsTemperatureRO, &
                                    1:op_maxNstepsDensityRO,     &
                                    1:op_nEnergyGroups,          &
                                    1:op_maxTablesRO             ), stat = status)
      if (status > 0) then
          call Driver_abortFlash ('[op_initTabulated] ERROR: op_RosselandTables allocation failed')
      end if
  end if
!
!
!    ...Read the tables and store the needed opacity values and all associated data
!       into the arrays.
!
!
  indexPA = 0
  indexPE = 0
  indexRO = 0

  do species = 1,op_totalSpecies

     needPATable  =      (op_absorptionKind (species) == OP_TABULAR_PA) &
                    .or. (op_emissionKind   (species) == OP_TABULAR_PA) &
                    .or. (op_transportKind  (species) == OP_TABULAR_PA)

     needPETable  =      (op_absorptionKind (species) == OP_TABULAR_PE) &
                    .or. (op_emissionKind   (species) == OP_TABULAR_PE) &
                    .or. (op_transportKind  (species) == OP_TABULAR_PE)

     needROTable  =      (op_absorptionKind (species) == OP_TABULAR_RO) &
                    .or. (op_emissionKind   (species) == OP_TABULAR_RO) &
                    .or. (op_transportKind  (species) == OP_TABULAR_RO)

     needTable    =       needPATable &
                    .or.  needPETable &
                    .or.  needROTable

     if (needPATable) then
         indexPA = indexPA + 1
         op_species2PATableIndex (species) = indexPA
     end if

     if (needPETable) then
         indexPE = indexPE + 1
         op_species2PETableIndex (species) = indexPE
     end if

     if (needROTable) then
         indexRO = indexRO + 1
         op_species2ROTableIndex (species) = indexRO
     end if

     if (needTable) then

         call op_readTables (op_tableKind (species), &
                             op_tableName (species), &
                             needPATable,            &
                             needPETable,            &
                             needROTable,            &
                             indexPA,                &
                             indexPE,                &
                             indexRO                 )
     end if

  end do
!
!
!    ...Determine the temperature bounds (in K) for each species/Table combination.
!       If no Table for a particular species was created, the temperature bound
!       will retain its original -1 value.
!
!
  do species = 1,op_totalSpecies

     indexPA = op_species2PATableIndex (species)
     indexPE = op_species2PETableIndex (species)
     indexRO = op_species2ROTableIndex (species)

     if (indexPA > 0) then

         nstepsTemperaturePA                = op_nstepsTemperaturePA (indexPA)

         if (op_useLogTables) then
             op_speciesMinTempPATable (species) = ten ** op_tableTemperaturePA  (                  1,indexPA)
             op_speciesMaxTempPATable (species) = ten ** op_tableTemperaturePA  (nstepsTemperaturePA,indexPA)
         else
             op_speciesMinTempPATable (species) = op_tableTemperaturePA  (                  1,indexPA)
             op_speciesMaxTempPATable (species) = op_tableTemperaturePA  (nstepsTemperaturePA,indexPA)
         end if
     end if

     if (indexPE > 0) then

         nstepsTemperaturePE                = op_nstepsTemperaturePE (indexPE)

         if (op_useLogTables) then
             op_speciesMinTempPETable (species) = ten ** op_tableTemperaturePE  (                  1,indexPE)
             op_speciesMaxTempPETable (species) = ten ** op_tableTemperaturePE  (nstepsTemperaturePE,indexPE)
         else
             op_speciesMinTempPETable (species) = op_tableTemperaturePE  (                  1,indexPE)
             op_speciesMaxTempPETable (species) = op_tableTemperaturePE  (nstepsTemperaturePE,indexPE)
         end if
     end if

     if (indexRO > 0) then

         nstepsTemperatureRO                = op_nstepsTemperatureRO (indexRO)

         if (op_useLogTables) then
             op_speciesMinTempROTable (species) = ten ** op_tableTemperatureRO  (                  1,indexRO)
             op_speciesMaxTempROTable (species) = ten ** op_tableTemperatureRO  (nstepsTemperatureRO,indexRO)
         else
             op_speciesMinTempROTable (species) = op_tableTemperatureRO  (                  1,indexRO)
             op_speciesMaxTempROTable (species) = op_tableTemperatureRO  (nstepsTemperatureRO,indexRO)
         end if
     end if

  end do
!
!
!    ...Write out the Opacity tables and associated data (if requested).
!
!
  if (op_writeOpacityInfo) then
      call op_writeTables ()
  end if
!
!
!    ...Set initialization status.
!
!
  op_initializedTabulated = .true.
!
!
!   ...Ready! 
!
!
  return
end subroutine op_initTabulated
