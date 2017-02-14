!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_finalizeTabulated
!!
!! NAME
!!
!!  op_finalizeTabulated
!!
!! SYNOPSIS
!!
!!  call op_finalizeTabulated ()
!!
!! DESCRIPTION
!!
!!  Finalizes the tabulated section of the opacity unit.
!!
!! ARGUMENTS
!!
!!***
subroutine op_finalizeTabulated ()

  use op_tabulatedData, ONLY : op_tableKind,                 &
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

  implicit none
!
!
!   ...Deallocate the arrays.
!
!
  if (allocated (op_tableKind)                ) deallocate (op_tableKind)
  if (allocated (op_tableName)                ) deallocate (op_tableName)
  if (allocated (op_nstepsDensityPA)          ) deallocate (op_nstepsDensityPA)
  if (allocated (op_nstepsDensityPE)          ) deallocate (op_nstepsDensityPE)
  if (allocated (op_nstepsDensityRO)          ) deallocate (op_nstepsDensityRO)
  if (allocated (op_nstepsTemperaturePA)      ) deallocate (op_nstepsTemperaturePA)
  if (allocated (op_nstepsTemperaturePE)      ) deallocate (op_nstepsTemperaturePE)
  if (allocated (op_nstepsTemperatureRO)      ) deallocate (op_nstepsTemperatureRO)
  if (allocated (op_species2PATableIndex)     ) deallocate (op_species2PATableIndex)
  if (allocated (op_species2PETableIndex)     ) deallocate (op_species2PETableIndex)
  if (allocated (op_species2ROTableIndex)     ) deallocate (op_species2ROTableIndex)
  if (allocated (op_speciesMaxTempPATable)    ) deallocate (op_speciesMaxTempPATable)
  if (allocated (op_speciesMaxTempPETable)    ) deallocate (op_speciesMaxTempPETable)
  if (allocated (op_speciesMaxTempROTable)    ) deallocate (op_speciesMaxTempROTable)
  if (allocated (op_speciesMinTempPATable)    ) deallocate (op_speciesMinTempPATable)
  if (allocated (op_speciesMinTempPETable)    ) deallocate (op_speciesMinTempPETable)
  if (allocated (op_speciesMinTempROTable)    ) deallocate (op_speciesMinTempROTable)
  if (allocated (op_tabulatedEnergyBoundaries)) deallocate (op_tabulatedEnergyBoundaries)
  if (allocated (op_tableDensityPA)           ) deallocate (op_tableDensityPA)
  if (allocated (op_tableDensityPE)           ) deallocate (op_tableDensityPE)
  if (allocated (op_tableDensityRO)           ) deallocate (op_tableDensityRO)
  if (allocated (op_tableTemperaturePA)       ) deallocate (op_tableTemperaturePA)
  if (allocated (op_tableTemperaturePE)       ) deallocate (op_tableTemperaturePE)
  if (allocated (op_tableTemperatureRO)       ) deallocate (op_tableTemperatureRO)
  if (allocated (op_PlanckAbsorptionTables)   ) deallocate (op_PlanckAbsorptionTables) 
  if (allocated (op_PlanckEmissionTables)     ) deallocate (op_PlanckEmissionTables)
  if (allocated (op_RosselandTables)          ) deallocate (op_RosselandTables)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_finalizeTabulated
