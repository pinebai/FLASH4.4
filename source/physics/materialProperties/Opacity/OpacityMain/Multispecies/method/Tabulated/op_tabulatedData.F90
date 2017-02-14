!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_tabulatedData
!!
!! NAME
!!
!!  op_tabulatedData
!!
!! SYNOPSIS
!!
!!  use op_tabulatedData
!!
!! DESCRIPTION
!!
!!  Defines and stores the local data for the tabulated opacities.
!!  
!!***
module op_tabulatedData
  
  implicit none

  logical, save :: op_initializedTabulated = .false.
  logical, save :: op_useLogTables
  real, save :: op_tableEnergyTolerance

  integer, save :: op_maxNstepsDensityPA
  integer, save :: op_maxNstepsDensityPE
  integer, save :: op_maxNstepsDensityRO
  integer, save :: op_maxNstepsTemperaturePA
  integer, save :: op_maxNstepsTemperaturePE
  integer, save :: op_maxNstepsTemperatureRO
  integer, save :: op_maxTablesPA
  integer, save :: op_maxTablesPE
  integer, save :: op_maxTablesRO

  character (len=80), allocatable, save :: op_tableKind      (:)
  character (len=80), allocatable, save :: op_tableName      (:)

  integer, allocatable, save :: op_nstepsDensityPA           (:)
  integer, allocatable, save :: op_nstepsDensityPE           (:)
  integer, allocatable, save :: op_nstepsDensityRO           (:)
  integer, allocatable, save :: op_nstepsTemperaturePA       (:)
  integer, allocatable, save :: op_nstepsTemperaturePE       (:)
  integer, allocatable, save :: op_nstepsTemperatureRO       (:)
  integer, allocatable, save :: op_species2PATableIndex      (:)
  integer, allocatable, save :: op_species2PETableIndex      (:)
  integer, allocatable, save :: op_species2ROTableIndex      (:)
  real,    allocatable, save :: op_speciesMaxTempPATable     (:)
  real,    allocatable, save :: op_speciesMaxTempPETable     (:)
  real,    allocatable, save :: op_speciesMaxTempROTable     (:)
  real,    allocatable, save :: op_speciesMinTempPATable     (:)
  real,    allocatable, save :: op_speciesMinTempPETable     (:)
  real,    allocatable, save :: op_speciesMinTempROTable     (:)
  real,    allocatable, save :: op_tabulatedEnergyBoundaries (:)

  real,    allocatable, save :: op_tableDensityPA            (:,:)
  real,    allocatable, save :: op_tableDensityPE            (:,:)
  real,    allocatable, save :: op_tableDensityRO            (:,:)
  real,    allocatable, save :: op_tableTemperaturePA        (:,:)
  real,    allocatable, save :: op_tableTemperaturePE        (:,:)
  real,    allocatable, save :: op_tableTemperatureRO        (:,:)

  real,    allocatable, save :: op_PlanckAbsorptionTables    (:,:,:,:)
  real,    allocatable, save :: op_PlanckEmissionTables      (:,:,:,:)
  real,    allocatable, save :: op_RosselandTables           (:,:,:,:)

end module op_tabulatedData

