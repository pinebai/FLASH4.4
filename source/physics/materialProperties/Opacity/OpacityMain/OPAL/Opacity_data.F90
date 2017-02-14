!!****if* source/physics/materialProperties/Opacity/OpacityMain/OPAL/Opacity_data
!!
!! NAME
!!
!!  Opacity_data
!!
!! SYNOPSIS
!!  use Opacity_data
!!
!! DESCRIPTION
!!
!!  Defines and stores local data for the OPAL implementation of the
!!  Opacity unit.
!!  
!!  
!! PARAMETERS
!!
!!  op_absoluteLowestTemperature : the lowest temperature (in Kelvin) that is still
!!                                 recognized by the opacity unit
!!  op_energyDifferenceTolerance : the epsilon value for the energies. If two energies
!!                                 E and E' obey |E - E'| < epsilon, then E and E' are
!!                                 considered equal.
!!  
!!***
module Opacity_data
  
  implicit none

  logical, save :: op_useOpacity
  logical, save :: op_useMGD
  logical, save :: op_useLogTables
  logical, save :: op_writeOpacityInfo
  integer, save :: op_nEnergyGroups
  integer, save :: op_totalElements
  integer, save :: op_opalNumHydrogenAbundances
  integer, save :: op_maxNelementsPerSpecies
  integer, save :: op_maxNspeciesPerElement
  integer, save :: op_globalMe
  real,    save :: op_cellMassDensity
  real,    save :: op_cellTemperature
  real,    save :: op_totalIonNumberDensity

  real,    parameter :: op_absoluteLowestTemperature = 1.0
  real,    parameter :: op_energyDifferenceTolerance = 1.E-8


  integer, allocatable, save :: op_absorptionKind           (:)
  integer, allocatable, save :: op_emissionKind             (:)
  integer, allocatable, save :: op_transportKind            (:)
  integer, allocatable, save :: op_elementNumberofSpecies   (:)
  integer, allocatable, save :: op_speciesNumberofElements  (:)
  real,    allocatable, save :: op_speciesWeights           (:)
  real,    allocatable, save :: op_speciesLowTempCutoff     (:)
  real,    allocatable, save :: op_ionNumberDensity         (:)
  real,    allocatable, save :: op_energyGroupBoundaries    (:)

  real,    allocatable, save :: op_speciesOpacities         (:,:)

  integer, save :: op_hydrogenMassFracVar
  real,    save :: op_hydrogenMassFrac

  real, save :: op_emitConst
  real, save :: op_absorbConst

end module Opacity_data
