!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/Opacity_data
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
!!  Defines and stores the local data for the main Opacity unit. Data specific
!!  for each opacity method (constant, tabulated, lowtemp, etc) is stored in the
!!  corresponding subdirectories /method/<subdirectory>.
!!  
!! PARAMETERS
!!
!!  op_maxAtomicNumber           : the maximum atomic number that can be treated
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
  logical, save :: op_ignoreLowTemp
  logical, save :: op_writeOpacityInfo
  integer, save :: op_nEnergyGroups
  integer, save :: op_totalElements
  integer, save :: op_totalSpecies
  integer, save :: op_maxNelementsPerSpecies
  integer, save :: op_maxNspeciesPerElement
  integer, save :: op_globalMe
  real,    save :: op_cellMassDensity
  real,    save :: op_cellTemperature
  real,    save :: op_totalIonNumberDensity

  integer, parameter :: op_maxAtomicNumber           = 100         ! up to Fermium Fm
  real,    parameter :: op_absoluteLowestTemperature = 1.0
  real,    parameter :: op_energyDifferenceTolerance = 1.E-8

  character (len=12), allocatable, save :: op_atomName      (:)
  real,               allocatable, save :: op_atomWeight    (:)

  integer, allocatable, save :: op_absorptionKind           (:)
  integer, allocatable, save :: op_emissionKind             (:)
  integer, allocatable, save :: op_transportKind            (:)
  integer, allocatable, save :: op_elementNumberofSpecies   (:)
  integer, allocatable, save :: op_speciesNumberofElements  (:)
  integer, allocatable, save :: op_element2AtomicNumber     (:)
  real,    allocatable, save :: op_speciesWeights           (:)
  real,    allocatable, save :: op_speciesLowTempCutoff     (:)
  real,    allocatable, save :: op_cellSpeciesMassFractions (:)
  real,    allocatable, save :: op_ionNumberDensity         (:)
  real,    allocatable, save :: op_energyGroupBoundaries    (:)

  integer, allocatable, save :: op_elements2Species         (:,:)
  integer, allocatable, save :: op_species2Elements         (:,:)
  real,    allocatable, save :: op_species2FractionElements (:,:)
  real,    allocatable, save :: op_speciesOpacities         (:,:)

end module Opacity_data
