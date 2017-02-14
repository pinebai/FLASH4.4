!!****if* source/physics/utilities/PlasmaState/PlasmaStateMain/PlasmaState_data
!!
!! NAME
!!
!!  PlasmaState_data
!!
!! SYNOPSIS
!!  use PlasmaState_data
!!
!! DESCRIPTION
!!
!!  Defines and stores the local data of the PlasmaState unit.
!!  
!! PARAMETERS
!!
!!  usePlasmaState               : flags whether the PlasmaState unit is being used at all
!!  
!!***
module PlasmaState_data
  
  implicit none

  logical, save :: pls_usePlasmaState
  logical, save :: pls_useMGD
  logical, save :: pls_useLogTables
  logical, save :: pls_ignoreLowTemp
  logical, save :: pls_writePlasmaStateInfo
  integer, save :: pls_nEnergyGroups
  integer, save :: pls_totalElements
  integer, save :: pls_totalSpecies
  integer, save :: pls_maxNelementsPerSpecies
  integer, save :: pls_maxNspeciesPerElement
  integer, save :: pls_meshMe
  integer, save :: pls_globalMe
  real,    save :: pls_cellMassDensity
  real,    save :: pls_cellTemperature
  real,    save :: pls_totalIonNumberDensity

  integer, parameter :: pls_maxAtomicNumber           = 100         ! up to Fermium Fm
  real,    parameter :: pls_absoluteLowestTemperature = 1.0
  real,    parameter :: pls_energyDifferenceTolerance = 1.E-8

!!$  character (len=12), allocatable, save :: pls_atomName      (:)
!!$  real,               allocatable, save :: pls_atomWeight    (:)

  integer, allocatable, save :: pls_elementNumberofSpecies   (:)
  integer, allocatable, save :: pls_speciesNumberofElements  (:)
  integer, allocatable, save :: pls_element2AtomicNumber     (:)
  integer, allocatable, save :: pls_element2AtomicWeight     (:) !KW
  real,    allocatable, save :: pls_speciesWeights           (:)
  real,    allocatable, save :: pls_cellSpeciesMassFractions (:)
  real,    allocatable, save :: pls_ionNumberDensity         (:)

  integer, allocatable, save :: pls_elements2Species         (:,:)
  integer, allocatable, save :: pls_species2Elements         (:,:)
  real,    allocatable, save :: pls_species2FractionElements (:,:)

end module PlasmaState_data
