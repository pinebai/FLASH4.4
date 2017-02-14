!!****if* source/physics/Diffuse/DiffuseMain/Diffuse_data
!!
!! NAME
!!
!!  Diffuse_data
!!
!! SYNOPSIS
!!  use Diffuse_data
!!
!! DESCRIPTION
!!
!!  Defines and stores the local data for Unit Diffuse. All the variables
!!  defined here are initialized in Diffuse_init() by calling 
!!  RuntimeParameters_get subroutine. These data variables are for 
!!  Unit Scope ->  Diffuse, and used by all three major diffuse subroutines,
!!  Diffuse_species(for mass diffuzivity), Diffuse_therm (for 
!!  thermal Conductivity), and Diffuse_visc (for viscosity).
!!  
!!***


Module Diffuse_data

  implicit none

#include "constants.h"

  integer, save :: diff_meshMe, diff_meshNumProcs, diff_meshComm
  logical, save :: useDiffuseTherm,useDiffuseVisc,useDiffuseMagneticResistivity,&
                   useDiffuseSpecies,useDiffuse
  logical, save :: useDiffuseComputeDtTherm,useDiffuseComputeDtVisc,useDiffuseComputeDtMagnetic,&
                   useDiffuseComputeDtSpecies
  logical, save :: diff_useEleCond
  logical, save :: diff_useIonCond
  integer, save :: thermal_diff_method
  real,    save :: dt_diff_factor, diffusion_cutoff_density
  logical, save :: diff_geometricMeanDiff
  integer, save :: diff_geometry
  integer, save, dimension(MDIM) :: diff_dirGeom

  type diff_dbgContext_t
     integer :: component
     integer :: group
  end type diff_dbgContext_t

  ! Structure that holds context information on the current operation,
  ! for debugging
  type(diff_dbgContext_t), save :: diff_dbgContext

  real, save :: diff_scaleFactThermFlux
  real, save :: diff_singleSpeciesZ
  real, save :: diff_singleSpeciesA

  integer, save :: diff_eleFlMode ! Flux limiter mode for electrons
  integer, save :: diff_ionFlMode ! Flux limiter mode for ions
  real, save :: diff_eleFlCoef ! Electron flux limiter coefficient
  real, save :: diff_ionFlCoef ! Ion flux limiter coefficient

  ! Physical constants:
  real, save :: diff_asol ! Radiation constant
  real, save :: diff_mele ! Electron mass
  real, save :: diff_speedlt ! Speed of light
  real, save :: diff_boltz ! Boltzmann constant
  real, save :: diff_avo ! Avogadros number

end Module Diffuse_data
