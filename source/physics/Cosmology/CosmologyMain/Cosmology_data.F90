!!****if* source/physics/Cosmology/CosmologyMain/Cosmology_data
!!
!! NAME
!!  Cosmology_data
!!
!! SYNOPSIS
!!
!!  use Cosmology_data
!!
!! DESCRIPTION 
!!  
!!  Holds all the Cosmology data that is needed by the Cosmology unit     
!!
!! ARGUMENTS
!!
!!  none    
!!
!!
!!***
Module Cosmology_data
  integer, save :: csm_meshMe, csm_meshNumProcs
  real, save    :: csm_smlrho, csm_smalle, csm_smallp, csm_eintSwitch
  logical, save :: csm_eintSwitchExist, csm_computeRedshiftOnly
  real, save    :: csm_hubble,csm_omega,csm_lambda,csm_curv,csm_newton
  real, save    :: csm_c,csm_maxScaleChange, csm_redshiftFinal,csm_scaleFactor
  real, save    :: csm_oldScaleFactor,csm_baryon
end Module Cosmology_data

