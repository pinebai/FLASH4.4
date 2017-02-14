!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/PowerLaw-gray/Conductivity_data
!!
!! NAME
!!
!!  Conductivity_data
!!
!! SYNOPSIS
!!
!!  use Conductivity_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for the Powerlaw gray Conductivity implementation.
!!
!!
!!***
Module Conductivity_data

  implicit none

  logical, save::  cond_useConductivity !see diff_* in Diffuse_data instead:  , useRaddiffusion, useElecond
  
  integer, save::  cond_meshMe

  real, save :: cond_TemperatureExponent, cond_K0, cond_alpha, Raddiff_K0r, Raddiff_TemperatureExponent

end Module Conductivity_data
