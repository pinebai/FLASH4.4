!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/Constant-diff/Conductivity_data
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
!!  Stores the local data for the constant conductivity implementation.
!!
!!
!!***
Module Conductivity_data

  implicit none
  
  logical, save::  cond_useConductivity

  real, save :: cond_diffConstant

end Module Conductivity_data
