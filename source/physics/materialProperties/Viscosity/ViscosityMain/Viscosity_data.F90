!!****if* source/physics/materialProperties/Viscosity/ViscosityMain/Viscosity_data
!!
!! NAME
!!
!!  Viscosity_data
!!
!! SYNOPSIS
!!
!!  use Viscosity_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for the viscosity implementation.
!!
!!
!!***
Module Viscosity_data

  implicit none

  logical, save::  visc_useViscosity

  integer, save :: visc_whichCoefficientIsConst
  real, save :: visc_diffNu, visc_diffMu
! The following three are currently not used in the Constant
! Viscosity implementation.
  real, save :: viscTempLow, viscTempHigh
  real, save :: viscSuppressFactor 

end Module Viscosity_data
