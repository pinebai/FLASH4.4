!!****if* source/physics/materialProperties/MagneticResistivity/MagneticResistivityMain/Constant/MagneticResistivity_data
!!
!! NAME
!!  MagneticResistivity_data
!!
!! SYNOPSIS
!!  MagneticResistivity_data()
!!
!! DESCRIPTION
!!  A placeholder for the magnetic resistivity data
!!
!! ARGUMENTS
!!  No arguments
!!
!!***

module MagneticResistivity_data

  implicit none

  logical, save      :: mag_useMagneticResistivity
  
  real,         save :: mResistivity
  real,    PARAMETER :: c = 29979245800. ! CGS units
  character(4), save :: mUnit

end module MagneticResistivity_data
