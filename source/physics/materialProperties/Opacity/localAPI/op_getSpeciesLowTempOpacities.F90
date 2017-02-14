!!****if* source/physics/materialProperties/Opacity/localAPI/op_getSpeciesLowTempOpacities
!!
!! NAME
!!
!!  op_getSpeciesLowTempOpacities
!!
!! SYNOPSIS
!!
!!  call op_getSpeciesLowTempOpacities (integer (in)  :: species,
!!                                      real    (in)  :: speciesTemperature,
!!                                      integer (in)  :: speciesEnergyGroup)
!!
!! DESCRIPTION
!!
!!  Computes absorption, emission and transport opacities from tabulated low temperature
!!  opacity values for a particular (temperature, energy group, species) triple.
!!  Extracts via interpolation low temperature Planck and Rosseland opacities for a
!!  specific temperature, energy group and species and copies them to the absorption,
!!  emission and transport opacities. The interpolation method is linear.
!!
!!  The mapping between the low temperature Planck and Rosseland opacities and the
!!  absorption, emission and transport opacities is as follows:
!!
!!              absorption opacity = low temperature Planck
!!              emission   opacity = low temperature Planck
!!              transport  opacity = low temperature Rosseland
!!
!! ARGUMENTS
!!
!!  species            : The species index
!!  speciesTemperature : The species temperature (in K)
!!  speciesEnergyGroup : The species energy group index
!!
!!***
subroutine op_getSpeciesLowTempOpacities (species,            &
                                          speciesTemperature, &
                                          speciesEnergyGroup  )

  implicit none

  integer, intent (in)  :: species
  integer, intent (in)  :: speciesEnergyGroup
  real,    intent (in)  :: speciesTemperature

  return
end subroutine op_getSpeciesLowTempOpacities
