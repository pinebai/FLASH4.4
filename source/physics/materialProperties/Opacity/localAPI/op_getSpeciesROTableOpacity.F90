!!****if* source/physics/materialProperties/Opacity/localAPI/op_getSpeciesROTableOpacity
!!
!! NAME
!!
!!  op_getSpeciesROTableOpacity
!!
!! SYNOPSIS
!!
!!  call op_getSpeciesROTableOpacity (integer (in)  :: species,
!!                                    real    (in)  :: speciesTemperature,
!!                                    real    (in)  :: speciesDensity,
!!                                    integer (in)  :: speciesEnergyGroup,
!!                                    real    (out) :: opacityRO)
!!
!! DESCRIPTION
!!
!!  Extracts via interpolation a Rosseland opacity value for
!!  the specified energy group and species. The interpolation method
!!  is bilinear.
!!
!! ARGUMENTS
!!
!!   species            : The species index
!!   speciesTemperature : The species temperature
!!   speciesDensity     : The species density
!!   speciesEnergyGroup : The species energy group index
!!   opacityRO          : The value of the determined Rosseland opacity (in cm^2/g)
!!
!!***
subroutine op_getSpeciesROTableOpacity (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        speciesEnergyGroup, &
                                        opacityRO           )

  implicit none

  integer, intent (in)  :: species
  integer, intent (in)  :: speciesEnergyGroup
  real,    intent (in)  :: speciesTemperature
  real,    intent (in)  :: speciesDensity
  real,    intent (out) :: opacityRO

  opacityRO = 0.0

  return
end subroutine op_getSpeciesROTableOpacity
