!!****if* source/physics/materialProperties/Opacity/localAPI/op_getSpeciesPETableOpacity
!!
!! NAME
!!
!!  op_getSpeciesPETableOpacity
!!
!! SYNOPSIS
!!
!!  call op_getSpeciesPETableOpacity (integer (in)  :: species,
!!                                    real    (in)  :: speciesTemperature,
!!                                    real    (in)  :: speciesDensity,
!!                                    integer (in)  :: speciesEnergyGroup,
!!                                    real    (out) :: opacityPE)
!!
!! DESCRIPTION
!!
!!  Extracts via interpolation a Planck Emission opacity value for
!!  the specified energy group and species. The interpolation method
!!  is bilinear.
!!
!! ARGUMENTS
!!
!!   species            : The species index
!!   speciesTemperature : The species temperature
!!   speciesDensity     : The species density
!!   speciesEnergyGroup : The species energy group index
!!   opacityPE          : The value of the determined Planck Emission opacity (in cm^2/g)
!!
!!***
subroutine op_getSpeciesPETableOpacity (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        speciesEnergyGroup, &
                                        opacityPE           )

  implicit none

  integer, intent (in)  :: species
  integer, intent (in)  :: speciesEnergyGroup
  real,    intent (in)  :: speciesTemperature
  real,    intent (in)  :: speciesDensity
  real,    intent (out) :: opacityPE

  opacityPE = 0.0

  return
end subroutine op_getSpeciesPETableOpacity
