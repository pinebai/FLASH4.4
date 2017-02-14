!!****if* source/physics/materialProperties/Opacity/localAPI/op_getSpeciesPATableOpacity
!!
!! NAME
!!
!!  op_getSpeciesPATableOpacity
!!
!! SYNOPSIS
!!
!!  call op_getSpeciesPATableOpacity (integer (in)  :: species,
!!                                    real    (in)  :: speciesTemperature,
!!                                    real    (in)  :: speciesDensity,
!!                                    integer (in)  :: speciesEnergyGroup,
!!                                    real    (out) :: opacityPA)
!!
!! DESCRIPTION
!!
!!  Extracts via interpolation a Planck Absorption opacity value for
!!  the specified energy group and species. The interpolation method
!!  is bilinear.
!!
!! ARGUMENTS
!!
!!   species            : The species index
!!   speciesTemperature : The species temperature
!!   speciesDensity     : The species density
!!   speciesEnergyGroup : The species energy group index
!!   opacityPA          : The value of the determined Planck Absorption opacity (in cm^2/g)
!!
!!***
subroutine op_getSpeciesPATableOpacity (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        speciesEnergyGroup, &
                                        opacityPA           )

  implicit none
  
  integer, intent (in)  :: species
  integer, intent (in)  :: speciesEnergyGroup
  real,    intent (in)  :: speciesTemperature
  real,    intent (in)  :: speciesDensity
  real,    intent (out) :: opacityPA

  opacityPA = 0.0

  return
end subroutine op_getSpeciesPATableOpacity
