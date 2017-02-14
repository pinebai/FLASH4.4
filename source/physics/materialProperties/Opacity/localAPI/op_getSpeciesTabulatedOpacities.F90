!!****if* source/physics/materialProperties/Opacity/localAPI/op_getSpeciesTabulatedOpacities
!!
!! NAME
!!
!!  op_getSpeciesTabulatedOpacities
!!
!! SYNOPSIS
!!
!!  call op_getSpeciesTabulatedOpacities (integer (in) :: species,
!!                                        real    (in) :: speciesTemperature,
!!                                        real    (in) :: speciesDensity,
!!                                        integer (in) :: speciesEnergyGroup,
!!                                        logical (in) :: needPATable,
!!                                        logical (in) :: needPETable,
!!                                        logical (in) :: needROTable)
!!
!! DESCRIPTION
!!
!!  Computes absorption, emission and transport opacities from tabulated
!!  values for a particular cell and species within that cell.
!!
!! ARGUMENTS
!!
!!   species            : The species index
!!   speciesTemperature : The species temperature
!!   speciesDensity     : The species density
!!   speciesEnergyGroup : The species energy group index
!!   needPATable        : Knob to activate opacity extraction from the Planck Absorption tables
!!   needPETable        : Knob to activate opacity extraction from the Planck Emission tables
!!   needROTable        : Knob to activate opacity extraction from the Rosseland tables
!!
!!***
subroutine op_getSpeciesTabulatedOpacities (species,            &
                                            speciesTemperature, &
                                            speciesDensity,     &
                                            speciesEnergyGroup, &
                                            needPATable,        &
                                            needPETable,        &
                                            needROTable         )

  implicit none
  
  logical, intent (in) :: needPATable
  logical, intent (in) :: needPETable
  logical, intent (in) :: needROTable
  integer, intent (in) :: species
  integer, intent (in) :: speciesEnergyGroup
  real,    intent (in) :: speciesTemperature
  real,    intent (in) :: speciesDensity

  return
end subroutine op_getSpeciesTabulatedOpacities
