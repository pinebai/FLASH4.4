!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_getSpeciesTabulatedOpacities
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
!!  values for a particular (temperature, density, energy group, species)
!!  quadruple.
!!
!! ARGUMENTS
!!
!!  species            : The species index
!!  speciesTemperature : The species temperature
!!  speciesDensity     : The species density
!!  speciesEnergyGroup : The species energy group index
!!  needPATable        : Knob to activate opacity extraction from the Planck Absorption tables
!!  needPETable        : Knob to activate opacity extraction from the Planck Emission tables
!!  needROTable        : Knob to activate opacity extraction from the Rosseland tables
!!
!!***
subroutine op_getSpeciesTabulatedOpacities (species,            &
                                            speciesTemperature, &
                                            speciesDensity,     &
                                            speciesEnergyGroup, &
                                            needPATable,        &
                                            needPETable,        &
                                            needROTable         )

  use Opacity_data, ONLY : op_absorptionKind,     &
                           op_emissionKind,       &
                           op_transportKind,      &
                           op_speciesOpacities

  use op_interface, ONLY : op_getSpeciesPATableOpacity, &
                           op_getSpeciesPETableOpacity, &
                           op_getSpeciesROTableOpacity

  implicit none

#include "Opacity.h"
#include "constants.h"

  logical, intent (in) :: needPATable
  logical, intent (in) :: needPETable
  logical, intent (in) :: needROTable
  integer, intent (in) :: species
  integer, intent (in) :: speciesEnergyGroup
  real,    intent (in) :: speciesTemperature
  real,    intent (in) :: speciesDensity

  real :: opacityPA
  real :: opacityPE
  real :: opacityRO
!
!
!   ...Extract only the necessary opacities from the tables.
!
!
  if (needPATable) then

      call op_getSpeciesPATableOpacity (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        speciesEnergyGroup, &
                                        opacityPA           )
  end if

  if (needPETable) then

      call op_getSpeciesPETableOpacity (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        speciesEnergyGroup, &
                                        opacityPE           )
  end if

  if (needROTable) then

      call op_getSpeciesROTableOpacity (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        speciesEnergyGroup, &
                                        opacityRO           )
  end if
!
!
!   ...Store the absorption opacity for the current species.
!
!
  if (     op_absorptionKind (species) == OP_TABULAR_PA) then
      op_speciesOpacities (ABSORPTION,species) = opacityPA
  else if (op_absorptionKind (species) == OP_TABULAR_PE) then
      op_speciesOpacities (ABSORPTION,species) = opacityPE
  else if (op_absorptionKind (species) == OP_TABULAR_RO) then
      op_speciesOpacities (ABSORPTION,species) = opacityRO
  end if
!
!
!   ...Store the emission opacity for the current species.
!
!
  if (     op_emissionKind (species) == OP_TABULAR_PA) then
      op_speciesOpacities (  EMISSION,species) = opacityPA
  else if (op_emissionKind (species) == OP_TABULAR_PE) then
      op_speciesOpacities (  EMISSION,species) = opacityPE
  else if (op_emissionKind (species) == OP_TABULAR_RO) then
      op_speciesOpacities (  EMISSION,species) = opacityRO
  end if
!
!
!   ...Store the transport opacity for the current species.
!
!
  if (     op_transportKind (species) == OP_TABULAR_PA) then
      op_speciesOpacities ( TRANSPORT,species) = opacityPA
  else if (op_transportKind (species) == OP_TABULAR_PE) then
      op_speciesOpacities ( TRANSPORT,species) = opacityPE
  else if (op_transportKind (species) == OP_TABULAR_RO) then
      op_speciesOpacities ( TRANSPORT,species) = opacityRO
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_getSpeciesTabulatedOpacities
