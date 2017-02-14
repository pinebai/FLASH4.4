!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/Opacity
!!
!! NAME
!!
!!  Opacity
!!
!! SYNOPSIS
!!
!!  call Opacity (real    (in)  :: soln(:),
!!                integer (in)  :: ngrp,
!!                real    (out) :: opacityAbsorption,
!!                real    (out) :: opacityEmission,
!!                real    (out) :: opacityTransport )
!!
!! DESCRIPTION
!!
!!  Computes absorption, emission and transport opacities for a particular cell.
!!
!! ARGUMENTS
!!
!!   soln              : The solution vector for the cell
!!   ngrp              : The energy group number
!!   opacityAbsorption : the absorption opacity (in 1/cm)
!!   opacityEmission   : the emission opacity (in 1/cm)
!!   opacityTransport  : the transport opacity (in 1/cm)
!!
!!***
subroutine Opacity (soln, ngrp, opacityAbsorption, opacityEmission, opacityTransport)

  use Driver_interface,  ONLY : Driver_abortFlash

  use Opacity_data,      ONLY : op_useOpacity,               &
                                op_ignoreLowTemp,            &
                                op_absorptionKind,           &
                                op_emissionKind,             &
                                op_transportKind,            &
                                op_ionNumberDensity,         &
                                op_totalSpecies,             &
                                op_nEnergyGroups,            &
                                op_cellMassDensity,          &
                                op_cellTemperature,          &
                                op_cellSpeciesMassFractions, &
                                op_speciesLowTempCutoff

  use op_tabulatedData,  ONLY : op_speciesMinTempPATable, &
                                op_speciesMinTempPETable, &
                                op_speciesMinTempROTable

  use op_interface,      ONLY : op_getSpeciesConstantOpacities,  &
                                op_getSpeciesLowTempOpacities,   &
                                op_getSpeciesTabulatedOpacities, &
                                op_computeMultispeciesOpacities, &
                                op_computeIonNumberDensities
  implicit none
  
#include "Flash.h"  
#include "constants.h"
#include "Opacity.h"
  
  integer, intent (in)  :: ngrp
  real,    intent (out) :: opacityAbsorption
  real,    intent (out) :: opacityEmission
  real,    intent (out) :: opacityTransport

  real,    intent (in), dimension (:) :: soln

  logical :: needConstant
  logical :: needLowTemp
  logical :: needLowTempTable
  logical :: needPALowTemp
  logical :: needPELowTemp
  logical :: needROLowTemp
  logical :: needPATable
  logical :: needPETable
  logical :: needROTable
  logical :: needTable
  logical :: needConstcm2g

  integer :: species

  real    :: speciesLowTempCutoff
  real    :: speciesMinTempPATable
  real    :: speciesMinTempPETable
  real    :: speciesMinTempROTable
  real    :: speciesTemperature
  real    :: speciesDensity
!
!
!   ...Check, if opacities need to be evaluated at all.
!
!
  if (.not.op_useOpacity) then
       return
  end if
!
!
!   ...Check the energy group and abort, if out of range.
!
!
  if ( (ngrp < 1) .or. (ngrp > op_nEnergyGroups) ) then
      call Driver_abortFlash ('[Opacity] ERROR: no corresponding energy group')
  end if
!
!
!   ...Store all cell specific properties into internal data.
!
!
  op_cellMassDensity = soln (DENS_VAR)
  op_cellTemperature = soln (TELE_VAR)

  do species = 1,op_totalSpecies
     op_cellSpeciesMassFractions (species) = soln (SPECIES_BEGIN - 1 + species)
  end do
!
!
!   ...Generate the ion number densities for all the species in current cell.
!
!
  call op_computeIonNumberDensities ()
!
!
!   ...Loop over each species and decide between constant, low temperature (if not
!      rejected by user) or tabulated opacities. Extract and store the individual
!      species opacities.
!
!
  do species = 1,op_totalSpecies

     needConstant =      (op_absorptionKind (species) == OP_CONSTANT) &
                    .or. (op_emissionKind   (species) == OP_CONSTANT) &
                    .or. (op_transportKind  (species) == OP_CONSTANT)

     needConstcm2g =      (op_absorptionKind (species) == OP_CONSTCM2G) &
                     .or. (op_emissionKind   (species) == OP_CONSTCM2G) &
                     .or. (op_transportKind  (species) == OP_CONSTCM2G)

     needPATable  =      (op_absorptionKind (species) == OP_TABULAR_PA) &
                    .or. (op_emissionKind   (species) == OP_TABULAR_PA) &
                    .or. (op_transportKind  (species) == OP_TABULAR_PA)

     needPETable  =      (op_absorptionKind (species) == OP_TABULAR_PE) &
                    .or. (op_emissionKind   (species) == OP_TABULAR_PE) &
                    .or. (op_transportKind  (species) == OP_TABULAR_PE)

     needROTable  =      (op_absorptionKind (species) == OP_TABULAR_RO) &
                    .or. (op_emissionKind   (species) == OP_TABULAR_RO) &
                    .or. (op_transportKind  (species) == OP_TABULAR_RO)

     needTable    =       needPATable &
                    .or.  needPETable &
                    .or.  needROTable

     if (needConstant) then
         call op_getSpeciesConstantOpacities (species, op_cellMassDensity)
     end if

     if (needConstcm2g) then
         call op_getSpeciesConstcm2gOpacities (species)
     end if

     if (needTable) then

         speciesTemperature    = op_cellTemperature

         if (.not.op_ignoreLowTemp) then

              speciesLowTempCutoff  = op_speciesLowTempCutoff  (species)
              speciesMinTempPATable = op_speciesMinTempPATable (species)
              speciesMinTempPETable = op_speciesMinTempPETable (species)
              speciesMinTempROTable = op_speciesMinTempROTable (species)

              needPALowTemp    = needPATable .and. (speciesTemperature < speciesMinTempPATable)
              needPELowTemp    = needPETable .and. (speciesTemperature < speciesMinTempPETable)
              needROLowTemp    = needROTable .and. (speciesTemperature < speciesMinTempROTable)

              needLowTempTable = needPALowTemp .or. needPELowTemp .or. needROLowTemp
              needLowTemp      = needLowTempTable .or. (speciesTemperature < speciesLowTempCutoff)

              if (needLowTemp) then
                  call op_getSpeciesLowTempOpacities (species,speciesTemperature,ngrp)
                  CYCLE
              end if

         end if

         speciesDensity = op_ionNumberDensity (species)

         call op_getSpeciesTabulatedOpacities (species,            &
                                               speciesTemperature, &
                                               speciesDensity,     &
                                               ngrp,               &
                                               needPATable,        &
                                               needPETable,        &
                                               needROTable         )
     end if

  end do
!
!
!   ...Calculate the total cell opacities.
!
!
  call op_computeMultispeciesOpacities (opacityAbsorption, &
                                        opacityEmission,   &
                                        opacityTransport   )
!
!
!   ...Ready!
!
!
  return
end subroutine Opacity

