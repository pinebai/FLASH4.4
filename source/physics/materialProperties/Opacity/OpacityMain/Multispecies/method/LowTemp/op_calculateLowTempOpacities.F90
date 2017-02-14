!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_calculateLowTempOpacities
!!
!! NAME
!!
!!  op_calculateLowTempOpacities
!!
!! SYNOPSIS
!!
!!  call op_calculateLowTempOpacities ()
!!
!! DESCRIPTION
!!
!!  This routine constructs the low temperature Planck and Rosseland energy group opacities for all
!!  temperature grid/species/energy group triple combinations and places them in the corresponding
!!  tables. Currently it considers the photoelectronic (PE) and Klein-Nishina (KN) cross sections only.
!!  The PE data is extracted from the Biggs & Lighthill tables. 
!!
!! ARGUMENTS
!!
!!***
subroutine op_calculateLowTempOpacities ()

  use Opacity_data,     ONLY : op_nEnergyGroups,             &
                               op_totalElements,             &
                               op_totalSpecies,              &
                               op_elementNumberofSpecies,    &
                               op_elements2Species,          &
                               op_species2FractionElements,  &
                               op_energyGroupBoundaries

  use op_lowTempData,   ONLY : op_ignoreKleinNishina,    &
                               op_minLowTempTemperature, &
                               op_maxLowTempTemperature, &
                               op_maxNstepsLowTemp,      &
                               op_tableLowTemp,          &
                               op_PlanckLowTempTables,   &
                               op_RosselandLowTempTables

  use op_numericsData,  ONLY : zero,         &
                               op_Boltzmann, &
                               op_keV2erg,   &
                               op_erg2keV

  use Driver_interface, ONLY : Driver_abortFlash

  use op_interface,     ONLY : op_BiggsGroupOpacity, &
                               op_KleinGroupOpacity

  implicit none

# include "Opacity.h"

  character (len=9) Planck, Rosseland

  integer :: element
  integer :: group
  integer :: n,s
  integer :: nGrid
  integer :: nSpecies
  integer :: species

  real    :: Elower,Eupper
  real    :: fraction
  real    :: opacityBiggsPlanck, opacityBiggsRosseland
  real    :: opacityKleinPlanck, opacityKleinRosseland
  real    :: opacityPlanck, opacityRosseland
  real    :: Temperature
  real    :: Tlow,Thigh,Tstep
!
!
!   ...Initialize the low temperature arrays.
!
!
  op_tableLowTemp           = zero
  op_PlanckLowTempTables    = zero
  op_RosselandLowTempTables = zero
!
!
!   ...Establish the low temperature grid (linear grid with equal temperature difference
!      between grid points).
!
!
  Tlow  = op_minLowTempTemperature
  Thigh = op_maxLowTempTemperature
  nGrid = op_maxNstepsLowTemp
  Tstep = (Thigh - Tlow) / real (nGrid - 1)

  do n = 1,nGrid
     op_tableLowTemp (n) = Tlow + (n - 1) * Tstep
  end do
!
!
!   ...Outer loop over all energy groups, middle loop over all atomic elements,
!      inner loop over all temperatures on the temperature grid.
!
!
  Planck    = "Planck"
  Rosseland = "Rosseland"

  do group = 1,op_nEnergyGroups

     Elower = op_energyGroupBoundaries (group)         ! in eV
     Eupper = op_energyGroupBoundaries (group+1)       ! in eV

     do element = 1,op_totalElements

        nSpecies = op_elementNumberofSpecies (element)

        do n = 1,nGrid

           Temperature = op_tableLowTemp (n)

           call op_BiggsGroupOpacity (Planck,                        &
                                      element,                       &
                                      Temperature,                   &
                                      Elower,                        &
                                      Eupper,                        &
                                               opacityBiggsPlanck    )

           call op_BiggsGroupOpacity (Rosseland,                     &
                                      element,                       &
                                      Temperature,                   &
                                      Elower,                        &
                                      Eupper,                        &
                                               opacityBiggsRosseland )

           if (.not. op_ignoreKleinNishina) then

                call op_KleinGroupOpacity (Planck,                        &
                                           element,                       &
                                           Temperature,                   &
                                           Elower,                        &
                                           Eupper,                        &
                                                    opacityKleinPlanck    )

                call op_KleinGroupOpacity (Rosseland,                     &
                                           element,                       &
                                           Temperature,                   &
                                           Elower,                        &
                                           Eupper,                        &
                                                    opacityKleinRosseland )
           else
                opacityKleinPlanck    = zero
                opacityKleinRosseland = zero
           end if

           opacityPlanck    = opacityBiggsPlanck    + opacityKleinPlanck
           opacityRosseland = opacityBiggsRosseland + opacityKleinRosseland

           do s = 1,nSpecies

              species  = op_elements2Species         (element,s)
              fraction = op_species2FractionElements (species,element)

              op_PlanckLowTempTables    (n,species,group) =   op_PlanckLowTempTables (n,species,group)    &
                                                                  + fraction * opacityPlanck
              op_RosselandLowTempTables (n,species,group) =   op_RosselandLowTempTables (n,species,group) &
                                                                  + fraction * opacityRosseland
           end do

        end do
     end do
  end do
!
!
!   ...Ready! 
!
!
  return
end subroutine op_calculateLowTempOpacities
