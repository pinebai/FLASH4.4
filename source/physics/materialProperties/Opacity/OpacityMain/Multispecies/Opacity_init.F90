!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/Opacity_init
!!
!! NAME
!!
!!  Opacity_init
!!
!! SYNOPSIS
!!
!!  call Opacity_init ()
!!
!! DESCRIPTION
!!
!!  Initializes the data for the Opacity unit.
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!
!!  useOpacity
!!  rt_useMGD
!!  opacity_writeOpacityInfo
!!  opacity_ignoreLowTemp
!!  rt_mgdNumGroups
!!
!!***
subroutine Opacity_init ()

  use Opacity_data

  use op_numericsData,             ONLY : zero

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use Driver_interface,            ONLY : Driver_abortFlash, &
                                          Driver_getMype

  use op_interface,                ONLY : op_initNumerics,                &
                                          op_initConstant,                &
                                          op_initConstcm2g,               &
                                          op_initTabulated,               &
                                          op_initIntegrate,               &
                                          op_initLowTemp,                 &
                                          op_setAtomNames,                &
                                          op_setAtomWeights,              &
                                          op_setEnergyGroupBoundaries,    &
                                          op_setLowTemperatureBoundaries, &
                                          op_setSpeciesElementsData,      &
                                          op_writeOpacity

  use Timers_interface, ONLY : Timers_start, Timers_stop 
  implicit none

# include "constants.h"
# include "Opacity.h"
  
  integer :: species

  call Driver_getMype(GLOBAL_COMM,op_globalMe)

!
!
!    ...Check, if the opacity unit is needed at all.
!
!
  call RuntimeParameters_get ("useOpacity",   op_useOpacity)

  if (.not.op_useOpacity) then
       return
  end if

  call RuntimeParameters_get("rt_useMGD", op_useMGD)

  !!if(.not. op_useMGD) then
  !!   call Driver_abortFlash("[Opacity_init] useOpacity is .true. but rt_useMGD is .false.")
  !!end if

  call Timers_start("Opacity_init")
!
!
!    ...Get the write and low temperature action parameters.
!
!
  call RuntimeParameters_get ("opacity_writeOpacityInfo",   op_writeOpacityInfo)
  call RuntimeParameters_get ("opacity_ignoreLowTemp",      op_ignoreLowTemp)
!
!
!    ...Initialize the numeric subunit to have all numerical constants available.
!
!
  call op_initNumerics ()
!
!
!    ...Sets all the data related to species and elements.
!
!
  call op_setSpeciesElementsData ()
!
!
!    ...Set internal parameters:
!
!          i) # of radiation energy groups
!
!
  call RuntimeParameters_get("rt_mgdNumGroups", op_nEnergyGroups)

  if (op_nEnergyGroups < 1) then
      call Driver_abortFlash ('[Opacity_init] ERROR: no energy groups found')
  end if
!
!
!    ...Allocate all needed arrays.
!
!
  allocate (op_atomName                 (1:100))
  allocate (op_atomWeight               (1:100))
  allocate (op_absorptionKind           (1:op_totalSpecies))
  allocate (op_emissionKind             (1:op_totalSpecies))
  allocate (op_transportKind            (1:op_totalSpecies))
  allocate (op_cellSpeciesMassFractions (1:op_totalSpecies))
  allocate (op_ionNumberDensity         (1:op_totalSpecies))
  allocate (op_energyGroupBoundaries    (1:op_nEnergyGroups + 1))

  allocate (op_speciesOpacities         (ABSORPTION:TRANSPORT , 1:op_totalSpecies))
!
!
!    ...Inititalize those arrays with zeros, which will contain the same data throughout the run.
!
!
  op_absorptionKind = 0
  op_emissionKind   = 0
  op_transportKind  = 0
!
!
!    ...Set the names and weights of all possible atoms and the energy group boundaries.
!
!
  call op_setAtomNames   ()
  call op_setAtomWeights ()
  call op_setEnergyGroupBoundaries ()
!
!
!    ...Inititalize the other subunits. The low temperature subunit can only be inititalized
!       after the low temperature boundaries have been set, which in turn requires
!       initialization of the tabulated subunit, so the order matters.
!
!
  call op_initConstant  ()
  call op_initConstcm2g ()
  call op_initTabulated ()
  call op_initIntegrate ()

  if (.not.op_ignoreLowTemp) then
       call op_setLowTemperatureBoundaries ()
       call op_initLowTemp   ()
  end if

  ! Check for errors:
  

  if(  any(op_absorptionKind == 0) .or. &
       any(op_emissionKind   == 0) .or. &
       any(op_transportKind  == 0) ) then
     ! If you are here, then there one of the opacity models for a
     ! species is not specified properly. Make sure the opacity model
     ! for each species is specified in the flash.par

     if (op_globalMe == MASTER_PE) then
        do species = 1, op_totalSpecies
           print '(4i6)', species, & 
                op_absorptionKind(species), &
                op_emissionKind(species), &
                op_transportKind(species)
        end do
     end if

     call Driver_abortFlash("[opacity_init] bad opacity kind")
  end if

!
!
!    ...Write the opacity unit data to file (if requested).
!
!
  if (op_writeOpacityInfo) then
      call op_writeOpacity ()
  end if
!
!
!    ...Ready!
!
!
  call Timers_stop("Opacity_init")
  return
end subroutine Opacity_init
