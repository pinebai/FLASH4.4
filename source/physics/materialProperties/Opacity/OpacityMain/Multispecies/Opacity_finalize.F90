!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/Opacity_finalize
!!
!! NAME
!!
!!  Opacity_finalize
!!
!! SYNOPSIS
!!
!!  call Opacity_finalize ()
!!
!! DESCRIPTION
!!
!!  Cleans up the Opacity unit.
!!
!! ARGUMENTS
!!
!!***
subroutine Opacity_finalize ()

  use Opacity_data,   ONLY : op_useOpacity,               &
                             op_atomName,                 &
                             op_absorptionKind,           &
                             op_emissionKind,             &
                             op_transportKind,            &
                             op_elementNumberofSpecies,   &
                             op_speciesNumberofElements,  &
                             op_element2AtomicNumber,     &
                             op_speciesWeights,           &
                             op_speciesLowTempCutoff,     &
                             op_cellSpeciesMassFractions, &
                             op_ionNumberDensity,         &
                             op_energyGroupBoundaries,    &
                             op_elements2Species,         &
                             op_species2Elements,         &
                             op_species2FractionElements, &
                             op_speciesOpacities

  use op_interface,   ONLY : op_finalizeNumerics,  &
                             op_finalizeConstant,  &
                             op_finalizeConstcm2g, &
                             op_finalizeTabulated, &
                             op_finalizeIntegrate, &
                             op_finalizeLowTemp

  implicit none
!
!
!    ...Check, if the opacity unit was used at all.
!
!
  if (.not.op_useOpacity) then
       return
  end if
!
!
!    ...Remove the allocated arrays.
!
!
  if (allocated (op_atomName)                ) deallocate (op_atomName)
  if (allocated (op_absorptionKind)          ) deallocate (op_absorptionKind)
  if (allocated (op_emissionKind)            ) deallocate (op_emissionKind)
  if (allocated (op_transportKind)           ) deallocate (op_transportKind)
  if (allocated (op_elementNumberofSpecies)  ) deallocate (op_elementNumberofSpecies)
  if (allocated (op_speciesNumberofElements) ) deallocate (op_speciesNumberofElements)
  if (allocated (op_element2AtomicNumber)    ) deallocate (op_element2AtomicNumber)
  if (allocated (op_speciesWeights)          ) deallocate (op_speciesWeights)
  if (allocated (op_speciesLowTempCutoff)    ) deallocate (op_speciesLowTempCutoff)
  if (allocated (op_cellSpeciesMassFractions)) deallocate (op_cellSpeciesMassFractions)
  if (allocated (op_ionNumberDensity)        ) deallocate (op_ionNumberDensity)
  if (allocated (op_energyGroupBoundaries)   ) deallocate (op_energyGroupBoundaries)
  if (allocated (op_elements2Species)        ) deallocate (op_elements2Species)
  if (allocated (op_species2Elements)        ) deallocate (op_species2Elements)
  if (allocated (op_species2FractionElements)) deallocate (op_species2FractionElements)
  if (allocated (op_speciesOpacities)        ) deallocate (op_speciesOpacities)
!
!
!    ...Remove the allocated arrays in the subunits.
!
!
  call op_finalizeNumerics  ()
  call op_finalizeConstant  ()
  call op_finalizeConstcm2g ()
  call op_finalizeTabulated ()
  call op_finalizeIntegrate ()
  call op_finalizeLowTemp   ()
!
!
!    ...Ready!
!
!
end subroutine Opacity_finalize
