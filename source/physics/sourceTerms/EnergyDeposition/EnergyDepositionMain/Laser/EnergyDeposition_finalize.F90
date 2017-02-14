!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/EnergyDeposition_finalize
!!
!! NAME
!!  
!!  EnergyDeposition_finalize
!!
!! SYNOPSIS
!! 
!!  call EnergyDeposition_finalize ()
!!
!! DESCRIPTION
!!
!!  Finalizes the EnergyDeposition unit.
!!
!! ARGUMENTS
!!
!!***
subroutine EnergyDeposition_finalize ()

  use EnergyDeposition_data, ONLY: ed_beams,                    &
                                   ed_cellCenters,              &
                                   ed_cellDensity,              &
                                   ed_cellEdges,                &
                                   ed_cellGradNele,             &
                                   ed_cellGradTele,             &
                                   ed_cellNele,                 &
                                   ed_cellTele,                 &
                                   ed_cellVolume,               &
                                   ed_cellZbar,                 &
                                   ed_energyProfileFileName,    &
                                   ed_energyProfileFileUnit,    &
                                   ed_gcMask,                   &
                                   ed_laserIONumberOfPositions, &
                                   ed_laserIORayPositions,      &
                                   ed_laserIORayPower,          &
                                   ed_laserIORayTags,           &
                                   ed_pulseNumberOfSections,    &
                                   ed_pulses,                   &
                                   ed_rayBlockID,               &
                                   ed_rayNumberBlockID,         &
                                   ed_rays,                     &
                                   ed_raysSaved

  use ed_commInterface,      ONLY : ed_commFinalize
  implicit none

  logical :: fileOpen
!
!
!   ...Deallocate everything that might still be allocated.
!
!
  if (allocated (ed_beams)                   ) deallocate (ed_beams)
  if (allocated (ed_cellCenters)             ) deallocate (ed_cellCenters)
  if (allocated (ed_cellDensity)             ) deallocate (ed_cellDensity)
  if (allocated (ed_cellEdges)               ) deallocate (ed_cellEdges)
  if (allocated (ed_cellGradNele)            ) deallocate (ed_cellGradNele)
  if (allocated (ed_cellGradTele)            ) deallocate (ed_cellGradTele)
  if (allocated (ed_cellNele)                ) deallocate (ed_cellNele)
  if (allocated (ed_cellTele)                ) deallocate (ed_cellTele)
  if (allocated (ed_cellVolume)              ) deallocate (ed_cellVolume)
  if (allocated (ed_cellZbar)                ) deallocate (ed_cellZbar)
  if (allocated (ed_laserIONumberOfPositions)) deallocate (ed_laserIONumberOfPositions)
  if (allocated (ed_laserIORayPositions)     ) deallocate (ed_laserIORayPositions)
  if (allocated (ed_laserIORayPower)         ) deallocate (ed_laserIORayPower)
  if (allocated (ed_laserIORayTags)          ) deallocate (ed_laserIORayTags)
  if (allocated (ed_pulseNumberOfSections)   ) deallocate (ed_pulseNumberOfSections)
  if (allocated (ed_pulses)                  ) deallocate (ed_pulses)
  if (allocated (ed_rayBlockID)              ) deallocate (ed_rayBlockID)
  if (allocated (ed_rayNumberBlockID)        ) deallocate (ed_rayNumberBlockID)
  if (allocated (ed_rays)                    ) deallocate (ed_rays)
  if (allocated (ed_raysSaved)               ) deallocate (ed_raysSaved)

  if (allocated (ed_gcMask)                  ) deallocate (ed_gcMask)
!
!
!   ...Close the energy profile printout file (if open).
!
!
  inquire (unit   = ed_energyProfileFileUnit, &
           name   = ed_energyProfileFileName, &
           opened = fileOpen                  )

  if (fileOpen) then
      close (ed_energyProfileFileUnit)
  end if
!
!
!    ...Ready!
!
!  
  call ed_commFinalize()
  return
end subroutine EnergyDeposition_finalize
