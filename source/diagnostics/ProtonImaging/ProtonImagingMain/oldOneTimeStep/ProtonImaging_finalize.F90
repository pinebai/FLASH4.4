!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonImaging_finalize
!!
!! NAME
!!  
!!  ProtonImaging_finalize
!!
!! SYNOPSIS
!! 
!!  call ProtonImaging_finalize ()
!!
!! DESCRIPTION
!!
!!  Finalizes the ProtonImaging unit. Deallocates all arrays and closes all
!!  proton imaging detector files.
!!
!! ARGUMENTS
!!
!!***
subroutine ProtonImaging_finalize ()

  use ProtonImaging_data, ONLY: pi_beams,                    &
                                pi_cellBfield,               &
                                pi_cellEdgesX,               &
                                pi_cellEdgesY,               &
                                pi_cellEdgesZ,               &
                                pi_cellEfield,               &
                                pi_detectorFilesID,          &
                                pi_detectors,                &
                                pi_IOprotonPointCount,       &
                                pi_IOprotonPoints,           &
                                pi_IOprotonTags,             &
                                pi_mpiScreenProtonType,      &
                                pi_protonBlockID,            &
                                pi_protonNumberBlockID,      &
                                pi_protons,                  &
                                pi_randomNumberSeedArray,    &
                                pi_screenProtonBucketCount,  &
                                pi_screenProtonBuckets,      &
                                pi_screenProtonCountOffsets, &
                                pi_screenProtonCountProcs,   &
                                pi_screenProtons,            &
                                pi_xCircle,                  &
                                pi_yCircle,                  &
                                pi_xSphere,                  &
                                pi_ySphere,                  &
                                pi_zSphere

  use pi_interface,      ONLY : pi_statisticalFinalize

  implicit none

  integer :: error
!
!
!   ...Deallocate everything that might still be allocated.
!
!
  if (allocated (pi_beams)                      ) deallocate (pi_beams)
  if (allocated (pi_cellBfield)                 ) deallocate (pi_cellBfield)
  if (allocated (pi_cellEdgesX)                 ) deallocate (pi_cellEdgesX)
  if (allocated (pi_cellEdgesY)                 ) deallocate (pi_cellEdgesY)
  if (allocated (pi_cellEdgesZ)                 ) deallocate (pi_cellEdgesZ)
  if (allocated (pi_cellEfield)                 ) deallocate (pi_cellEfield)
  if (allocated (pi_detectorFilesID)            ) deallocate (pi_detectorFilesID)
  if (allocated (pi_detectors)                  ) deallocate (pi_detectors)
  if (allocated (pi_IOprotonPointCount)         ) deallocate (pi_IOprotonPointCount)
  if (allocated (pi_IOprotonPoints)             ) deallocate (pi_IOprotonPoints)
  if (allocated (pi_IOprotonTags)               ) deallocate (pi_IOprotonTags)
  if (allocated (pi_protonBlockID)              ) deallocate (pi_protonBlockID)
  if (allocated (pi_protonNumberBlockID)        ) deallocate (pi_protonNumberBlockID)
  if (allocated (pi_protons)                    ) deallocate (pi_protons)
  if (allocated (pi_screenProtonBucketCount)    ) deallocate (pi_screenProtonBucketCount)
  if (allocated (pi_screenProtonBuckets)        ) deallocate (pi_screenProtonBuckets)
  if (allocated (pi_screenProtonCountOffsets)   ) deallocate (pi_screenProtonCountOffsets)
  if (allocated (pi_screenProtonCountProcs)     ) deallocate (pi_screenProtonCountProcs)
  if (allocated (pi_screenProtons)              ) deallocate (pi_screenProtons)
  if (allocated (pi_xCircle)                    ) deallocate (pi_xCircle)
  if (allocated (pi_yCircle)                    ) deallocate (pi_yCircle)
  if (allocated (pi_xSphere)                    ) deallocate (pi_xSphere)
  if (allocated (pi_ySphere)                    ) deallocate (pi_ySphere)
  if (allocated (pi_zSphere)                    ) deallocate (pi_zSphere)
!
!
!   ...Finalize the statistical environment.
!
!
  call pi_statisticalFinalize ()
!
!
!   ...Free the screen proton datatype handle.
!
!
  call MPI_Type_Free (pi_mpiScreenProtonType, error)
!
!
!    ...Ready!
!
!  
  return
end subroutine ProtonImaging_finalize
