!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDetection/pi_setupScreenProtons
!!
!! NAME
!!
!!  pi_setupScreenProtons
!!
!! SYNOPSIS
!!
!!  call pi_setupScreenProtons ()
!!
!! DESCRIPTION
!!
!!  Sets up the screen protons, which means to allocate the needed screen proton array.
!!  The screen protons are used to accumulate protons on the detector screen(s).
!!  The routine also sets up and commits a mpi type structure corresponding to the
!!  screen protons, which will be used when sending screen protons between processors.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine pi_setupScreenProtons ()

  use ProtonImaging_data,  ONLY : pi_maxProtonCount,      &
                                  pi_mpiScreenProtonType, &
                                  pi_screenProtons
  
  use Driver_interface,    ONLY : Driver_abortFlash

  implicit none

#include "ProtonImaging.h"
 include "Flash_mpi.h"

  integer :: error
!
!
!     ...Allocate the screen proton array.
!
!
  allocate (pi_screenProtons (1:SCREEN_ATTRCOUNT,1:pi_maxProtonCount))
!
!
!     ...Create and commit the screen proton mpi type structure.
!
!
  call MPI_Type_Contiguous (SCREEN_ATTRCOUNT,                   &
                            FLASH_REAL,                         &
                                        pi_mpiScreenProtonType, &
                                        error                   )

  call MPI_Type_Commit (pi_mpiScreenProtonType, error)
!
!
!     ...Ready!
!
!
  return
end subroutine pi_setupScreenProtons
