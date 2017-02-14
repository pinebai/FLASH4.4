!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDisk/pi_readDiskProtons
!!
!! NAME
!!
!!  pi_readDiskProtons
!!
!! SYNOPSIS
!!
!!  call pi_readDiskProtons (integer, intent (in)    :: blockCount,
!!                           integer, intent (in)    :: blockList (:),
!!                           logical, intent (inout) :: moreOnDisk)
!!
!! DESCRIPTION
!!
!!  Reads in a batch of disk protons from the old disk proton file. This routine
!!  calls the appropriate subroutines according to the domain grid geometry specified.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  moreOnDisk     : if true, there are more disk protons on the old disk proton file
!!
!!***

subroutine pi_readDiskProtons (blockCount, blockList, moreOnDisk)

  use Driver_interface,    ONLY : Driver_abortFlash

  use pi_interface,        ONLY : pi_readDiskProtons3DRec

  use ProtonImaging_data,  ONLY : pi_globalComm,        &
                                  pi_globalMe,          &
                                  pi_globalProtonCount, &
                                  pi_gridGeometry,      &
                                  pi_monitorFileUnit,   &
                                  pi_protonCount

  implicit none

#include "constants.h"
#include "Flash.h"
#include "ProtonImaging.h"

  include "Flash_mpi.h"

  integer :: error

  integer, intent (in)    :: blockCount
  integer, intent (in)    :: blockList (1:blockCount)
  logical, intent (inout) :: moreOnDisk
!
!
!     ...Select the appropriate subroutine.
!
!
  select case (pi_gridGeometry)

    case (GRID_3DCARTESIAN)
      call pi_readDiskProtons3DRec (blockCount, blockList, moreOnDisk)
    case default
      call Driver_abortFlash ('[pi_readDiskProtons] ERROR: unsupported/unknown geometry')

  end select
!
!
!     ...Inform the master processor how many protons were read back in globally.
!
!
  call MPI_Reduce (pi_protonCount,       &
                   pi_globalProtonCount, &
                   1,                    &
                   MPI_INTEGER,          &
                   MPI_SUM,              &
                   MASTER_PE,            &
                   pi_globalComm,        &
                   error                 )

 if (pi_globalMe == MASTER_PE) then
     write (pi_monitorFileUnit,'(a,i8,a)') ' loaded  ', pi_globalProtonCount,' Domain Disk   Protons'
 end if
!
!
!    ...Ready!
!
!
  return
end subroutine pi_readDiskProtons
