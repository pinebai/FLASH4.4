!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonsCreate/pi_createProtons
!!
!! NAME
!!
!!  pi_createProtons
!!
!! SYNOPSIS
!!
!!  call pi_createProtons (integer, intent (in) :: blockCount, 
!!                         integer, intent (in) :: blockList (:),
!!                         real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Generates protons and places them in their initial blocks. This routine calls the
!!  appropriate subroutines according to the domain grid geometry specified.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  timeSimulation : current simulation time
!!
!!***

subroutine pi_createProtons (blockCount, blockList, timeSimulation)

  use Driver_interface,    ONLY : Driver_abortFlash

  use pi_interface,        ONLY : pi_createProtons3DRec

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

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeSimulation
!
!
!     ...Select the appropriate subroutine.
!
!
  select case (pi_gridGeometry)

    case (GRID_3DCARTESIAN)
      call pi_createProtons3DRec (blockCount, blockList, timeSimulation)
    case default
      call Driver_abortFlash ('[pi_createProtons] ERROR: unsupported/unknown geometry')

  end select
!
!
!     ...Inform the master processor how many protons were created globally.
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
     write (pi_monitorFileUnit,'(a,i8,a)') ' created ', pi_globalProtonCount,' Domain Beam   Protons'
 end if
!
!
!    ...Ready!
!
!
  return
end subroutine pi_createProtons
