!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhTreeWalkEnd
!!
!! NAME
!!
!!  Gravity_bhTreeWalkEnd
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhTreeWalkEnd()
!!
!! DESCRIPTION
!!
!!  Called at the end of the Tree Walk.
!!
!! ARGUMENTS
!!
!!
!!***

subroutine Gravity_bhTreeWalkEnd()
  use Gravity_data, ONLY : grv_poisson_max, grv_sink_max, &
    grv_meshComm, grv_meshMe
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer :: ierr
  real :: tot_poisson_max, tot_sink_max
  character(len=MAX_STRING_LENGTH) :: strBuff

  call MPI_Reduce(grv_poisson_max,tot_poisson_max,1,FLASH_REAL,MPI_MAX,MASTER_PE,grv_meshComm, ierr)
  call MPI_Reduce(grv_sink_max,tot_sink_max,1,FLASH_REAL,MPI_MAX,MASTER_PE,grv_meshComm, ierr)

  if (grv_meshMe == MASTER_PE) then
    write (strBuff, '("max gacc comp. Poisson: ", e12.5, ", Sinks: ", e12.5)') &
    & tot_poisson_max, tot_sink_max
    call Logfile_stamp( strBuff, "[BHTree]")
  endif

  ! reset vars for measurement of maximum grav. acceleration
  ! it is done here, because these variables are set in Gravity_accelOneRow
  ! which is outside of the tree solver main loop (Grid_solvePoisson)
  grv_poisson_max = 0.0
  grv_sink_max = 0.0

  return
end subroutine Gravity_bhTreeWalkEnd

