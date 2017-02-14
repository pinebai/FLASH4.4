!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/gr_ptMoveMappedData
!!
!! NAME
!!  gr_ptMoveMappedData
!!
!! SYNOPSIS
!!
!!  gr_ptMoveMappedData(integer,intent(IN) :: varGrid, &
!!                       integer,intent(IN) :: bufferSize, &
!!                       real,dimension(bufferSize),intent(INOUT) :: sendBuf, &
!!                       integer,intent(INOUT) :: sendCount, &
!!                       real,dimension(bufferSize),intent(INOUT) :: recvBuf)
!!
!! DESCRIPTION
!!
!! Routine which manages the communication of smeared grid cells between processors. 
!!
!! This is a STUB.
!! 
!! ARGUMENTS
!!               varGrid:   Index of gridded variable to receive interpolated
!!                              quantity
!!               bufferSize:  The size of the sendBuf and recvBuf arrays
!!               sendBuf:  An array used to store data intended for another processor
!!               sendCount:  The number of data elementes to be sent to another
!!                           processor
!!               recvBuf:  An array containing the data just receieved from another 
!!                         processor
!!
!! PARAMETERS
!! 
!!***

subroutine gr_ptMoveMappedData(varGrid,bufferSize,sendBuf,sendCount,recvBuf)

  use Grid_data, ONLY : gr_globalMe, gr_meshNumProcs
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY: Logfile_stampMessage
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_ptInterface, ONLY : gr_ptPackUnpackData

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "gr_ptMapToMesh.h"

  integer,intent(IN) :: varGrid
  integer,intent(IN) :: bufferSize
  real,dimension(bufferSize),intent(INOUT) :: sendBuf
  integer,intent(INOUT) :: sendCount
  real,dimension(bufferSize),intent(INOUT) :: recvBuf

end subroutine gr_ptMoveMappedData
