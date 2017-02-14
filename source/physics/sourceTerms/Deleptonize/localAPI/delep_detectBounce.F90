!!****if* source/physics/sourceTerms/Deleptonize/localAPI/delep_detectBounce
!!
!! NAME
!!  
!!  delep_detectBounce 
!!
!!
!! SYNOPSIS
!! 
!!  call delep_detectBounce (integer(IN) :: blockCount,
!!             integer(IN) :: blockList(blockCount),
!!             real(IN)    :: dt,
!!             real(IN)    :: time)
!!
!!  
!!  
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!  blockCount : number of blocks to operate on
!!  blockList  : list of blocks to operate on
!!  dt         : current timestep
!!  time       : current time
!!
!!***

subroutine delep_detectBounce (blockCount,blockList,dt,time)
!
!==============================================================================
!
#include "Flash.h"
#include "constants.h"

  use Deleptonize_data, ONLY : delep_postBounce, delep_meshComm, delep_meshMe
 
  implicit none
  include "Flash_mpi.h"
  
  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN)::blockList
  real,intent(IN) :: dt,time
  
end subroutine delep_detectBounce
