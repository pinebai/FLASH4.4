!!****if* source/physics/Gravity/GravityMain/Constant/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow 
!!
!! SYNOPSIS
!!
!!  call Gravity_accelOneRow(integer(IN)  :: pos(2),
!!                           integer(IN)  :: sweepDir,
!!                           integer(IN)  :: blockID,
!!                           integer(IN)  :: numCells,
!!                           real(INOUT)  :: grav(numCells),
!!                           integer(IN),optional :: potentialIndex,
!!                           integer(IN),optional :: extraAccelVars(MDIM))
!!
!! DESCRIPTION
!!
!!  This routine computes the gravitational acceleration for a row
!!  of cells in a specified direction in a given block.
!!
!! ARGUMENTS
!!
!!  pos      :  Row indices transverse to the sweep direction
!!  sweepDir :    The sweep direction:  allowed values are 
!!              SWEEP_X, SWEEP_Y, and SWEEP_Z. These values are defined
!!              in constants.h.
!!  blockID  :  The local identifier of the block to work on
!!  numCells :  Number of cells to update in grav()
!!  grav     :   Array to receive result
!!  potentialIndex :  optional, not applicable in constant gravity
!!  extraAccelVars :  optional, ignored in this implementation
!!                    
!! 
!!***

subroutine Gravity_accelOneRow (pos, sweepDir, blockID, numCells, grav, &
                                potentialIndex, extraAccelVars)

!==============================================================================
!

  use Gravity_data, ONLY : useGravity, grv_vector

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, dimension(2), intent(IN) :: pos
  integer,INTENT(in) :: sweepDir
  integer,INTENT(in) :: blockID
  integer,INTENT(in) :: numCells
  real,dimension(numCells),INTENT(inout) :: grav
  integer,intent(IN),optional :: potentialIndex
  integer,intent(IN),OPTIONAL :: extraAccelVars(MDIM)

  real :: grv_val


  if (useGravity) then
     grv_val = grv_vector(sweepDir)
  
     grav(1:numCells) = grv_val
  end if


!
!==============================================================================
!
  return
end subroutine Gravity_accelOneRow
