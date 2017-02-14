!!****f* source/physics/Gravity/Gravity_accelListOfBlocks
!!
!! NAME
!!
!!  Gravity_accelListOfBlocks  
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelListOfBlocks(integer(IN) :: blockCount,
!!                         integer(:)(IN) :: blockList,
!!                         integer(:)(IN) :: component,
!!                         integer(IN)    :: accelIndex,
!!                         integer(IN),optional    :: potentialIndex)
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration on all mesh blocks.  Either a single component
!!  of the acceleration or all three can be computed.
!!
!! ARGUMENTS
!!
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which to calculate acceleration.
!!   component    : The directional component of the acceleration to compute.
!!                  Permitted values are IAXIS, JAXIS, KAXIS, and ALLDIR.
!!                  These constants are defined in constants.h.
!!   accelIndex   : grid variable to store the acceleration, e.g., GRAC
!!   potentialIndex :   Variable to take as potential if present.  If not
!!                      present, GPOT is assumed by default.
!!
!! NOTES
!!
!!   This routine can be used as a wrapper to Gravity_accelOneRow.  Each implementation
!!   of the Gravity unit has a version of Gravity_accelOneRow, but this wrapper remains
!!   constant.
!!
!!***

subroutine Gravity_accelListOfBlocks (blockCount,blockList,component, &
     accelIndex, potentialIndex)

!==============================================================================

  implicit none

  integer,intent(IN)                      :: blockCount
  integer,dimension(blockCount), intent(IN)     :: blockList
  integer, INTENT(in) ::  component
  integer, INTENT(in) ::  accelIndex
  integer,intent(IN),optional :: potentialIndex
!==============================================================================


  return
end subroutine Gravity_accelListOfBlocks
