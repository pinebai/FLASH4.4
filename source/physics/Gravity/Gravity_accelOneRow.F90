!!****f* source/physics/Gravity/Gravity_accelOneRow
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
!!  grav()   :   Array to receive result
!!  potentialIndex :   Variable # to take as potential if present
!!                    ( this argument is applicable only in self-gravity
!!                      formulated as Poisson's equation)
!!  extraAccelVars      -  if specified,  Variables from which extra accelerations
!!                         are taken. Used to identify the UNK variables
!!                         that contain sink-on-gas accelerations when
!!                         sink particles are used.
!!
!! NOTES
!!
!!  The following applies to Gravity implementations which solve the
!!  Poisson problem for self-gravity:
!!  If certain variables declared by the sink particles inplementation are
!!  declared, it is assumed that sink particles are in use.
!!  The sets of variables to make this determination are
!!    o  those given by extraAccelVars   extraAccelVars if present;
!!    o  {SGXO_VAR, SGYO_VAR, SGZO_VAR}  if potentialIndex is GPOL_VAR;
!!    o  {SGAX_VAR, SGAY_VAR, SGAZ_VAR}  otherwise.
!!  If it is assumed that sink particles are in use, then the acceleration
!!  returned in the grav array will have the appropriate sink particle
!!  acceleration component added to the acceleration computed by differencing
!!  the potential variable given by potentialIndex.
!!***

subroutine Gravity_accelOneRow (pos, sweepDir, blockID, numCells, grav, &
                                potentialIndex, extraAccelVars)

!===============================================================================

  implicit none
#include "constants.h"

  integer, intent(IN) :: sweepDir,blockID,numCells
  integer, dimension(2),INTENT(in) ::pos
  real, dimension(numCells),INTENT(inout) :: grav
  integer,intent(IN),optional :: potentialIndex
  integer, intent(IN),OPTIONAL      :: extraAccelVars(MDIM)
!======================================================================

  return
end subroutine Gravity_accelOneRow
