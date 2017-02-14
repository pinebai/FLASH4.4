!!****if* source/physics/Gravity/GravityMain/Gravity_accelOneBlock
!!
!! NAME
!!
!!  Gravity_accelOneBlock
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelOneBlock(integer, intent(in) :: blockID, 
!!                        integer, intent(in) :: ngcellcomp,
!!                        real(NDIM,:,:,:)),intent(out) :: gvec, 
!!                        integer, intent(in),optional :: potentialIndex)
!!                      
!!                      
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration for this block.  Include ngcell layers outside
!!  block interior.
!!
!!  This routine computes the gravitational acceleration for
!!  zones in a given block. First-order
!!  finite-volume differencing is used everywhere.  It is assumed
!!  here that the requisite number of guard cells have peen appropriately
!!  filled for the variable containting the gravitational potential.
!!
!!  Dean Townsley 2008
!!  Contributed to Flash Center at the University of Chicago 2008
!!
!! ARGUMENTS
!!
!!  blockID            -  The local identifier of the block to work on
!!  gvec(NDIM,:,:,:)   -  Array to receive gravitational acceleration
!!                        as as NDIM-dimensional vector.  It is assumed
!!                        the the space provided is the size of the block
!!                        plus all guard cells.
!!  ngcellcomp         -  Number of layers outside of block interior to
!!                        compute gravity
!!  potentialIndex     -  if specified,  Variable # to take as potential.
!!                        Default is GPOT_VAR for the potential stored in the
!!                        gpot slot of unk, which should correspond to the
!!                        potential at the current timestep.
!!
!!
!!***


subroutine Gravity_accelOneBlock ( blockID, ngcellcomp, gvec, potentialIndex)

  use Driver_interface, only : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(in)             :: blockID,  ngcellcomp
  real, intent(out)               :: gvec(:,:,:,:)
  integer, intent(in),optional    :: potentialIndex

  !==================================================

  gvec = 0.0
  call Driver_abortFlash("Gravity_accelOneBlock not implemented")

  return
   
end subroutine Gravity_accelOneBlock
