!!****if* source/physics/Gravity/GravityMain/PlanePar/grv_accelOneZone
!!
!! NAME
!!
!!  grv_accelOneZone
!!
!! SYNOPSIS
!!
!!  grv_accelOneZone(    real(IN) :: xZone, 
!!                       real(IN) :: yZone,
!!                       real(IN) :: zZone,
!!                       real(:)(OUT) :: gravZone)
!!
!! DESCRIPTION
!!
!!  A standalone function which computes g for a single zone given the
!!  coordinates in xZone,yZone and zZone
!!
!!  The gravitational acceleration is returned as a vector, giving the 
!!  acceleration in each of the coordinate directions
!!
!!  This version is for a plane parallel atmosphere
!!
!! ARGUMENTS
!!
!!  xZone, yZone, zZone:   coordinates of zone
!!  gravZone:              acceleration vector
!!
!!***

subroutine grv_accelOneZone(xZone, yZone, zZone, gravZone)

!==============================================================================

  use Gravity_data, ONLY: grv_ptdirn, grv_ptxpos, grv_newton, grv_ptmass

  implicit none

#include "Flash.h"
#include "constants.h"
        
  real, INTENT(in) :: xZone, yZone, zZone
  real, INTENT(out), DIMENSION(MDIM) :: gravZone
  
  real :: height, iheight
  
  gravZone(:) = 0.0
  
  if (grv_ptdirn == 1) then
     height = xZone - grv_ptxpos
     
  else if (grv_ptdirn == 2) then
     height = yZone - grv_ptxpos
     
  else if (grv_ptdirn == 3) then
     height = zZone - grv_ptxpos

  endif
  
  iheight = 1.0/height
  gravZone(grv_ptdirn) = -grv_newton*grv_ptmass*iheight*iheight
  
  return
end subroutine grv_accelOneZone



