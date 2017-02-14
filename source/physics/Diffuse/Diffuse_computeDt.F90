!!****f* source/physics/Diffuse/Diffuse_computeDt
!!
!! NAME
!!  
!!  Diffuse_computeDt
!!
!!
!! SYNOPSIS
!! 
!!  Diffuse_computeDt ( integer(IN) : blockID, 
!!                  real(IN):  xCenter(:), 
!!                  real(IN):  xLeft(:), 
!!                  real(IN):  xRight(:), 
!!                  real(IN): dx(:), 
!!                  real(IN): uxgrid(:),
!!                  real(IN):  yCenter(:), 
!!                  real(IN):  yLeft(:), 
!!                  real(IN):  yRight(:), 
!!                  real(IN): dy(:), 
!!                  real(IN): uygrid(:), 
!!                  real(IN):  zCenter(:), 
!!                  real(IN):  zLeft(:), 
!!                  real(IN):  zRight(:), 
!!                  real(IN): dz(:), 
!!                  real(IN): uzgrid(:), 
!!                  real,pointer :  solnData(:,:,:,:),   
!!                  real,(INOUT):   dt_check, 
!!                  integer(INOUT): dt_minloc(:) )
!!  
!! DESCRIPTION
!!
!!  Computes the timestep limiter for diffusion source term solver.
!! 
!!  The current implementation may be very conservative, especially with
!!  respect to the viscosity term.  Users may want to change the implementation
!!  to be less conservative, and/or tweak the time step by tweaking the
!!  dt_diff_factor runtime parameter.
!!
!! ARGUMENTS
!!
!!  blockID        local block ID
!!  xCenter         X coordinates at the center of the cell
!!  xLeft           X coordinates at the left edge of the cell
!!  xRight          X coordinates at the right edge of the cell
!!  yCenter         Y coordinates at the center of the cell
!!  yLeft           Y coordinates at the left edge of the cell
!!  yRight          Y coordinates at the right edge of the cell
!!  zCenter         Z coordinates at the center of the cell
!!  zLeft           Z coordinates at the left edge of the cell
!!  zRight          Z coordinates at the right edge of the cell
!!  d*              deltas in each {*=x, y z} directions
!!  u*grid          velocity of grid expansion in {*=x, y z} directions
!!  solnData        the physical, solution data from grid
!!  dt_check        variable to hold timestep constraint
!!  dt_minloc(5)    array to hold limiting zone info:  zone indices
!!
!!***


subroutine Diffuse_computeDt (blockID, &
                              xCenter,xLeft,xRight, dx, uxgrid, &
                              yCenter,yLeft,yRight, dy, uygrid, &
                              zCenter,zLeft,zRight, dz, uzgrid, &
                              blkLimits,blkLimitsGC,        &
                              solnData,   &
                              dt_check, dt_minloc )


  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(IN) :: blockID
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real,INTENT(INOUT)    :: dt_check
  integer,INTENT(INOUT)    :: dt_minloc(5)
  real, pointer :: solnData(:,:,:,:) 

  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: &
       xCenter,xLeft,xRight
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: &
       yCenter,yLeft,yRight
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) ::&
        zCenter,zLeft,zRight
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: &
       dx, uxgrid
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: &
       dy, uygrid
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: &
       dz, uzgrid
 
  return
end subroutine Diffuse_computeDt

