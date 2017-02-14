!!****f* source/physics/Hydro/Hydro_computeDt
!!
!! NAME
!!
!!  Hydro_computeDt
!!
!!
!! SYNOPSIS
!!
!!  Hydro_computeDt(integer(IN) :: blockID, 
!!  
!!                  real(IN) ::  x(:), 
!!                  real(IN) :: dx(:), 
!!                  real(IN) :: uxgrid(:),
!!                  real(IN) :: y(:), 
!!                  real(IN) :: dy(:), 
!!                  real(IN) :: uygrid(:), 
!!                  real(IN) ::  z(:), 
!!                  real(IN) :: dz(:), 
!!                  real(IN) :: uzgrid(:), 
!!                  real,pointer ::  solnData(:,:,:,:),   
!!                  real,(INOUT) ::   dtCheck, 
!!                  integer(INOUT) :: dtMinLoc(:),
!!                  real(INOUT), optional :: extraInfo)
!!
!! DESCRIPTION
!!
!!  Computes the timestep limiter for the hydrodynamical solver.  For pure
!!  hydrodynamics, the Courant-Fredrichs-Lewy criterion is used.  The sound
!!  speed is computed and together with the velocities, is used to constrain
!!  the timestep such that no information can propagate more than one zone
!!  per timestep.
!!
!!
!! ARGUMENTS
!!
!!  blockID --      local block ID
!!   
!!  x, y, z --      coordinates
!!  dx, dy, dz --   deltas in each {x, y z} directions
!!  uxgrid, uygrid, uzgrid-- velocity of grid expansion in {x, y z} directions
!!  solnData --     the physical, solution data from grid
!!  dtCheck --     variable to hold timestep constraint
!!  dtMinLoc(5) -- array to hold location of cell responsible for minimum dt:
!!                 dtMinLoc(1) = i index
!!                 dtMinLoc(2) = j index
!!                 dtMinLoc(3) = k index
!!                 dtMinLoc(4) = blockID
!!                 dtMinLoc(5) = hy_meshMe
!!  extraInfo   -  Driver_computeDt can provide extra info to the caller
!!                 using this argument.
!!
!!***

subroutine Hydro_computeDt (blockID,  &
                           x, dx, uxgrid, &
                           y, dy, uygrid, &
                           z, dz, uzgrid, &
                           blkLimits,blkLimitsGC,        &
                           solnData,   &
                           dtCheck, dtMinLoc, &
                           extraInfo)
     
  
#include "Flash.h"
#include "constants.h"

  implicit none


  integer, intent(IN) :: blockID 
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real,INTENT(INOUT)    :: dtCheck
  integer,INTENT(INOUT)    :: dtMinLoc(5)
  real, pointer, dimension(:,:,:,:) :: solnData
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: x, dx, uxgrid
  real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: y, dy, uygrid
  real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: z, dz, uzgrid
#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: x, dx, uxgrid
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: y, dy, uygrid
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: z, dz, uzgrid
#endif
  real, OPTIONAL,intent(INOUT) :: extraInfo

  return
end subroutine Hydro_computeDt

