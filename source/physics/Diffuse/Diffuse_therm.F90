!!****f* source/physics/Diffuse/Diffuse_therm
!!
!! NAME
!!
!!  Diffuse_therm
!!
!!
!! SYNOPSIS
!!
!! call Diffuse_therm(integer(IN) :: sweepDir, 
!!                    integer(IN) :: igeom, 
!!                    integer(IN) :: blockID,
!!                    integer(IN) :: numCells,
!!                    integer(IN) :: blkLimits(2,MDIM),
!!                    integer(IN) :: blkLimitsGC(2,MDIM),
!!                    real(IN)    :: leftCoords(MAXCELLS or numCells),
!!                    real(IN)    :: rightCoords(MAXCELLS or numCells),
!!                    real(INOUT) :: temp_flx(NFLUXES,:,:,:),
!!                    real(IN)    :: areaLeft(:,:,:))
!!
!!
!! DESCRIPTION
!!
!!  Diffuse_therm alters the energy flux to account for heat losses through
!!  thermal diffusion.  This is an explicit method, so a timestep limiter
!!  will be required.  Stability is guaranteed for 
!! 
!!                           dx**2 
!!                 dt < .5* -------
!!                             D
!!
!!  where D is the diffusion coefficient, related to the isochoric conductivity (sigma)
!!  and specific heat at constant volume (c_v) by
!!
!!               D = sigma/(rho*c_v)
!!
!!  Fluxes in PPM are stored at the zone boundaries, the temperature, etc. are
!!  at the zone centers.  So the flux {F = -sigma * grad(T)} is
!! 
!!                      sigma_i + sigma_{i-1}     T_i - T_{i-1}
!!             F_i = - ----------------------- * --------------- 
!!                               2                      dx
!!
!!  and then the contribution to the energy equation (evolution of rho*E),
!!  where E is the energy/gram is - div(F) or
!! 
!!                        F_{i+1} - F_i
!!             dE_i  = - ---------------      
!!                             dx    
!!
!! 
!!  together, these give
!! 
!!                     sigma_{i+1} + sigma_i     T_{i+1} - T_i
!!            dE_i =  ----------------------- * ---------------   -
!!                               2                   dx**2
!!
!!                     sigma_i + sigma_{i-1}     T_i - T_{i-1}
!!                    ----------------------- * ---------------   -
!!                               2                   dx**2
!!
!!  which reduces to the standard expression for diffusion if sigma is constant.
!!                     
!!
!! NOTES 
!!
!!    This routine is used by the PPM implementation of the Hydro unit.
!!    Other implementations of Hydro, in particular MHD and RHD, may have their
!!    own mechanisms for handling diffusive effects that bypass Diffuse_therm
!!    (and the Diffuse code unit in general), or they may lack support for
!!    diffusive effects.
!!
!!    This routine computes the heat fluxes and adds them to the energy flux
!!    returned from hydro_1d.  It is to be called after calling hydro_1d on a 
!!    block.  The updated energy fluxes are then used in hy_updateSoln to 
!!    produce the updated energy.
!!
!!    Boundary conditions on the temperature are handled through the GridBoundaryConditions
!!    subunit, not in this routine.  Ex: setting a reflecting boundary will reflect the
!!    temperature into the guard cells, which are used in computing the flux
!!    at the block boundaries.
!!
!!
!!
!! ARGUMENTS
!!
!!  sweepDir    -   the current sweep direction
!!
!!  igeom       -   the geometry flag for the current sweep direction
!!
!!  blockID     -   the block number to operate on
!!  numCells    -   the number of Cells along the sweep direction
!!  leftCoords  -   Coordinates of the left edge of the zones
!!  rightCoords -   Coordinates of the right edge of the zones
!!
!!  temp_flx    -   Temporary storage for fluxes in sweep direction
!!  areaLeft    -   Cell face areas at the left (smaller) cell side
!!  blkLimits -   endpoints of block indices without including gcells
!!  blkLimitsGC -   endpoints of block indices including gcells
!!
!! HISTORY
!!
!!  Apparently this subroutine used to be called therm_explicit.
!!***

subroutine Diffuse_therm(sweepDir, igeom, blockID,numCells,blkLimits,blkLimitsGC,&
                           leftCoords,rightCoords,&
                           temp_flx, areaLeft)


  
  implicit none
#include "constants.h"  
#include "Flash.h"
  
  integer, intent(IN) :: sweepDir, igeom, blockID, numCells
  integer, dimension(2,MDIM), intent (IN) :: blkLimitsGC, blkLimits
  
#ifdef FIXEDBLOCKSIZE
  real, intent(INOUT), DIMENSION(NFLUXES,                   &
                               GRID_ILO_GC:GRID_IHI_GC,     &
                               GRID_JLO_GC:GRID_JHI_GC,     &
                               GRID_KLO_GC:GRID_KHI_GC) ::  &
                               temp_flx
  real, intent(IN), DIMENSION(                   &
                               GRID_ILO_GC:GRID_IHI_GC,     &
                               GRID_JLO_GC:GRID_JHI_GC,     &
                               GRID_KLO_GC:GRID_KHI_GC) ::  &
                               areaLeft
  real,intent(IN), DIMENSION(MAXCELLS) :: leftCoords ,rightCoords
#else
  real, intent(INOUT), DIMENSION(NFLUXES,                   &
                               blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
                               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
                               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  &
                               temp_flx
  real, intent(IN), DIMENSION(                   &
                               blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
                               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
                               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  &
                               areaLeft
  
  real, intent(IN), DIMENSION(numCells) :: leftCoords ,rightCoords
#endif

  return
end subroutine Diffuse_therm

      

