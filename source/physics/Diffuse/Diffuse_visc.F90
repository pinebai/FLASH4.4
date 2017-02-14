!!****f* source/physics/Diffuse/Diffuse_visc
!!
!! NAME
!!
!!  Diffuse_visc
!!
!!
!! SYNOPSIS
!!
!! call Diffuse_visc (integer(IN) :: sweepDir, 
!!                    integer(IN) :: igeom, 
!!                    integer(IN) :: blockID,
!!                    integer(IN) :: numCells,
!!                    integer(IN) :: blkLimits(2,MDIM),
!!                    integer(IN) :: blkLimitsGC(2,MDIM),
!!                    real(IN)    :: leftCoords(MAXCELLS or numCells),
!!                    real(IN)    :: rightCoords(MAXCELLS or numCells),
!!                    real(INOUT) :: temp_flx(NFLUXES,:,:,:),
!!                    real(IN)    :: areaLeft(:,:,:),
!!                    real(IN)    :: secondCoord(MAXCELLS or numCells),
!!                    real(IN)    :: thirdCoord(MAXCELLS or numCells))
!!
!!
!! DESCRIPTION
!!
!!
!! Diffuse_visc alters the velocity fluxes for a block to account for viscosity.
!!
!! Energy fluxes (for E_FLUX and, if defined, EINT_FLUX) are also adjusted
!! accordingly.
!!
!! This is an explicit method, so a timestep limiter will be required.  
!! Stability is guaranteed for 
!! 
!!                           dx**2 
!!                 dt < .5* -------
!!                            nu
!!
!! Fluxes in PPM are stored at the zone boundaries, the temperature, etc. are
!! at the zone centers.  So the velocity fluxes (more exactly, momentum fluxes)
!! acquire terms like  F = -nu * grad(v), which is approximated as
!! 
!!                       nu_i + nu_{i-1}      v(xyz)_i - v(xyz)_{i-1}
!!         F(xyz)_i = - -----------------  * ------------------------   .
!!                              2                      dx
!!
!!
!! This routine computes the velocity and energy fluxes from
!! viscosity and adds them to the corresponding fluxes returned from
!! hydro_1d.  It is to be called after calling hydro_1d on a block.
!! The updated energy fluxes are then used in hy_ppm_updateSoln to
!! produce the updated velocities and energy.
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
!!  temp_flx    -   Temporary storage for fluxes along sweep direction
!!  areaLeft    -   Cell face areas at the left (smaller) cell side
!!
!!  secondCoord -  for an x sweep: y coordinate; for a y sweep: x coord; for a z sweep: x coord
!!  thirdCoord  -  for an x sweep: z coordinate; for a y sweep: z coord; for a z sweep: y coord
!!  blkLimits -   endpoints of block indices without including gcells
!!  blkLimitsGC -   endpoints of block indices including gcells
!!
!!
!! NOTES
!!
!!    This routine is used by the PPM implementation of the Hydro unit.
!!    Other implementations of Hydro, in particular MHD and RHD, may have their
!!    own mechanisms for handling diffusive effects that bypass Diffuse_visc
!!    (and the Diffuse code unit in general), or they may lack support for
!!    diffusive effects.
!!
!! SEE ALSO
!!
!!  Diffuse_therm
!!
!! HISTORY
!!
!!  Apparently this subroutine used to be called visc_explicit.
!!  Visc_explicit started out as a very lightly modified version of therm_explicit.
!!
!!***

subroutine Diffuse_visc(sweepDir, igeom, blockID,numCells,blkLimits,blkLimitsGC,&
                           leftCoords,rightCoords,&
                           temp_flx, areaLeft, secondCoord, thirdCoord)


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
  real, intent(IN), DIMENSION(MAXCELLS) :: leftCoords ,rightCoords, &
       secondCoord, thirdCoord
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
  real, intent(IN), DIMENSION(numCells) :: leftCoords ,rightCoords, &
       secondCoord, thirdCoord
#endif

  return
end subroutine Diffuse_visc
