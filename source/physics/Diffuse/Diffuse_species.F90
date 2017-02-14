!!****f* source/physics/Diffuse/Diffuse_species
!!
!! NAME
!!
!!  Diffuse_species
!!
!!
!! SYNOPSIS
!!
!! Diffuse_species(integer(IN) :sweepDir, 
!!                 integer(IN) : igeom, 
!!                 integer(IN) : blockID,
!!                 integer(IN) : numCells,
!!                 integer(:,:)(IN):blkLimits, 
!!                 integer(:,:)(IN):blkLimitsGC
!!                 real(:)(IN):leftCoords,
!!                 real(:)(IN):rightCoords,
!!                 real(:,:)(INOUT) :temp_flx, 
!!                 real(:,:)(INOUT) :temp_fly, 
!!                 real(:,:)(INOUT) :temp_flz)
!!
!!
!! DESCRIPTION
!!
!!  Diffuse_species alters the mass-fraction flux to account for 
!!  material diffusion.  This is an explicit method, so a timestep limiter
!!  will be required.  Stability is guaranteed for
!!
!!                           dx**2
!!                 dt < .5* -------
!!                             D
!!
!!  where D is the diffusion coefficient.
!!
!!  Fluxes in PPM are stored at the zone boundaries; mass fractions are 
!!  at the zone centers.  So the flux {F = -sigma * grad(T)} is
!!
!!                      D_i + D_{i-1}     T_i - T_{i-1}
!!             F_i = - -------------- * ---------------
!!                           2                  dx
!!
!! NOTES
!!
!!    This routine computes the species fluxes and adds them to the flux
!!    returned from hydro.  It is to be called after calling hydro_1d on a
!!    block.  The updated energy fluxes are then used in update_soln to
!!    produce the updated energy.
!!
!!    Boundary conditions are handled through tot_bnd,
!!    not in this routine.  Ex: setting a reflecting boundary will reflect the
!!    species into the guard cells, which are used in computing the flux
!!    at the block boundaries.
!!
!!
!! ARGUMENTS
!!
!!  sweepDir        the current sweep direction
!!
!!  igeom           the geometry flag for the current sweep direction
!!
!!  blockID        the block number to operate on
!!  numCells        the number of Cells along the sweep direction
!!  leftCoords      Coordinates of the left edge of the zones
!!  rightCoords     Coordinates of the right edge of the zones
!!
!!  tempflx         Temporary storage for flux along first direction
!!  tempfly         Temporary storage for flux along second direction
!!  tempflz         Temporary storage for flux along third direction
!!
!!***

subroutine Diffuse_species(sweepDir, igeom, blockID,numCells,blkLimits,blkLimitsGC,&
                           leftCoords,rightCoords,&
                           temp_flx, temp_fly, temp_flz)


    
  
  implicit none
#include "constants.h"  
#include "Flash.h"
  
  integer, intent(IN) :: sweepDir, igeom, blockID, numCells
  integer, dimension(2,MDIM), intent (IN) :: blkLimitsGC, blkLimits
  
  
  real, intent(INOUT), DIMENSION(NFLUXES,                   &
                               blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
                               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
                               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  &
                               temp_flx, temp_fly, temp_flz
  
  real, intent(IN), DIMENSION(numCells) :: leftCoords ,rightCoords


  
  return
end subroutine Diffuse_species
