!!****if* source/physics/Hydro/HydroMain/split/PPM/hy_ppm_geom
!!
!! NAME
!!
!!  hy_ppm_geom
!!
!!
!! SYNOPSIS
!!
!!  call hy_ppm_geom(integer(IN) :: numIntCells,
!!               integer(IN) :: numCells, 
!!               integer(IN) :: j, 
!!               integer(IN) :: k, 
!!               integer(IN) :: xyzswp, 
!!               integer(IN) :: igeom, 
!!               real(OUT)   :: areal(numCells), 
!!               real(OUT)   :: arear(numCells), 
!!               real(OUT)   :: area(numCells), 
!!               real(IN)    :: dx(numCells), 
!!               real(OUT)   :: dvol(numCells), 
!!               real(IN)    :: xl(numCells), 
!!               real(IN)    :: xr(numCells),
!!               real(IN)    :: radialCoord(numCells),
!!               real(IN)    :: thirdCoord(numCells))
!!
!! DESCRIPTION
!!
!!  Compute the geometrical factors needed for the gradient and 
!!  divergence terms in the hydro equations.  Given xl, xr, and dx, compute
!!  the areas at the left and right faces of the zone (areal and arear), and 
!!  the average area of the zone face, (area) and the volume of the zone.
!!
!!  
!! NOTES
!!
!!  This routine only works for 2-d cylindrical (r,z) and 1, 2, and 3-d 
!!  Cartesian geometries.   DEV: Verify this claim!
!!
!!
!! ARGUMENTS
!!
!!  numIntCells -- The number of zones in the current sweep direction
!!
!!  numCells --    The maximum number of zones in any sweep direction
!!
!!  j, k --        indices of the zones in the two coordinate directions that
!!                 are not the current sweep direction (in coordinate order)
!!
!!  xyzswp --      the current sweep direction for the directionally split scheme
!!
!!  igeom --       the geometry of the current direction
!!
!!  areal, arear-- The area at the left and right boundaries
!! 
!!  area --        The average cross-sectional area normal to the sweep direction
!!
!!  dx --          The zone width in the sweep direction
!!
!!  dvol --        The differential volume of the zone used in a divergence
!!
!!  xl, xr --      The coordinates at the left and right edges of the current zone
!!
!!  radialCoord -  cell-centered radial coordinates for the cells in current block
!!
!!  thirdCoord -   cell-centered values of third coordinate (second transversal
!!                 direction to current sweep) for cells in the current block
!!
!!
!!***                


subroutine hy_ppm_geom (numIntCells, numCells,j, k, xyzswp, igeom,   &
                 areal, arear, area, dx, dvol, xl, xr,        &
                 radialCoord, thirdCoord)
  use Driver_interface, ONLY : Driver_abortFlash

  
  implicit none
  
#include "constants.h"
#include "Flash.h"

  integer, INTENT(in)              :: j, k, xyzswp, igeom, numIntCells, numCells
  real, DIMENSION(numCells),INTENT(in)    :: dx,xl, xr
  real, DIMENSION(numCells),INTENT(in) :: radialCoord, thirdCoord
  real, DIMENSION(numCells),INTENT(out) :: areal, arear, area, dvol

  integer :: i

  integer :: numIntCells_max

!!$            the radial coordinate.

  if ( (igeom .eq. RAD_CYL .or. igeom .eq. RAD_SPH) .and. &
       &  (xyzswp .ne. SWEEP_X) ) then
     write (*,*) 'fatal in hy_ppm_geom:  trying to do a radial sweep in y or z'
     call Driver_abortFlash("Fatal error in hy_ppm_geom: trying to do a radial sweep in y or z")
     stop
  endif


!!$            We ought to test here that the user has supplied a consistent
!!$            set of coordinates.  For now we assume that the coordinates
!!$            are consistent.
!!$
!!$-------------------------------------------------------------------------------






  numIntCells_max = numIntCells + 2*NGUARD

  select case (igeom)
     
  case (XYZ)        ! Cartesian or planar geometry
     
     do i = 1, numIntCells_max
        areal(i) = 1.e00
        arear(i) = 1.e00
        area(i)  = 1.e00
        dvol(i)  = dx(i)
     enddo
     
  case (RAD_CYL)            ! Cylindrical geometry, radial coord.
     
     do i = 1, numIntCells_max
        areal(i) = abs (xl(i))
        arear(i) = abs (xr(i))
        area(i)  = 0.5e00 * (arear(i) + areal(i))
        dvol(i)  = area(i) * dx(i)
     enddo
     
  case (RAD_SPH)            ! Spherical geometry, radial coord.
     
     do i = 1, numIntCells_max
        areal(i) = xl(i) * xl(i)
        arear(i) = xr(i) * xr(i)
        dvol(i)  = (xr(i) * arear(i) - xl(i) * areal(i)) / 3.e00
        area(i)  = dvol(i) / dx(i)
     enddo
     
  case (PHI_CYL)            ! Cylindrical geometry, angular coord.
     
     do i = 1, numIntCells_max
        areal(i) = 1.e00
        arear(i) = 1.e00
        area(i)  = 1.e00
        dvol(i)  = area(i) * dx(i) * radialCoord(j)
     enddo
     
  case (THETA)            ! Spherical geometry, theta coord.
     
     do i = 1, numIntCells_max
        areal(i) = sin (xl(i))
        arear(i) = sin (xr(i))
        area(i)  = sin (0.5e00 * (xl(i) + xr(i)))
        dvol(i)  = area(i) * dx(i) * radialCoord(j)
     enddo
     
  case (PHI_SPH)            ! Spherical geometry, phi coord.

!!$ Whether it's the SWEEP_X or the SWEEP_Y sweep, 
!!$ the theta coordinate (indexed by k) will be in thirdCoord.
     
     do i = 1, numIntCells_max
        areal(i) = 1.e00
        arear(i) = 1.e00
        area(i)  = 1.e00
        dvol(i)  = dx(i) * radialCoord(j) * sin(thirdCoord(k))
     enddo
        
  end select


  return
end subroutine hy_ppm_geom
