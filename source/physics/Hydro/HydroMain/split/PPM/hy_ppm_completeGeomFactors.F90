!!****if* source/physics/Hydro/HydroMain/split/PPM/hy_ppm_completeGeomFactors
!!
!! NAME
!!
!!  hy_ppm_completeGeomFactors
!!
!!
!! SYNOPSIS
!!
!!  hy_ppm_completeGeomFactors(integer(IN) :: numIntCells4,
!!       integer(IN) :: numCells, 
!!       integer(IN) :: igeom, 
!!       real(IN)    :: dx(numCells), 
!!       real(IN)    :: r,
!!       real(OUT)   :: dvol(numCells), 
!!       real(IN)    :: cvol(numCells), 
!!       real(OUT)   :: area(numCells), 
!!       real(IN)    :: areaLeft(numCells))
!!
!! DESCRIPTION
!!
!!  Complete the geometrical factors needed for gradient and 
!!  divergence terms in the hydro equations.  Given xl, xr, and dx, compute
!!  the areas at the left and right faces of the zone (areal and arear), and 
!!  the average area of the zone face, (area) and the volume of the zone.
!!
!!  
!! NOTES
!!
!!  Using Grid_getBlkData to get geometry factors and then completing them with
!!  this routine is an alternative to calling geom.
!!
!!  Geom computes cell areas and volumes in a hybrid form: only some factors
!!  are included. For example, areas are always 1 in Cartestian geometry.
!!  The alternative way using Grid_getBlkData and this subroutine, on the
!!  other hand, represents areas in real physical units: areas are consistent
!!  with the units to express the physical extent of the domain and of cell
!!  lengths.
!!
!! ARGUMENTS
!!
!!  numIntCells4 - The number of interior cells in the current sweep direction plus 4
!!
!!  numCells --    The maximum number of cells in any sweep direction
!!
!!  igeom --       the geometry of the current direction
!!
!!  dx --          The cell lenght in the sweep direction. Note that angle coordinates
!!                 are given as such, not as arc lengths.
!!
!!  r  --          cell-centered radial coordinate for the cells in current sweep
!!
!!  dvol --        cell volumes, to be used by caller
!!
!!  cvol --        cell volumes from calling Grid_getBlkData(...,CELL_VOLUME,...)
!!
!!  area --        average cross-scectional areas normal too sweep direction, to be used by caller
!!
!!  areaLeft --    left cell face areas from calling Grid_getBlkData(...,CELL_FACEAREA,...)
!!
!! SEE ALSO
!!
!!  geom
!!  Grid_getBlkData
!!
!!***                


subroutine hy_ppm_completeGeomFactors (numIntCells4, numCells, igeom,   &
                 dx, r, dvol, cvol,        &
                 area, areaLeft)

  
  implicit none
  
#include "constants.h"

  integer, INTENT(in)                   :: igeom, numIntCells4, numCells
  real, DIMENSION(numCells),INTENT(in)  :: dx
  real, INTENT(in)                      :: r !        the radial coordinate.
  real, DIMENSION(numCells),intent(out) :: area, dvol
  real, DIMENSION(numCells),intent(in)  :: areaLeft, cvol

  integer :: i




!!$            Maybe we ought to test here that the user has supplied a consistent
!!$            set of coordinates.  For now we assume that the coordinates
!!$            are consistent. There is some testing in geom that could be repeated here.


  dvol(5:numIntCells4)  = cvol(5:numIntCells4)



  select case (igeom)
     
  case (XYZ,PHI_CYL,PHI_SPH)    ! Cartesian or planar geometry
                                ! or Cylindrical geometry, angular coord.
                                ! or Spherical geometry, phi coord.
     area(5:numIntCells4) = areaLeft(5:numIntCells4)
     
  case (RAD_CYL)            ! Cylindrical geometry, radial coord.
     do i = 5,numIntCells4
        area(i) = 0.5 * (areaLeft(i)+areaLeft(i+1))
     enddo
     
  case (RAD_SPH)            ! Spherical geometry, radial coord.
     do i = 5, numIntCells4
        area(i)  = dvol(i) / dx(i)
     enddo
     
  case (THETA)            ! Spherical geometry, theta coord.
     do i = 5, numIntCells4
        area(i)  = dvol(i) / ( dx(i) * r )
     enddo
     
        
  end select


end subroutine hy_ppm_completeGeomFactors
