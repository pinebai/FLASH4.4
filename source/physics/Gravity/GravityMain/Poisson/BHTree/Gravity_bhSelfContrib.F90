!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhSelfContrib
!!
!! NAME
!!
!!  Gravity_bhSelfContrib
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhSelfContrib(
!!                  real(in)       :: node(:),
!!                  real(in)       :: cellsize(MDIM),
!!                  integer(in)    :: blockno,
!!                  integer(in)    :: point(MDIM),
!!                  integer(in)    :: blkLimits(2,MDIM),
!!                  real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Calculates contribution of the cell at the point-of-calculation to 
!!  the gravitational potential and adds it to the solnData array (with 
!!  index grv_defaultGpotVar).
!!
!! ARGUMENTS
!!
!!  node        : array of the node of the tree, whose
!!                contribution is added to the gravitational potential
!!  cellsize    : physical size of the cell (in each dimension)
!!  blockno     : number of block into which the point-of-calculation belongs
!!  point       : indeces of the point-of-calculation in the block
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!!***

subroutine Gravity_bhSelfContrib(node, cellsize, blockno, point, blkLimits, solnData)
  use Gravity_data, ONLY : useGravity, grv_bhPiGhalf, grv_bhSContrGeom, grv_bhIM, &
    grv_defaultGpotVar, grv_bhA1Dir
  implicit none
#include "constants.h"
#include "FortranLangFeatures.fh"
  real, dimension(:), intent(IN) :: node
  real, dimension(MDIM), intent(IN) :: cellsize
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData

  if (.not. useGravity) return

  ! add the contribution to the potential
  !print *, "SC start: ", cellsize, blockno, point
  !print *, "SC1: ",  grav_twopiG, solnData(grv_bhDensVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
  !& , grav_scontr_geom
  !print *, "SC2: ", cellsize, grv_bhTreeA2Dir, grv_bhTreeA3Dir
  !print *, "SC3: ", cellsize(grv_bhTreeA2Dir), cellsize(grv_bhTreeA3Dir)

  ! grav. pot. in the centre of a tri-axial ellipsoid (Galactic Dynamics):
  ! Phi = - 2*pi*G*a2*a3*F(theta,k)/sin(theta)
  ! where F in incomplete elliptic integral
  ! theta = acos(a3/a1) and k = sqrt((a1^2 - a2^2)/(a1^2-a3^2))
  ! and a1, a2 and a3 are the ellipsoid semi-axes (a1>a2>a3)
  ! a[2,3] = 0.5*cellsize(grv_bhTreeA[2,3]Dir) (therefore 2*pi*G -> 0.5*pi*G)
  solnData(grv_defaultGpotVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
  & = solnData(grv_defaultGpotVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
  & - grv_bhPiGhalf * grv_bhSContrGeom * node(grv_bhIM) &
  & / cellsize(grv_bhA1Dir)

  ! alternative, using density in solnData instead of the node mass:
  !solnData(grv_defaultGpotVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
  !& = solnData(grv_defaultGpotVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
  !& - grav_piGhalf * solnData(grv_bhDensVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
  !& * grav_scontr_geom * cellsize(grv_bhTreeA2Dir) * cellsize(grv_bhTreeA3Dir)


  return
end subroutine Gravity_bhSelfContrib

