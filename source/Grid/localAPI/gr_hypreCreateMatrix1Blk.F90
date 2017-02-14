!!****if* source/Grid/localAPI/gr_hypreCreateMatrix1Blk
!!
!!  NAME 
!!
!!   gr_hypreCreateMatrix1Blk
!!
!!  SYNOPSIS
!!
!!   call gr_hypreCreateMatrix1Blk(integer(IN)           :: iVar,
!!                             integer(IN)           :: iFactorB,
!!                             integer(IN)           :: iFactorA,
!!                             integer(IN)           :: bcTypes(6),
!!                             real(IN)              :: bcValues(2,6),
!!                             real(IN)              :: dt,
!!                             real(IN)              :: alpha,
!!                             integer(IN)           :: blockCount,
!!                             integer(IN)           :: blockList(blockCount),
!!                             logical(IN)           :: JacobiMatrix)
!!
!!
!!  DESCRIPTION
!!      This routine is used by some implementations of  gr_hypreCreateMatrix
!!      do do most of the per-block work.
!!
!!
!! ARGUMENTS
!! 
!!   iVar         : Variable on which the diffusion operation is performed (e.g., TEMP_VAR)
!!   iFactorB     : a factors in the equation with spatial variation.
!!   iFactorA     : a factors in the equation with spatial variation.
!!   bcTypes      : Boundary condition types.  Should be chosen from the constants
!!                  GRID_PDE_BND_* defined in Grid_interface.F90, or the special value VACUUM
!!                  defined in constants.h.
!!                  Presently OUTFLOW and VACUUM are supported, DIRICHLET less well tested.
!!   bcValues     : Values of iVar,iFactorB (!DEV: ??) on boundary (currently used for DIRICHLET and GIVENGRAD BCs).                        
!!   dt           : The time step.
!!   alpha        : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank-Nicolson
!!   blockCount   : The number of blocks in the list.   
!!   blockList    : The list of blocks on which the solution must be updated.        
!!   JacobiMatrix : TRUE, computes A and FALSE computes M.
!!
!! NOTES
!!
!!  This is a stub implementation.
!!
!! SIDE EFFECTS
!!
!!   On return, various HYPRE objects have been modified.
!!
!! SEE ALSO
!!
!!  gr_hypreCreateMatrix
!!***

subroutine gr_hypreCreateMatrix1Blk(iFactorB, dt, &
     theta, &
     blkLimits, datasize, solnVec, nentries, faces, mylevel, var, faceAreas, &
     del, negh_del, xflux,yflux,zflux, &
     lb, blockID)
  
  implicit none
 
#include "constants.h"
  
  integer, intent(IN) :: iFactorB
  real,    intent(IN) :: dt
  real,    intent(IN) :: theta
  real, intent(IN), dimension(MDIM)     :: del
  real, intent(IN), dimension(2*MDIM, MDIM) :: negh_del
  integer, intent(IN), dimension(2,MDIM):: blkLimits 
  integer, intent(IN) :: datasize(MDIM)
  integer, intent(IN) ::  mylevel
  integer, intent(IN) ::  var
  integer, intent(IN) ::  nentries
  integer, dimension(2,MDIM), intent(IN) :: faces 
  real, intent(IN), dimension(:,:,:,:) :: xflux, yflux, zflux
  real, intent(IN) :: faceAreas  (:,:,:,:)  
  integer, intent(IN) :: lb, blockID
      
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec

end subroutine gr_hypreCreateMatrix1Blk
