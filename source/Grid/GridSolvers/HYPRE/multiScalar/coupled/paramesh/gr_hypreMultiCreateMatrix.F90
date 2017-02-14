!!****if* source/Grid/GridSolvers/HYPRE/multiScalar/coupled/paramesh/gr_hypreMultiCreateMatrix
!!
!!  NAME 
!!
!!   gr_hypreMultiCreateMatrix
!!
!!  SYNOPSIS
!!
!!   call gr_hypreMultiCreateMatrix(integer(IN)           :: iVar,
!!                             integer(IN)           :: iFactorB,
!!                         OPTIONAL,integer(IN)           :: iFactorA,
!!                             integer(IN)           :: bcTypes(6),
!!                             real(IN)              :: bcValues(2,6),
!!                             real(IN)              :: dt,
!!                             real(IN)              :: theta,
!!                             integer(IN)           :: blockCount,
!!                             integer(IN)           :: blockList(blockCount),
!!                             logical(IN)           :: JacobiMatrix)
!!
!!
!!  DESCRIPTION 
!!      This routine computes one of the matrices A and B, depending on the
!!      logical input argument 'JacobiMatrix':
!!          Ax = b, where A is the matrix to be inverted
!!          B = MX, where M is a matrix whose product with iVar produces RHS B.
!!
!!      A*(df/dt) + C*f = div(B*grad(f)) + D
!!      f -> Variable to be diffused.
!!      C,D are optional factors (not implemented here, the caller should add them later.)
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
!!   theta        : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank-Nicolson
!!   blockCount   : The number of blocks in the list.   
!!   blockList    : The list of blocks on which the solution must be updated.        
!!   JacobiMatrix : TRUE, computes A and FALSE computes M.
!!
!! SIDE EFFECTS
!!
!!   On return, the elements of the HYPRE matrix A that represent the second-derivative
!!   operator (div B grad) term have been defined, and it is ready for use.
!!   The gr_hypreData module variable gr_hypreMatA holds the
!!   handle for the HYPRE Solver object.
!!
!! NOTES
!!
!!   This routine does not actually 'create' the matrix object in the sense of HYPRE.
!!   It expects a matrix object already created and initialized; that is done,
!!   together with initialization of the grid object, in gr_hypreSetupGrid.
!!
!!   Currently, gr_hypreMultiCreateMatrix is called from Grid_advanceDiffusion with
!!   JacobiMatrix==.FALSE. only when the implicitness parameter theta passed to
!!   Grid_advanceDiffusion is 0. (KW 2012-12-05, corrected 2014-12-05)
!!
!! SEE ALSO
!!
!!  Grid_interface
!!***

subroutine gr_hypreMultiCreateMatrix(iVar, iFactorB, iFactorA, bcTypes, bcValues, dt, &
     theta, blockCount, blockList, JacobiMatrix)
  
  use gr_hypreData,     ONLY : gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, &
                               gr_hypreRefineMIN, &
                               gr_hypreNVars
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkRefineLevel
  use Timers_interface, ONLY : Timers_start, Timers_stop 
  
  implicit none
 
#include "Flash.h" 
#include "constants.h"
#include "HYPREf.h"    
  
  !!-----------------------------------------------------------------------
  !!         ARGUMENTS
  !!-----------------------------------------------------------------------
  integer, intent(IN) :: iVar
  integer, intent(IN) :: iFactorB
  integer, intent(IN) :: iFactorA
  integer, intent(IN) :: bcTypes(6)
  real,    intent(IN) :: bcValues(2,6)
  real,    intent(IN) :: dt
  real,    intent(IN) :: theta
  integer, intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  logical, intent(IN) :: JacobiMatrix
      
  !!-----------------------------------------------------------------------
  !!         LOCAL VARIABLES.
  !!-----------------------------------------------------------------------  
  integer :: ierr
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  integer :: datasize(MDIM)
  integer ::  mypart,mylevel
  integer ::  var
  integer :: blockID
  integer ::  nentries
  integer,allocatable ::  stencilIndices(:)
  integer :: i, lb

  character(len=32) :: matfile
  real, allocatable :: BoxVal(:)
  
  call Timers_start("gr_hypreMultiCreateMatrix")  
  
  
  nentries = gr_hypreNVars
  allocate(stencilIndices(nentries))
  do i = 1, nentries
     stencilIndices(i) = i-1
  enddo

  mypart = 0  !! part iterator 
  var    = 0  !! var iterator.
  
  
  do lb = 1, blockCount     
     blockID = blockList(lb)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     call Grid_getBlkRefineLevel(blockID,mylevel)

     datasize  (1:MDIM)= blkLimits  (HIGH,1:MDIM)-blkLimits  (LOW,1:MDIM)+1
     
     allocate(BoxVal(nentries*product(datasize(1:NDIM))))
     
     mypart = mylevel - gr_hypreRefineMIN                  
     
     BoxVal = 0.0
     
     call HYPRE_SStructMatrixSetBoxValues(gr_hypreMatA, mypart, gr_hypreLower(lb,1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), var, nentries, stencilIndices(1:nentries), BoxVal(:), ierr)


     deallocate (BoxVal)

  end do !! block

 deallocate(stencilIndices)

 call Timers_stop("gr_hypreMultiCreateMatrix")

 return

end subroutine gr_hypreMultiCreateMatrix
