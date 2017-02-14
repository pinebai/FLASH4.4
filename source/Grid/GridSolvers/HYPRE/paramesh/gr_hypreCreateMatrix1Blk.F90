!!****if* source/Grid/GridSolvers/HYPRE/paramesh/gr_hypreCreateMatrix1Blk
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
!!   alpha        : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank-Nicolson
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
!!   Currently, gr_hypreCreateMatrixFcB is called from Grid_advanceDiffusion with
!!   JacobiMatrix==.FALSE. only when the implicitness parameter theta passed to
!!   Grid_advanceDiffusion is 0. (KW 2012-12-05, corrected 2014-12-05)
!!
!! SEE ALSO
!!
!!  Grid_interface
!!***

!!REORDER(4): solnVec

subroutine gr_hypreCreateMatrix1Blk(iFactorB, dt, &
     theta, &
     blkLimits, datasize, solnVec, nentries, faces, mylevel, var, faceAreas, &
     del, negh_del, xflux,yflux,zflux, &
     lb, blockID)
  
  use gr_hypreData,     ONLY : gr_hypreSolverType,gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, &
                               gr_hypreRefineMIN, gr_hypreRefineMAX, gr_hypreNeghLevels, &
                               gr_hypreSurrBlkSum
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_fillGuardCells, Grid_getBlkBC, &
    Grid_getBlkCornerID, Grid_getCellCoords, Grid_getFluxData, Grid_getBlkData, Grid_getDeltas, Grid_getBlkRefineLevel
  use Timers_interface, ONLY : Timers_start, Timers_stop 
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
                               GRID_PDE_BND_NEUMANN,   &
                               GRID_PDE_BND_DIRICHLET
  
  
  implicit none
 
#include "Flash.h" 
#include "constants.h"
#include "Flash_mpi.h"
#include "HYPREf.h"    
  
  !!-----------------------------------------------------------------------
  !!         ARGUMENTS
  !!-----------------------------------------------------------------------
  integer, intent(IN) :: iFactorB
  real,    intent(IN) :: dt
  real,    intent(IN) :: theta
  real, intent(IN), dimension(MDIM)     :: del
  real, intent(IN), dimension(2*MDIM, MDIM) :: negh_del
!  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  integer, intent(IN), dimension(2,MDIM):: blkLimits 
!  integer :: datasize(MDIM), datasizeGC(MDIM)
  integer, intent(IN) :: datasize(MDIM)
  integer, intent(IN) ::  mylevel
  integer, intent(IN) ::  var
  integer, intent(IN) ::  nentries
  integer, dimension(2,MDIM), intent(IN) :: faces 
  real, intent(IN), dimension(:,:,:,:) :: xflux, yflux, zflux
  real, intent(IN) :: faceAreas  (:,:,:,:)  
  integer, intent(IN) :: lb, blockID
      
  !!-----------------------------------------------------------------------
  !!         LOCAL VARIABLES.
  !!-----------------------------------------------------------------------  
  integer :: ierr, pos(NDIM)
  real, dimension(MDIM)     :: delph, delmh
  real :: CDiv 
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer ::  mypart
  integer ::  stencil_indices(7)  
  real    ::  values(7), graph_value(7) 
  integer :: i, j, k, ii
  real   :: condiph, condimh
  real   :: condjph, condjmh
  real   :: condkph, condkmh
  integer :: eachNegh,numNegh 

  integer :: dir
  integer :: numGraph, iter
  real, allocatable :: BoxVal(:)
  integer, parameter :: nFluxVars = 2**(NDIM-1)
  
  do i = 1, nentries
     stencil_indices(i) = i-1
  enddo
  
     
     allocate(BoxVal(nentries*product(datasize(1:NDIM))))
     
     mypart = mylevel - gr_hypreRefineMIN                  
     
     
     iter = 1
     
     
     do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)                                   
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)                    
              
              !!-----------------------------------------------------------------------
              !!         SWAP INDICES FOR FLASH-HYPRE COMPATIBILITY.
              !!-----------------------------------------------------------------------
#if NDIM == 1
              pos(1) = gr_hypreLower(lb,1) + (i-blkLimits(LOW,IAXIS))
#endif
#if NDIM == 2                                 
              pos(1) = gr_hypreLower(lb,1) + (j-blkLimits(LOW,JAXIS))
              pos(2) = gr_hypreLower(lb,2) + (i-blkLimits(LOW,IAXIS))
#endif
#if NDIM ==3
              pos(1) = gr_hypreLower(lb,1) + (k-blkLimits(LOW,KAXIS))
              pos(2) = gr_hypreLower(lb,2) + (j-blkLimits(LOW,JAXIS))
              pos(3) = gr_hypreLower(lb,3) + (i-blkLimits(LOW,IAXIS))
#endif         
              
              !!-----------------------------------------------------------------------
              !!         STENCILED VALUES TO BE FED INTO HYPRE MATRIX.
              !!-----------------------------------------------------------------------
              values = 0.0                           
              
              delmh = del
              delph = del
              
              condimh = 0.
              condiph = 0.
              
              numGraph = 0
              
              !! i-1,j,k
              BoxVal(iter+1) = 0.0
              if ((i /= blkLimits(LOW, IAXIS)) .or. (faces(1,IAXIS) == NOT_BOUNDARY)) then                                         
                 if (i == blkLimits(LOW, IAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,LEFT_EDGE, 1+K2D, 1+K3D)) then  !! F/C boundary.                                         
                    numNegh = gr_hypreSurrBlkSum(lb) % regionInfo(LEFT_EDGE,1+K2D,1+K3D) % numNegh    
                    numGraph = numNegh                    
                    do eachNegh = 1, numNegh                                                            
                       ii = 2*NDIM+1+eachNegh-1                       
                       delmh(IAXIS) =  0.5*(negh_del(1,IAXIS) + del(IAXIS))
                       
                       if (numNegh > 1) then !! we are on coarse cell, use the average computed from fine cell (2 of them)
                          condimh = xflux(eachNegh,i,j,k)
                       else
                          if (NDIM > 1) then
                             !! we are on fine cell and looking at coarse block, use regular averaging.
                             condimh = 0.5*(solnVec(iFactorB,i-1,j,k)+ solnVec(iFactorB,i,j,k))*faceAreas(i,j,k,IAXIS)
                          else
                             condimh = xflux(eachNegh,i,j,k)
                          end if
                       end if
                       
                       graph_value(1) = -condimh*theta*dt/delmh(IAXIS)
                       call HYPRE_SStructMatrixSetValues(gr_hypreMatA, mypart, pos(1:NDIM), var, 1, ii, graph_value,ierr)
                    end do

                    if (numNegh > 1) then
                       condimh = sum(xflux(1:nFluxVars,i,j,k))
                    end if
                 else
                    !! STENCILED RELATIONSHIP.
                    delmh(IAXIS) = del(IAXIS)
                    condimh = 0.5*(solnVec(iFactorB,i-1,j,k)+ solnVec(iFactorB,i,j,k))*faceAreas(i,j,k,IAXIS)
                    BoxVal(iter+1) = -condimh*theta*dt/delmh(IAXIS)
                 end if
              end if

              !! i+1,j,k
              BoxVal(iter+2) = 0.0
              if (i /= blkLimits(HIGH, IAXIS) .or. (faces(2,IAXIS) == NOT_BOUNDARY)) then
                 if (i == blkLimits(HIGH, IAXIS) .and. mylevel /= gr_hypreNeghLevels(lb, RIGHT_EDGE, 1+K2D, 1+K3D)) then !! F/C boundary.
                    numNegh = gr_hypreSurrBlkSum(lb) % regionInfo(RIGHT_EDGE, 1+K2D, 1+K3D) % numNegh
                    numGraph = numNegh
                    do eachNegh = 1, numNegh
                       ii= 2*NDIM+1+eachNegh-1
                       delph(IAXIS) =  0.5*(negh_del(2,IAXIS) + del(IAXIS))
                       if (numNegh > 1) then !! we are on coarse cell, use fine cell area on fluxes.
                          condiph = xflux(eachNegh,i+1,j,k)
                       else
                          if (NDIM > 1) then
                             condiph = 0.5*(solnVec(iFactorB,i,j,k)  + solnVec(iFactorB,i+1,j,k))*faceAreas(i+1,j,k,IAXIS)
                          else
                             condiph = xflux(eachNegh,i+1,j,k)
                          end if
                       end if
                       graph_value(1) =  -condiph*theta*dt/delph(IAXIS)
                       call HYPRE_SStructMatrixSetValues(gr_hypreMatA, mypart, pos(1:NDIM), var, 1, ii, graph_value,ierr)
                    end do

                    if (numNegh > 1) then
                       condiph = sum(xflux(1:nFluxVars,i+1,j,k)) !! actual flux, goes towards computing diag.
                    end if
                 else
                    !! STENCILED RELATIONSHIP.
                    delph(IAXIS) =  del(IAXIS)
                    condiph = 0.5*(solnVec(iFactorB,i,j,k)  + solnVec(iFactorB,i+1,j,k))*faceAreas(i+1,j,k,IAXIS)
                    BoxVal(iter+2) =  -condiph*theta*dt/(delph(IAXIS))
                 end if
              end if

              BoxVal(iter) = theta*((Condimh/delmh(1))+(Condiph/delph(1)))*dt

#if NDIM >= 2
              condjmh = 0.
              condjph = 0.

              !! i,j-1,k
              BoxVal(iter+3) = 0.0
              if ((j /= blkLimits(LOW, JAXIS)) .or. (faces(1,JAXIS) == NOT_BOUNDARY)) then
                 if (j == blkLimits(LOW, JAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,1+K1D,LEFT_EDGE,1+K3D)) then  !! F/C boundary.
                    numNegh = gr_hypreSurrBlkSum(lb) % regionInfo(1+K1D,LEFT_EDGE,1+K3D) % numNegh
                    do eachNegh = 1, numNegh
                       delmh(JAXIS) =  0.5*(negh_del(3,JAXIS) + del(JAXIS))
                       ii=2*NDIM+1+eachNegh-1+numGraph
                       if (numNegh > 1) then !! we are on coarse cell, use fine cell area on fluxes.
                          condjmh = yflux(eachNegh,i,j,k)
                       else
                          condjmh = 0.5*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*faceAreas(i,j,k,JAXIS)
                       end if
                       graph_value(1) = -condjmh*theta*dt/delmh(JAXIS)
                       call HYPRE_SStructMatrixSetValues(gr_hypreMatA, mypart, pos(1:NDIM), var, 1, ii, graph_value,ierr)
                    end do

                    numGraph = numGraph + numNegh

                    if (numNegh > 1) then
                       condjmh = sum(yflux(1:nFluxVars,i,j,k))
                    end if

                 else
                    !! STENCILED RELATIONSHIP.
                    delmh(JAXIS) = del(JAXIS)
                    condjmh = 0.5*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*faceAreas(i,j,k,JAXIS)
                    BoxVal(iter+3) = -condjmh*theta*dt/(delmh(JAXIS))
                 end if
              end if

              !! i,j+1,k
              BoxVal(iter+4) = 0.0
              if ((j /= blkLimits(HIGH, JAXIS)) .or. (faces(2,JAXIS) == NOT_BOUNDARY)) then
                 if (j ==  blkLimits(HIGH, JAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,1+K1D,RIGHT_EDGE,1+K3D)) then !! F/C boundary.
                    numNegh = gr_hypreSurrBlkSum(lb) % regionInfo(1+K1D,RIGHT_EDGE,1+K3D) % numNegh
                    do eachNegh = 1, numNegh
                       delph(JAXIS) =  0.5*(negh_del(4,JAXIS) + del(JAXIS))
                       ii=2*NDIM+1+eachNegh-1+numGraph
                       if (numNegh > 1) then !! we are on coarse cell, use fine cell area on fluxes.
                          condjph = yflux(eachNegh,i,j+1,k)
                       else
                          condjph = 0.5*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*faceAreas(i,j+1,k,JAXIS)
                       end if
                       graph_value(1) = -condjph*theta*dt/delph(JAXIS)
                       call HYPRE_SStructMatrixSetValues(gr_hypreMatA, mypart, pos(1:NDIM), var, 1, ii, graph_value,ierr)
                    end do

                    numGraph = numGraph + numNegh

                    if (numNegh > 1) then
                       condjph = sum(yflux(1:nFluxVars,i,j+1,k))
                    end if
                 else
                    !! STENCILED RELATIONSHIP.
                    delph(JAXIS) = del(JAXIS)
                    condjph = 0.5*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*faceAreas(i,j+1,k,JAXIS)
                    BoxVal(iter+4) = -condjph*theta*dt/(delph(JAXIS))
                 end if
              end if

              BoxVal(iter) = BoxVal(iter) +  theta*(Condjmh/delmh(2)+Condjph/delph(2))*dt !! diag

#if NDIM == 3
              condkmh = 0.
              condkph = 0.
              !! i,j,k-1
              BoxVal(iter+5) = 0.0
              if ((k /= blkLimits(LOW, KAXIS)) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) then
                 if (k == blkLimits(LOW, KAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,1+K1D,1+K2D,LEFT_EDGE)) then  !! F/C boundary.
                    numNegh = gr_hypreSurrBlkSum(lb) % regionInfo(1+K1D,1+K2D,LEFT_EDGE) % numNegh
                    do eachNegh = 1, numNegh
                       delmh(KAXIS) =  0.5*(negh_del(5,KAXIS) + del(KAXIS))
                       ii=2*NDIM+1+eachNegh-1+numGraph
                       if (numNegh > 1) then !! we are on coarse cell, use fine cell area on fluxes.
                          condkmh = zflux(eachNegh,i,j,k)
                       else
                          condkmh = 0.5*(solnVec(iFactorB,i,j,k-1)+ solnVec(iFactorB,i,j,k))*faceAreas(i,j,k,KAXIS)
                       end if
                       graph_value(1) = -condkmh*theta*dt/delmh(KAXIS)
                       call HYPRE_SStructMatrixSetValues(gr_hypreMatA, mypart, pos(1:NDIM), var, 1, ii, graph_value,ierr)
                    end do

                    numGraph = numGraph + numNegh

                    if (numNegh > 1) then
                       condkmh = sum(zflux(1:nFluxVars,i,j,k))
                    end if

                 else
                    !! STENCILED RELATIONSHIP.
                    delmh(KAXIS) = del(KAXIS)
                    condkmh = 0.5*(solnVec(iFactorB,i,j,k-1)+ solnVec(iFactorB,i,j,k))*faceAreas(i,j,k,KAXIS)
                    BoxVal(iter+5) = -condkmh*theta*dt/(delmh(KAXIS))
                 end if
              end if

              !! i,j,k+1
              BoxVal(iter+6) = 0.0
              if ((k /= blkLimits(HIGH,KAXIS)) .or. (faces(2,KAXIS) == NOT_BOUNDARY)) then
                 if (k ==  blkLimits(HIGH,KAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,1+K1D,1+K2D,RIGHT_EDGE)) then !! F/C boundary.
                    numNegh = gr_hypreSurrBlkSum(lb) % regionInfo(1+K1D,1+K2D,RIGHT_EDGE) % numNegh
                    do eachNegh = 1, numNegh
                       delph(KAXIS) =  0.5*(negh_del(6,KAXIS) + del(KAXIS))
                       ii=2*NDIM+1+eachNegh-1+numGraph
                       if (numNegh > 1) then !! we are on coarse cell, use fine cell area on fluxes.
                          condkph = zflux(eachNegh,i,j,k+1)
                       else
                          condkph = 0.5*(solnVec(iFactorB,i,j,k+1)+ solnVec(iFactorB,i,j,k))*faceAreas(i,j,k+1,KAXIS)
                       end if
                       graph_value(1) = -condkph*theta*dt/delph(KAXIS)
                       call HYPRE_SStructMatrixSetValues(gr_hypreMatA, mypart, pos(1:NDIM), var, 1, ii, graph_value,ierr)
                    end do

                    numGraph = numGraph + numNegh

                    if (numNegh > 1) then
                       condkph = sum(zflux(1:nFluxVars,i,j,k+1))
                    end if
                 else
                    !! STENCILED RELATIONSHIP.
                    delph(KAXIS) = del(KAXIS)
                    condkph = 0.5*(solnVec(iFactorB,i,j,k+1)+ solnVec(iFactorB,i,j,k))*faceAreas(i,j,k+1,KAXIS)
                    BoxVal(iter+6) = -condkph*theta*dt/(delph(KAXIS))
                 end if
              end if

              BoxVal(iter) = BoxVal(iter) +  theta*(Condkmh/delmh(KAXIS)+Condkph/delph(KAXIS))*dt !! diag
#endif
#endif
              iter = iter + nentries

           end do
        end do
     end do

     call HYPRE_SStructMatrixSetBoxValues(gr_hypreMatA, mypart, gr_hypreLower(lb,1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), var, nentries, stencil_indices(1:nentries), BoxVal(:), ierr)

     deallocate(BoxVal)

 return

end subroutine gr_hypreCreateMatrix1Blk
