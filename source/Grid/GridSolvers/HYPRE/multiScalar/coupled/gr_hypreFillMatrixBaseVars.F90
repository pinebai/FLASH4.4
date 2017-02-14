!!****if* source/Grid/GridSolvers/HYPRE/multiScalar/coupled/gr_hypreFillMatrixBaseVars
!!
!!  NAME 
!!
!!   gr_hypreFillMatrixBaseVars
!!
!!  SYNOPSIS
!!
!!   call gr_hypreFillMatrixBaseVars(integer(IN)           :: iVar,
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

subroutine gr_hypreFillMatrixBaseVars(iVar, iFactorB, iFactorA, bcTypes, bcValues, dt, &
     alpha, blockCount, blockList, JacobiMatrix)
  
  use gr_hypreLocalInterface, ONLY: gr_hypreGetFaceBFcB, gr_hypreApplyBcToFace
  use gr_hypreData,     ONLY : gr_hypreSolverType,gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, &
                               gr_hypreRefineMIN, gr_hypreRefineMAX, gr_hypreNeghLevels, &
                               gr_hypreSurrBlkSum
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_ascGetBlkPtr, Grid_ascReleaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_fillGuardCells, Grid_getBlkBC, &
    Grid_getBlkCornerID, Grid_getCellCoords, Grid_getFluxData, Grid_getBlkData, Grid_getDeltas, Grid_getBlkRefineLevel
  use Timers_interface, ONLY : Timers_start, Timers_stop 
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
                               GRID_PDE_BND_NEUMANN,   &
                               GRID_PDE_BND_DIRICHLET
  
  
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
  real,    intent(IN) :: alpha
  integer, intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  logical, intent(IN) :: JacobiMatrix
      
  !!-----------------------------------------------------------------------
  !!         LOCAL VARIABLES.
  !!-----------------------------------------------------------------------  
  integer :: ierr, pos(NDIM)
  real, dimension(MDIM)     :: del , delph, delmh
  real, dimension(2*MDIM, MDIM) :: negh_del
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  integer :: datasize(MDIM), datasizeGC(MDIM)
  integer ::  mypart,mylevel
  integer ::  var
  integer :: blockID
  integer ::  nentries,  stencil_indices(7)  
  real    ::  values(7), graph_value(7) 
  integer :: i, j, k, ii, lb
  real   :: condiph, condimh
  real   :: condjph, condjmh
  real   :: condkph, condkmh
  real   :: theta
  integer, dimension(2,MDIM):: faces 
  integer :: eachNegh,numNegh 

  real, allocatable, dimension(:,:,:,:) :: xflux, yflux, zflux
  real :: dirichlet_multiplier
  integer :: dir
  real, allocatable :: faceAreas  (:,:,:,:)  
  real, POINTER, DIMENSION(:,:,:,:) :: facBptrX, facBptrY, facBptrZ
  real, POINTER, DIMENSION(:,:,:)   :: facBptr1
  character(len=32) :: matfile
  integer :: numGraph, iter
  real, allocatable :: BoxVal(:)
  integer, parameter :: nFluxVars = 2**(NDIM-1)
  
  call Timers_start("gr_hypreFillMatrixBaseVars")  
  
  if(JacobiMatrix) then
     theta = alpha
     dirichlet_multiplier = 1.0
  else
     theta = alpha - 1.0
     dirichlet_multiplier = 0.0     
  end if
  
  nentries = 2*NDIM + 1 
  do i = 1, nentries
     stencil_indices(i) = i-1
  enddo

  nullify(facBptrY)
  nullify(facBptrZ)

  mypart = 0  !! part iterator 
  var    = 0  !! var iterator.
  
  if (blockCount > 0) then
     call Grid_getBlkIndexLimits(blockList(1),blkLimits,blkLimitsGC)         
     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1          
     allocate(xflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(yflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(zflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
  else
     datasize = 0.
  end if
  


  if(blockCount == 0) then
     call Timers_start("gr_hypreApplyBcToFace")     
     call Timers_stop("gr_hypreApplyBcToFace")     
  end if
  
  do lb = 1, blockCount     
     blockID = blockList(lb)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     call Grid_getBlkPtr(blockID, solnVec)
     call Grid_getDeltas(blockID, del)
     call Grid_getBlkBC (blockID, faces)     
     call Grid_getBlkRefineLevel(blockID,mylevel)

     call Grid_ascGetBlkPtr(blockID,facBptrX,FACEX)
#if NDIM >= 2
     call Grid_ascGetBlkPtr(blockID,facBptrY,FACEY)
#if NDIM == 3
     call Grid_ascGetBlkPtr(blockID,facBptrZ,FACEZ)
#endif
#endif     

     datasize  (1:MDIM)= blkLimits  (HIGH,1:MDIM)-blkLimits  (LOW,1:MDIM)+1
     datasizeGC(1:MDIM)= blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
     
     allocate(BoxVal(nentries*product(datasize(1:NDIM))))
     
     mypart = mylevel - gr_hypreRefineMIN                  
     
     !!-----------------------------------------------------------------------
     !!         COMPUTE CELL FACE AREAS
     !!-----------------------------------------------------------------------     

     allocate(faceAreas(NDIM, blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),   &
                              blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),   &
                              blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))               
     
     call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, & 
          blkLimitsGC(LOW,:), faceAreas(IAXIS,:,:,:), datasizeGC)         
     
#if NDIM >= 2
     
     call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR, &
          blkLimitsGC(LOW,:), faceAreas(JAXIS,:,:,:), datasizeGC)        
     
#if NDIM == 3
     
     call Grid_getBlkData(blockID, CELL_FACEAREA, KLO_FACE, EXTERIOR, &
          blkLimitsGC(LOW,:), faceAreas(KAXIS,:,:,:), datasizeGC)

#endif
#endif
     
     if (mylevel > gr_hypreNeghLevels(lb,LEFT_EDGE,1+K2D,1+K3D)) then !!F/C
        negh_del(1,:) = del(:)*2.0
     else !! C/F
        negh_del(1,:) = del(:)/2.0
     end if
     
     if (mylevel > gr_hypreNeghLevels(lb,RIGHT_EDGE,1+K2D,1+K3D)) then !!F/C
        negh_del(2,:) = del(:)*2.0
     else !! C/F
        negh_del(2,:) = del(:)/2.0
     end if
     
     if (mylevel > gr_hypreNeghLevels(lb,1+K1D,LEFT_EDGE,1+K3D)) then !!F/C
        negh_del(3,:) = del(:)*2.0
     else !! C/F
        negh_del(3,:) = del(:)/2.0
     end if
     
     if (mylevel > gr_hypreNeghLevels(lb,1+K1D,RIGHT_EDGE,1+K3D)) then !!F/C
        negh_del(4,:) = del(:)*2.0
     else !! C/F
        negh_del(4,:) = del(:)/2.0
     end if
     
     if (mylevel > gr_hypreNeghLevels(lb,1+K1D,1+K2D,LEFT_EDGE)) then !!F/C
        negh_del(5,:) = del(:)*2.0
     else !! C/F
        negh_del(5,:) = del(:)/2.0
     end if
     
     if (mylevel > gr_hypreNeghLevels(lb,1+K1D,1+K2D,RIGHT_EDGE)) then !!F/C
        negh_del(6,:) = del(:)*2.0
     else !! C/F
        negh_del(6,:) = del(:)/2.0
     end if

     ! This assumes that calls to Grid_putFluxData and Grid_conserveFluxes have
     ! happened before.  Currently this is done in gr_hypreExchangeFacB. KW 2012-12-04
     call Grid_getFluxData(blockID,IAXIS,xflux,datasizeGC)
#if NDIM >= 2
     call Grid_getFluxData(blockID,JAXIS,yflux,datasizeGC)
#if NDIM == 3
     call Grid_getFluxData(blockID,KAXIS,zflux,datasizeGC)
#endif
#endif     
     
     iter = 1
     BoxVal = 0.0
     
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
              if ((i /= blkLimits(LOW, IAXIS)) .or. (faces(1,IAXIS) == NOT_BOUNDARY)) then                                         
                 call AssoMed(facBptr1,facBptrX,iFactorB) ! facBptr1(:,:,:) => facBptrX(iFactorB,:,:,:)
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
                             condimh = facBptr1(i,j,k)*faceAreas(IAXIS,i,j,k)
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
                    condimh = facBptr1(i,j,k)*faceAreas(IAXIS,i,j,k)
                    BoxVal(iter+1) = -condimh*theta*dt/delmh(IAXIS)
                 end if
              end if

              !! i+1,j,k
              if (i /= blkLimits(HIGH, IAXIS) .or. (faces(2,IAXIS) == NOT_BOUNDARY)) then
                 call AssoMed(facBptr1,facBptrX,iFactorB) ! facBptr1(:,:,:) => facBptrX(iFactorB,:,:,:)
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
                             condiph = facBptr1(i+1,j,k)*faceAreas(IAXIS,i+1,j,k)
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
                    condiph = facBptr1(i+1,j,k)*faceAreas(IAXIS,i+1,j,k)
                    BoxVal(iter+2) =  -condiph*theta*dt/(delph(IAXIS))
                 end if
              end if

              BoxVal(iter) = BoxVal(iter) + theta*((Condimh/delmh(1))+(Condiph/delph(1)))*dt

#if NDIM >= 2
              condjmh = 0.
              condjph = 0.

              !! i,j-1,k
              if ((j /= blkLimits(LOW, JAXIS)) .or. (faces(1,JAXIS) == NOT_BOUNDARY)) then
                 call AssoMed(facBptr1,facBptrY,iFactorB) ! facBptr1(:,:,:) => facBptrY(iFactorB,:,:,:)
                 if (j == blkLimits(LOW, JAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,1+K1D,LEFT_EDGE,1+K3D)) then  !! F/C boundary.
                    numNegh = gr_hypreSurrBlkSum(lb) % regionInfo(1+K1D,LEFT_EDGE,1+K3D) % numNegh
                    do eachNegh = 1, numNegh
                       delmh(JAXIS) =  0.5*(negh_del(3,JAXIS) + del(JAXIS))
                       ii=2*NDIM+1+eachNegh-1+numGraph
                       if (numNegh > 1) then !! we are on coarse cell, use fine cell area on fluxes.
                          condjmh = yflux(eachNegh,i,j,k)
                       else
                          condjmh = facBptr1(i,j,k)*faceAreas(JAXIS,i,j,k)
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
                    condjmh = facBptr1(i,j,k)*faceAreas(JAXIS,i,j,k)
                    BoxVal(iter+3) = -condjmh*theta*dt/(delmh(JAXIS))
                 end if
              end if

              !! i,j+1,k
              if ((j /= blkLimits(HIGH, JAXIS)) .or. (faces(2,JAXIS) == NOT_BOUNDARY)) then
                 call AssoMed(facBptr1,facBptrY,iFactorB) ! facBptr1(:,:,:) => facBptrY(iFactorB,:,:,:)
                 if (j ==  blkLimits(HIGH, JAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,1+K1D,RIGHT_EDGE,1+K3D)) then !! F/C boundary.
                    numNegh = gr_hypreSurrBlkSum(lb) % regionInfo(1+K1D,RIGHT_EDGE,1+K3D) % numNegh
                    do eachNegh = 1, numNegh
                       delph(JAXIS) =  0.5*(negh_del(4,JAXIS) + del(JAXIS))
                       ii=2*NDIM+1+eachNegh-1+numGraph
                       if (numNegh > 1) then !! we are on coarse cell, use fine cell area on fluxes.
                          condjph = yflux(eachNegh,i,j+1,k)
                       else
                          condjph = facBptr1(i,j+1,k)*faceAreas(JAXIS,i,j+1,k)
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
                    condjph = facBptr1(i,j+1,k)*faceAreas(JAXIS,i,j+1,k)
                    BoxVal(iter+4) = -condjph*theta*dt/(delph(JAXIS))
                 end if
              end if

              BoxVal(iter) = BoxVal(iter) +  theta*(Condjmh/delmh(2)+Condjph/delph(2))*dt !! diag

#if NDIM == 3
              condkmh = 0.
              condkph = 0.
              !! i,j,k-1
              if ((k /= blkLimits(LOW, KAXIS)) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) then
                 call AssoMed(facBptr1,facBptrZ,iFactorB) ! facBptr1(:,:,:) => facBptrZ(iFactorB,:,:,:)
                 if (k == blkLimits(LOW, KAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,1+K1D,1+K2D,LEFT_EDGE)) then  !! F/C boundary.
                    numNegh = gr_hypreSurrBlkSum(lb) % regionInfo(1+K1D,1+K2D,LEFT_EDGE) % numNegh
                    do eachNegh = 1, numNegh
                       delmh(KAXIS) =  0.5*(negh_del(5,KAXIS) + del(KAXIS))
                       ii=2*NDIM+1+eachNegh-1+numGraph
                       if (numNegh > 1) then !! we are on coarse cell, use fine cell area on fluxes.
                          condkmh = zflux(eachNegh,i,j,k)
                       else
                          condkmh = facBptr1(i,j,k)*faceAreas(KAXIS,i,j,k)
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
                    condkmh = facBptr1(i,j,k)*faceAreas(KAXIS,i,j,k)
                    BoxVal(iter+5) = -condkmh*theta*dt/(delmh(KAXIS))
                 end if
              end if

              !! i,j,k+1
              if ((k /= blkLimits(HIGH,KAXIS)) .or. (faces(2,KAXIS) == NOT_BOUNDARY)) then
                 call AssoMed(facBptr1,facBptrZ,iFactorB) ! facBptr1(:,:,:) => facBptrZ(iFactorB,:,:,:)
                 if (k ==  blkLimits(HIGH,KAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,1+K1D,1+K2D,RIGHT_EDGE)) then !! F/C boundary.
                    numNegh = gr_hypreSurrBlkSum(lb) % regionInfo(1+K1D,1+K2D,RIGHT_EDGE) % numNegh
                    do eachNegh = 1, numNegh
                       delph(KAXIS) =  0.5*(negh_del(6,KAXIS) + del(KAXIS))
                       ii=2*NDIM+1+eachNegh-1+numGraph
                       if (numNegh > 1) then !! we are on coarse cell, use fine cell area on fluxes.
                          condkph = zflux(eachNegh,i,j,k+1)
                       else
                          condkph = facBptr1(i,j,k+1)*faceAreas(KAXIS,i,j,k+1)
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
                    condkph = facBptr1(i,j,k+1)*faceAreas(KAXIS,i,j,k+1)
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


     call Timers_start("gr_hypreApplyBcToFace")
     dir = ILO_FACE
     do i = IAXIS, NDIM
        do j = LOW, HIGH
           if (faces(j,i) /= NOT_BOUNDARY) then
              call gr_hypreApplyBcToFace(blkLimits,blkLimitsGC,mypart,var,iFactorB,bcTypes(dir),dir, &
                   bcValues(:,dir), dt, theta, del(i), gr_hypreLower(lb,:), dirichlet_multiplier, faceAreas(i,:,:,:), solnVec, &
                   blockID, facBptrX, facBptrY, facBptrZ)
           end if
           dir = dir + 1
        end do
     end do
     call Timers_stop("gr_hypreApplyBcToFace")

     call Grid_releaseBlkPtr(blockID, solnVec)
     deallocate (faceAreas)
     deallocate (BoxVal)
     call Grid_ascReleaseBlkPtr(blockID,facBptrX,FACEX)
#if NDIM >= 2
     call Grid_ascReleaseBlkPtr(blockID,facBptrY,FACEY)
#if NDIM == 3
     call Grid_ascReleaseBlkPtr(blockID,facBptrZ,FACEZ)
#endif
#endif     

  end do !! block

  !!-----------------------------------------------------------------------
 !!         THIS IS A GLOBAL CALL.
 !!-----------------------------------------------------------------------
 call HYPRE_SStructMatrixAssemble(gr_hypreMatA, ierr)

 if (blockCount > 0) then
    deallocate (xflux)
    deallocate (yflux)
    deallocate (zflux)
 end if

 call Timers_stop("gr_hypreFillMatrixBaseVars")

 return

contains
#include "FortranLangFeatures.fh"
  subroutine AssoMed(pp, mm, varNo)
    real,POINTER_INTENT_OUT :: pp(:,:,:)
    real,POINTER_INTENT_IN  :: mm(:,:,:,:)
    integer,intent(in) :: varNo
    call AssoFin(pp,mm(varNo,:,:,:),lbound(mm,1),lbound(mm,2),lbound(mm,3),lbound(mm,4))
  end subroutine AssoMed

  subroutine AssoFin(pp, dd, lb1,lb2,lb3,lb4)
    real,POINTER_INTENT_OUT :: pp(:,:,:)
    integer, intent(in) :: lb1,lb2,lb3,lb4
    real,   intent(in),target :: dd(lb2:,lb3:,lb4:)
    pp => dd
  end subroutine AssoFin
end subroutine gr_hypreFillMatrixBaseVars
