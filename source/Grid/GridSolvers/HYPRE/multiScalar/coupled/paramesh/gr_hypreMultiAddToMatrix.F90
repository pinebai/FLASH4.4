!!****if* source/Grid/GridSolvers/HYPRE/multiScalar/coupled/paramesh/gr_hypreMultiAddToMatrix
!!
!!  NAME 
!!
!!   gr_hypreMultiAddToMatrix
!!
!!  SYNOPSIS
!!
!!   call gr_hypreMultiAddToMatrix(integer(IN)           :: numVars,
!!                             integer(IN)           :: iFactorB,
!!                             integer(IN)           :: iFactorA,
!!                             integer(IN)           :: bcTypes(6),
!!                             real(IN)              :: bcValues(2,numGroups,6),
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
!!          B = MX, where M is a matrix whose product with th previous state produces RHS B.
!!
!!      A*(df/dt) + C*f = div(B*grad(f)) + D
!!      f -> Variable to be diffused.
!!      C,D are optional factors (not implemented here, the caller should add them later.)
!!
!!
!! ARGUMENTS
!! 
!!   numVars      : Number of variables on which the coupled diffusion operation is performed
!!   iFactorB     : factors in the equation with spatial variation.
!!   iFactorA     : 
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
!!   Currently, gr_hypreMultiAddToMatrix is called from Grid_advanceDiffusion with
!!   JacobiMatrix==.FALSE. only when the implicitness parameter theta passed to
!!   Grid_advanceDiffusion is 0. (KW 2012-12-05, corrected 2014-12-05)
!!
!!   At time of call, the data of baseVarDesc mst still be valid if thDfRHS .NE. 0.0.
!!
!! SEE ALSO
!!
!!  Grid_interface
!!***

!!REORDER(4): solnVec

subroutine gr_hypreMultiAddToMatrix(numVars, firstHypreVar, iFactorB, &
     absorpCoeffDesc, &
     emissCoeffDesc,emissTermDesc, &
     unkVarsDesc, &
     iFactorA, bcTypes, bcValues, dt, &
     alpha, blockCount, blockList, JacobiMatrix)
  
  use gr_hypreLocalInterface, ONLY: gr_hypreGetFaceBFcB, gr_hypreApplyBcToFace
  use gr_hypreData,     ONLY : gr_hypreSolverType,gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, gr_hypreVecB, &
                               gr_hypreRefineMIN, gr_hypreRefineMAX, gr_hypreNeghLevels, &
                               gr_hypreSurrBlkSum
  use gr_hypreMultiData,ONLY : gr_dfsvBaseVarDesc, &
                               gr_dfsvThetaC,gr_dfsvThetaD
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_ascGetBlkPtr, Grid_ascReleaseBlkPtr, &
    Grid_genGetBlkPtr, Grid_genReleaseBlkPtr, &
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
  integer, intent(IN) :: numVars
  integer, intent(IN) :: firstHypreVar
  integer, intent(IN) :: iFactorB
  integer,intent(in),dimension(:) :: absorpCoeffDesc
  integer,intent(in),dimension(:) :: emissCoeffDesc, emissTermDesc
  integer,intent(in),dimension(:) :: unkVarsDesc
  integer, intent(IN) :: iFactorA
  integer, intent(IN) :: bcTypes(6)
  real,    intent(IN) :: bcValues(:,:,:)
  real,    intent(IN) :: dt
  real,    intent(IN) :: alpha
  integer, intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  logical, intent(IN) :: JacobiMatrix
      
  !!-----------------------------------------------------------------------
  !!         LOCAL VARIABLES.
  !!-----------------------------------------------------------------------  
  integer :: iVar
  integer :: ierr, pos(NDIM)
  real, dimension(MDIM)     :: del , delph, delmh
  real, dimension(2*MDIM, MDIM) :: negh_del
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  integer :: datasize(MDIM), datasizeGC(MDIM)
  integer ::  mypart,mylevel
  integer ::  var, iv, hypreVar
  integer :: blockID
  real, allocatable :: cellVolumes(:,:,:)
  integer ::  nentries, nentriesSameVar, stencil_indices(8)  
  integer,allocatable ::  varIndices(:)
  real    ::  values(7), graph_value(7) 
  integer :: i, j, k, ii, lb
  real   :: condiph, condimh
  real   :: condjph, condjmh
  real   :: condkph, condkmh
  real   :: theta,thetaC,thetaD
  real   :: thCfRHS, thDfRHS
  integer, dimension(2,MDIM):: faces 
  integer :: eachNegh,numNegh 

  real, allocatable, dimension(:,:,:,:) :: xflux, yflux, zflux
  real :: dirichlet_multiplier
  integer :: dir
  real, allocatable :: faceAreas  (:,:,:,:)  
  real, POINTER, DIMENSION(:,:,:,:) :: facBptrX, facBptrY, facBptrZ
  real, POINTER, DIMENSION(:,:,:)   :: facBptr1
  real, POINTER, DIMENSION(:,:,:,:) :: facCptr,  facDptr,  facDrhsPtr,   unPtr, un0Ptr
  real, POINTER, DIMENSION(:,:,:)   :: facCptr1, facDptr1,  facDrhsPtr1, unPtr1,un0Ptr1
  real    ::  factorC, factorD,             un, un0
  integer :: iFactorC,iFactorD, iFactorDrhs
  character(len=32) :: matfile
  integer :: numGraph, iter, iter0, iter00
  real, allocatable :: BoxVal(:,:), BoxVal0(:), BoxVal00(:), RHSVal(:,:)
  integer, parameter :: fpv = 2**(NDIM-1) ! fluxes per variable
  logical,SAVE :: firstCall = .TRUE.
  
  call Timers_start("gr_hypreMultiAddToMatrix")  
  
  if(JacobiMatrix) then
     theta = alpha
     dirichlet_multiplier = 1.0
  else
     theta = alpha - 1.0
     dirichlet_multiplier = 0.0     
  end if
  
  iFactorC = absorpCoeffDesc(VARDESC_VAR)
  iFactorD = emissCoeffDesc(VARDESC_VAR)
  iFactorDrhs = emissTermDesc(VARDESC_VAR)
  iVar = unkVarsDesc(VARDESC_VAR)

  thetaC = gr_dfsvThetaC
  thetaD = gr_dfsvThetaD
  thCfRHS = (thetaC - theta) / theta
  thDfRHS = thetaD / theta
  un0 = 0.0
  if (firstCall) then
     print*,'thetaC,thetaD,thCfRHS,thDfRHS:'
     print*,thetaC,thetaD,thCfRHS,thDfRHS
     firstCall = .FALSE.
  end if

  nentriesSameVar = 2*NDIM + 1 
  nentries = nentriesSameVar + 1
  do i = 1, nentries
     stencil_indices(i) = i-1
  enddo
  allocate(varIndices(numVars))
  do i = 1, numVars
     varIndices(i) = firstHypreVar+i-1
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
     call Grid_genGetBlkPtr(blockID, facCptr,    absorpCoeffDesc)
     call Grid_genGetBlkPtr(blockID, facDptr,    emissCoeffDesc)
     call Grid_genGetBlkPtr(blockID, facDrhsPtr, emissTermDesc)
     call Grid_genGetBlkPtr(blockID, unPtr,      unkVarsDesc)
     if (thDfRHS .NE. 0.0) call Grid_genGetBlkPtr(blockID, un0Ptr,gr_dfsvBaseVarDesc)
     
     datasize  (1:MDIM)= blkLimits  (HIGH,1:MDIM)-blkLimits  (LOW,1:MDIM)+1
     datasizeGC(1:MDIM)= blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

     allocate(cellVolumes(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
     call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR,          &
          blkLimits(LOW,:), cellVolumes,datasize)      
     
     allocate(BoxVal(nentries*product(datasize(1:NDIM)),numVars))
     allocate(BoxVal0(numVars*product(datasize(1:NDIM))        ))
     allocate(BoxVal00(       product(datasize(1:NDIM))        ))
     allocate(RHSVal(         product(dataSize(1:NDIM)),0:numVars))
     
     mypart = mylevel - gr_hypreRefineMIN                  
     
     !!-----------------------------------------------------------------------
     !!         COMPUTE CELL FACE AREAS
     !!-----------------------------------------------------------------------     

     allocate(faceAreas(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),   &
                        blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),   &
                        blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS), NDIM))               
     
     call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, & 
          blkLimitsGC(LOW,:), faceAreas(:,:,:,IAXIS), datasizeGC)         
     
#if NDIM >= 2
     
     call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR, &
          blkLimitsGC(LOW,:), faceAreas(:,:,:,JAXIS), datasizeGC)        
     
#if NDIM == 3
     
     call Grid_getBlkData(blockID, CELL_FACEAREA, KLO_FACE, EXTERIOR, &
          blkLimitsGC(LOW,:), faceAreas(:,:,:,KAXIS), datasizeGC)

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
     iter0 = 1
     iter00 = 1
     BoxVal = 0.0
     BoxVal0(:) = 0.0
     BoxVal00(:) = 0.0
     RHSVal(:,:) = 0.0
     
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
              
              do var=1,numVars
                 iv = var-1
                 hypreVar = firstHypreVar + iv
                 !! i-1,j,k
                 if ((i /= blkLimits(LOW, IAXIS)) .or. (faces(1,IAXIS) == NOT_BOUNDARY)) then                                         
                    call AssoMed(facBptr1,facBptrX,iFactorB+iv) ! facBptr1(:,:,:) => facBptrX(iFactorB+iv,:,:,:)
                    if (i == blkLimits(LOW, IAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,LEFT_EDGE, 1+K2D, 1+K3D)) then  !! F/C boundary.                                         
                       numNegh = gr_hypreSurrBlkSum(lb) % regionInfo(LEFT_EDGE,1+K2D,1+K3D) % numNegh    
                       numGraph = numNegh                    
                       do eachNegh = 1, numNegh                                                            
                          ii = nentries+eachNegh-1                       
                          delmh(IAXIS) =  0.5*(negh_del(1,IAXIS) + del(IAXIS))

                          if (numNegh > 1) then !! we are on coarse cell, use the average computed from fine cell (2 of them)
                             condimh = xflux(eachNegh+iv*fpv,i,j,k)
                          else
                             if (NDIM > 1) then
                                !! we are on fine cell and looking at coarse block, use regular averaging.
                                condimh = facBptr1(i,j,k)*faceAreas(i,j,k,IAXIS)
                             else
                                condimh = xflux(eachNegh+iv*fpv,i,j,k)
                             end if
                          end if

                          graph_value(1) = -condimh*theta*dt/delmh(IAXIS)
                          call HYPRE_SStructMatrixSetValues(gr_hypreMatA, mypart, pos(1:NDIM), hypreVar, 1, ii, graph_value,ierr)
                       end do

                       if (numNegh > 1) then
                          condimh = sum(xflux(1+iv*fpv:var*fpv,i,j,k))
                       end if
                    else
                       !! STENCILED RELATIONSHIP.
                       delmh(IAXIS) = del(IAXIS)
                       condimh = facBptr1(i,j,k)*faceAreas(i,j,k,IAXIS)
                       BoxVal(iter+1,var) = -condimh*theta*dt/delmh(IAXIS)
                    end if
                 end if

              !! i+1,j,k
              if (i /= blkLimits(HIGH, IAXIS) .or. (faces(2,IAXIS) == NOT_BOUNDARY)) then
                 call AssoMed(facBptr1,facBptrX,iFactorB) ! facBptr1(:,:,:) => facBptrX(iFactorB,:,:,:)
                 if (i == blkLimits(HIGH, IAXIS) .and. mylevel /= gr_hypreNeghLevels(lb, RIGHT_EDGE, 1+K2D, 1+K3D)) then !! F/C boundary.
                    numNegh = gr_hypreSurrBlkSum(lb) % regionInfo(RIGHT_EDGE, 1+K2D, 1+K3D) % numNegh
                    numGraph = numNegh
                    do eachNegh = 1, numNegh
                       ii= nentries+eachNegh-1
                       delph(IAXIS) =  0.5*(negh_del(2,IAXIS) + del(IAXIS))
                       if (numNegh > 1) then !! we are on coarse cell, use fine cell area on fluxes.
                          condiph = xflux(eachNegh+iv*fpv,i+1,j,k)
                       else
                          if (NDIM > 1) then
                             condiph = facBptr1(i+1,j,k)*faceAreas(i+1,j,k,IAXIS)
                          else
                             condiph = xflux(eachNegh+iv*fpv,i+1,j,k)
                          end if
                       end if
                       graph_value(1) =  -condiph*theta*dt/delph(IAXIS)
                       call HYPRE_SStructMatrixSetValues(gr_hypreMatA, mypart, pos(1:NDIM), hypreVar, 1, ii, graph_value,ierr)
                    end do

                    if (numNegh > 1) then
                       condiph = sum(xflux(1+iv*fpv:var*fpv,i+1,j,k)) !! actual flux, goes towards computing diag.
                    end if
                 else
                    !! STENCILED RELATIONSHIP.
                    delph(IAXIS) =  del(IAXIS)
                    condiph = facBptr1(i+1,j,k)*faceAreas(i+1,j,k,IAXIS)
                    BoxVal(iter+2,var) =  -condiph*theta*dt/(delph(IAXIS))
                 end if
              end if

              BoxVal(iter,var) = BoxVal(iter,var) + theta*((Condimh/delmh(1))+(Condiph/delph(1)))*dt

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
                       ii=nentries+eachNegh-1+numGraph
                       if (numNegh > 1) then !! we are on coarse cell, use fine cell area on fluxes.
                          condjmh = yflux(eachNegh+iv*fpv,i,j,k)
                       else
                          condjmh = facBptr1(i,j,k)*faceAreas(i,j,k,JAXIS)
                       end if
                       graph_value(1) = -condjmh*theta*dt/delmh(JAXIS)
                       call HYPRE_SStructMatrixSetValues(gr_hypreMatA, mypart, pos(1:NDIM), hypreVar, 1, ii, graph_value,ierr)
                    end do

                    numGraph = numGraph + numNegh

                    if (numNegh > 1) then
                       condjmh = sum(yflux(1+iv*fpv:var*fpv,i,j,k))
                    end if

                 else
                    !! STENCILED RELATIONSHIP.
                    delmh(JAXIS) = del(JAXIS)
                    condjmh = facBptr1(i,j,k)*faceAreas(i,j,k,JAXIS)
                    BoxVal(iter+3,var) = -condjmh*theta*dt/(delmh(JAXIS))
                 end if
              end if

              !! i,j+1,k
              if ((j /= blkLimits(HIGH, JAXIS)) .or. (faces(2,JAXIS) == NOT_BOUNDARY)) then
                 call AssoMed(facBptr1,facBptrY,iFactorB) ! facBptr1(:,:,:) => facBptrY(iFactorB,:,:,:)
                 if (j ==  blkLimits(HIGH, JAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,1+K1D,RIGHT_EDGE,1+K3D)) then !! F/C boundary.
                    numNegh = gr_hypreSurrBlkSum(lb) % regionInfo(1+K1D,RIGHT_EDGE,1+K3D) % numNegh
                    do eachNegh = 1, numNegh
                       delph(JAXIS) =  0.5*(negh_del(4,JAXIS) + del(JAXIS))
                       ii=nentries+eachNegh-1+numGraph
                       if (numNegh > 1) then !! we are on coarse cell, use fine cell area on fluxes.
                          condjph = yflux(eachNegh+iv*fpv,i,j+1,k)
                       else
                          condjph = facBptr1(i,j+1,k)*faceAreas(i,j+1,k,JAXIS)
                       end if
                       graph_value(1) = -condjph*theta*dt/delph(JAXIS)
                       call HYPRE_SStructMatrixSetValues(gr_hypreMatA, mypart, pos(1:NDIM), hypreVar, 1, ii, graph_value,ierr)
                    end do

                    numGraph = numGraph + numNegh

                    if (numNegh > 1) then
                       condjph = sum(yflux(1+iv*fpv:var*fpv,i,j+1,k))
                    end if
                 else
                    !! STENCILED RELATIONSHIP.
                    delph(JAXIS) = del(JAXIS)
                    condjph = facBptr1(i,j+1,k)*faceAreas(i,j+1,k,JAXIS)
                    BoxVal(iter+4,var) = -condjph*theta*dt/(delph(JAXIS))
                 end if
              end if

              BoxVal(iter,var) = BoxVal(iter,var) +  theta*(Condjmh/delmh(2)+Condjph/delph(2))*dt !! diag

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
                       ii=nentries+eachNegh-1+numGraph
                       if (numNegh > 1) then !! we are on coarse cell, use fine cell area on fluxes.
                          condkmh = zflux(eachNegh+iv*fpv,i,j,k)
                       else
                          condkmh = facBptr1(i,j,k)*faceAreas(i,j,k,KAXIS)
                       end if
                       graph_value(1) = -condkmh*theta*dt/delmh(KAXIS)
                       call HYPRE_SStructMatrixSetValues(gr_hypreMatA, mypart, pos(1:NDIM), hypreVar, 1, ii, graph_value,ierr)
                    end do

                    numGraph = numGraph + numNegh

                    if (numNegh > 1) then
                       condkmh = sum(zflux(1+iv*fpv:var*fpv,i,j,k))
                    end if

                 else
                    !! STENCILED RELATIONSHIP.
                    delmh(KAXIS) = del(KAXIS)
                    condkmh = facBptr1(i,j,k)*faceAreas(i,j,k,KAXIS)
                    BoxVal(iter+5,var) = -condkmh*theta*dt/(delmh(KAXIS))
                 end if
              end if

              !! i,j,k+1
              if ((k /= blkLimits(HIGH,KAXIS)) .or. (faces(2,KAXIS) == NOT_BOUNDARY)) then
                 call AssoMed(facBptr1,facBptrZ,iFactorB) ! facBptr1(:,:,:) => facBptrZ(iFactorB,:,:,:)
                 if (k ==  blkLimits(HIGH,KAXIS) .and. mylevel /= gr_hypreNeghLevels(lb,1+K1D,1+K2D,RIGHT_EDGE)) then !! F/C boundary.
                    numNegh = gr_hypreSurrBlkSum(lb) % regionInfo(1+K1D,1+K2D,RIGHT_EDGE) % numNegh
                    do eachNegh = 1, numNegh
                       delph(KAXIS) =  0.5*(negh_del(6,KAXIS) + del(KAXIS))
                       ii=nentries+eachNegh-1+numGraph
                       if (numNegh > 1) then !! we are on coarse cell, use fine cell area on fluxes.
                          condkph = zflux(eachNegh+iv*fpv,i,j,k+1)
                       else
                          condkph = facBptr1(i,j,k+1)*faceAreas(i,j,k+1,KAXIS)
                       end if
                       graph_value(1) = -condkph*theta*dt/delph(KAXIS)
                       call HYPRE_SStructMatrixSetValues(gr_hypreMatA, mypart, pos(1:NDIM), hypreVar, 1, ii, graph_value,ierr)
                    end do

                    numGraph = numGraph + numNegh

                    if (numNegh > 1) then
                       condkph = sum(zflux(1+iv*fpv:var*fpv,i,j,k+1))
                    end if
                 else
                    !! STENCILED RELATIONSHIP.
                    delph(KAXIS) = del(KAXIS)
                    condkph = facBptr1(i,j,k+1)*faceAreas(i,j,k+1,KAXIS)
                    BoxVal(iter+6,var) = -condkph*theta*dt/(delph(KAXIS))
                 end if
              end if

                 BoxVal(iter,var) = BoxVal(iter,var) +  theta*(Condkmh/delmh(KAXIS)+Condkph/delph(KAXIS))*dt !! diag
#endif
#endif
                 call AssoMed(facCptr1,facCptr,iFactorC+iv) ! facCptr1(:,:,:) => facCptr(iFactorC+iv,:,:,:)
                 factorC = facCptr1(i,j,k) * cellVolumes(i,j,k)
!*****
!!$                 BoxVal(iter+nentriesSameVar,var) = - thetaC*(factorC)*dt                !! 0 <- var  coupling
                 BoxVal0(iter0+iv) = - thetaC*(factorC)*dt                !! 0 <- var  coupling
!*****

!!$                 BoxVal(iter,var) = BoxVal(iter,var)  - BoxVal(iter+nentriesSameVar,var) !! mod. diag coupling
                 BoxVal(iter,var) = BoxVal(iter,var)  + thetaC*(factorC)*dt

                 call AssoMed(facDptr1,facDptr,iFactorD+iv) ! facDptr1(:,:,:) => facDptr(iFactorD+iv,:,:,:)
                 factorD = facDptr1(i,j,k) * cellVolumes(i,j,k)       !factorD ...
!*****
!!$                 BoxVal0(iter0+iv) =      - thetaD*(factorD)*dt             !!     var <- 0  coupling
                 BoxVal(iter+nentriesSameVar,var) =      - thetaD*(factorD)*dt             !!     var <- 0  coupling
!*****
!!$                 BoxVal00(iter00) =  BoxVal00(iter00) - BoxVal0(iter0+iv)   !! mod. 0  <- 0  coupling
                 BoxVal00(iter00) =  BoxVal00(iter00) + thetaD*(factorD)*dt    !! mod. 0  <- 0  coupling

                 call AssoMed(unPtr1,  unPtr,     iVar+iv)        ! unPtr1(:,:,:)   => unPtr     (iVar+iv,:,:,:)
                 un      = unPtr1(i,j,k)   * cellVolumes(i,j,k)       !U^(n) ...
                 call AssoMed(un0Ptr1, un0Ptr,    gr_dfsvBaseVarDesc(VARDESC_VAR)) ! un0Ptr1(:,:,:)   => baseVar
                 if (thDfRHS .NE. 0.0) un0     = un0Ptr1(i,j,k)  * cellVolumes(i,j,k)       !U_0^(n) ...
                 call AssoMed(facDrhsPtr1,facDrhsPtr,iFactorDrhs+iv) ! facDrhsPtr1(:,:,:) => facDrhsPtr(iFactorDrhs+iv,:,:,:)
                 factorD = (-thDfRHS*facDptr1(i,j,k)*un0) + facDrhsPtr1(i,j,k) * cellVolumes(i,j,k)       !termD, etc!
                 RHSVal(iter00,var) =              (factorD + thCfRHS*facCptr1(i,j,k)*un)*dt             !!     RHS(var)
!!$                 print*,blockID,i,un,un0,&
!!$                      -thDfRHS*facDptr1(i,j,k)*un0,&
!!$                      facDrhsPtr1(i,j,k) * cellVolumes(i,j,k),&
!!$                      +thCfRHS*facCptr1(i,j,k)*un,&
!!$                      RHSVal(iter00,var)  / dt
                 RHSVal(iter00,0) =  RHSVal(iter00,0) - RHSVal(iter00,var)  !! mod. RHS(0)
                 RHSVal(iter00,var) = RHSVal(iter00,var) + un               !!     RHS(var) += _1_

              end do
              iter = iter + nentries
              iter0 = iter0 + numVars
              iter00 = iter00 + 1
           end do
        end do
     end do

     do var=1,numVars
        iv = var-1
        hypreVar = firstHypreVar + iv
        call HYPRE_SStructMatrixSetBoxValues(gr_hypreMatA, mypart, gr_hypreLower(lb,1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), hypreVar, nentries, stencil_indices(1:nentries), BoxVal(:,var), ierr)
     end do
     call HYPRE_SStructMatrixSetBoxValues(gr_hypreMatA, mypart, gr_hypreLower(lb,1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), 0, numVars, varIndices(1:numVars), BoxVal0, ierr)
     call HYPRE_SStructMatrixAddToBoxValu(gr_hypreMatA, mypart, gr_hypreLower(lb,1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), 0, 1, stencil_indices(1), BoxVal00, ierr)

     deallocate (BoxVal)
     deallocate (BoxVal0)
     deallocate (BoxVal00)
     deallocate(cellVolumes)

     do var=1,numVars
        iv = var-1
        hypreVar = firstHypreVar + iv
        call HYPRE_SStructVectorAddToBoxValu(gr_hypreVecB, mypart, gr_hypreLower(lb, 1:NDIM), &
             gr_hypreUpper(lb,1:NDIM), hypreVar, RHSVal(:,var), ierr)
     end do
     call HYPRE_SStructVectorAddToBoxValu(gr_hypreVecB, mypart, gr_hypreLower(lb, 1:NDIM), &
             gr_hypreUpper(lb,1:NDIM), 0, RHSVal(:,0), ierr)
     deallocate (RHSVal)

     call Timers_start("gr_hypreApplyBcToFace")
     dir = ILO_FACE
     do i = IAXIS, NDIM
        do j = LOW, HIGH
           if (faces(j,i) /= NOT_BOUNDARY) then
              do var=1,numVars
                 iv = var-1
                 hypreVar = firstHypreVar + iv
                 call gr_hypreApplyBcToFace(blkLimits,blkLimitsGC,mypart,hypreVar,iFactorB+iv,bcTypes(dir),dir, &
                      bcValues(:,var,dir), dt, theta, del(i), gr_hypreLower(lb,:), &
                      dirichlet_multiplier, faceAreas(:,:,:,i), solnVec, &
                      blockID, facBptrX, facBptrY, facBptrZ)
              end do
           end if
           dir = dir + 1
        end do
     end do
     call Timers_stop("gr_hypreApplyBcToFace")

     call Grid_releaseBlkPtr(blockID, solnVec)
     deallocate (faceAreas)
     call Grid_ascReleaseBlkPtr(blockID,facBptrX,FACEX)
#if NDIM >= 2
     call Grid_ascReleaseBlkPtr(blockID,facBptrY,FACEY)
#if NDIM == 3
     call Grid_ascReleaseBlkPtr(blockID,facBptrZ,FACEZ)
#endif
#endif     
     call Grid_genReleaseBlkPtr(blockID, facCptr,    absorpCoeffDesc)
     call Grid_genReleaseBlkPtr(blockID, facDptr,    emissCoeffDesc)
     call Grid_genReleaseBlkPtr(blockID, facDrhsPtr, emissTermDesc)
     call Grid_genReleaseBlkPtr(blockID, unPtr,      unkVarsDesc)
     if (thDfRHS .NE. 0.0) call Grid_genReleaseBlkPtr(blockID, un0Ptr,gr_dfsvBaseVarDesc)

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
 deallocate(varIndices)

 call Timers_stop("gr_hypreMultiAddToMatrix")

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
end subroutine gr_hypreMultiAddToMatrix
