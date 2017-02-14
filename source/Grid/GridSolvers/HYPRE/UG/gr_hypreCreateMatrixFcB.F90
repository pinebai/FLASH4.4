!!****if* source/Grid/GridSolvers/HYPRE/UG/gr_hypreCreateMatrixFcB
!!
!!  NAME 
!!
!!   gr_hypreCreateMatrixFcB
!!
!!  SYNOPSIS
!!
!!   call gr_hypreCreateMatrixFcB(integer(IN)           :: iVar,
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
!!   JacobiMatrix : TRUE computes A; FALSE computes M.
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

subroutine gr_hypreCreateMatrixFcB(iVar, iFactorB, iFactorA, bcTypes, bcValues, dt, &
     alpha, blockCount, blockList, JacobiMatrix)
  
  use gr_hypreLocalInterface, ONLY: gr_hypreApplyBcToFace
  use gr_hypreData,     ONLY : gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, &
                               gr_hypreAnisoDiffusion
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_ascGetBlkPtr, Grid_ascReleaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_fillGuardCells, Grid_getBlkBC, &
    Grid_getBlkCornerID, Grid_getCellCoords, Grid_getBlkData, Grid_getDeltas
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
  integer :: ierr
  real, dimension(MDIM)     :: del
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  integer :: datasize(MDIM), datasizeGC(MDIM)
  integer, parameter ::  mypart = 0  !! HYPRE part
  integer ::  var
  integer ::  blockID
  integer ::  nentries, stencil_indices(19)
  integer :: i, j, k,  lb
  integer, dimension(2,MDIM):: faces 
  real :: condimh, condiph
  real :: condjmh, condjph
  real :: condkmh, condkph
  real   :: theta
  real :: dirichlet_multiplier  
  integer :: dir, ii
  real, allocatable :: faceAreas  (:,:,:,:)   
  real, POINTER, DIMENSION(:,:,:,:) :: facBptrX, facBptrY, facBptrZ
  real, POINTER, DIMENSION(:,:,:)   :: facBptr1
  real, allocatable :: BoxVal(:)
  
  call Timers_start("gr_hypreCreateMatrixFcB") 
  
  if(JacobiMatrix) then
     theta = alpha
     dirichlet_multiplier = 1.0
  else
     theta = alpha - 1.0
     dirichlet_multiplier = 0.0     
  end if
  
  if (.not. gr_hypreAnisoDiffusion .OR. (NDIM == 1)) then
     nentries = NDIM*2 + 1
  else
     if (NDIM == 2) then
        nentries = 9
     elseif (NDIM == 3) then
        nentries = 19 !(27-8)
     endif
  endif

  do i = 1, nentries
     stencil_indices(i) = i-1
  enddo

  nullify(facBptrY)
  nullify(facBptrZ)

  var    = 0  !! var iterator.
  
  do lb = 1, blockCount 
     
     blockID = blockList(lb)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     call Grid_getBlkPtr(blockID, solnVec)
     call Grid_getDeltas(blockID, del)
     call Grid_getBlkBC (blockID, faces)     

     call Grid_ascGetBlkPtr(blockID,facBptrX,FACEX)
#if NDIM >= 2
     call Grid_ascGetBlkPtr(blockID,facBptrY,FACEY)
#if NDIM == 3
     call Grid_ascGetBlkPtr(blockID,facBptrZ,FACEZ)
#endif
#endif     

     datasize  (1:MDIM)= blkLimits  (HIGH,1:MDIM)-blkLimits  (LOW,1:MDIM)+1
     datasizeGC(1:MDIM)= blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

     
     allocate(BoxVal(nentries*product(datasize(1:NDIM)))) !nentries * total number of grid cells per block 
     
     !!-----------------------------------------------------------------------
     !!         COMPUTE CELL FACE AREAS
     !!-----------------------------------------------------------------------     

     allocate(faceAreas(NDIM, blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                              blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
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
     
     ii = 1
     BoxVal = 0.0     
     
     do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
              
              condimh = 0.
              condiph = 0.
              
              !! i-1,j,k
              if ((i /= blkLimits(LOW, IAXIS)) .or. (faces(1,IAXIS) == NOT_BOUNDARY)) then                         
                 call AssoMed(facBptr1,facBptrX,iFactorB) ! facBptr1(:,:,:) => facBptrX(iFactorB,:,:,:)
                 condimh = facBptr1(i,j,k)*faceAreas(IAXIS,i,j,k)
                 BoxVal(ii+1) = -condimh*theta*dt/del(IAXIS)
              end if
              
              !! i+1,j,k
              if (i /= blkLimits(HIGH, IAXIS) .or. (faces(2,IAXIS) == NOT_BOUNDARY)) then                                                                  
                 call AssoMed(facBptr1,facBptrX,iFactorB) ! facBptr1(:,:,:) => facBptrX(iFactorB,:,:,:)
                 condiph = facBptr1(i+1,j,k)*faceAreas(IAXIS,i+1,j,k)
                 BoxVal(ii+2) =  -condiph*theta*dt/(del(IAXIS))
              end if
              
              BoxVal(ii) =  BoxVal(ii) + (theta*((Condimh/del(IAXIS))+(Condiph/del(IAXIS)))*dt)    
              
#if NDIM >= 2
              condjmh = 0.
              condjph = 0.
              
              !! i,j-1,k
              if ((j /= blkLimits(LOW, JAXIS)) .or. (faces(1,JAXIS) == NOT_BOUNDARY)) then                                    
                 call AssoMed(facBptr1,facBptrY,iFactorB) ! facBptr1(:,:,:) => facBptrY(iFactorB,:,:,:)
                 condjmh = facBptr1(i,j,k)*faceAreas(JAXIS,i,j,k)
                 BoxVal(ii+3) = -condjmh*theta*dt/(del(JAXIS))               
              end if
              
              !! i,j+1,k
              if ((j /= blkLimits(HIGH, JAXIS)) .or. (faces(2,JAXIS) == NOT_BOUNDARY)) then                      
                 call AssoMed(facBptr1,facBptrY,iFactorB) ! facBptr1(:,:,:) => facBptrY(iFactorB,:,:,:)
                 condjph = facBptr1(i,j+1,k)*faceAreas(JAXIS,i,j+1,k)
                 BoxVal(ii+4) = -condjph*theta*dt/(del(JAXIS))                                        
              end if
              
              BoxVal(ii) =  BoxVal(ii) +  theta*(Condjmh/del(JAXIS)+Condjph/del(JAXIS))*dt !! diag           
              
              
#if NDIM == 3
              condkmh = 0.
              condkph = 0.
              
              !! i,j,k-1
              if ((k /= blkLimits(LOW, KAXIS)) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) then                                    
                 call AssoMed(facBptr1,facBptrZ,iFactorB) ! facBptr1(:,:,:) => facBptrZ(iFactorB,:,:,:)
                 condkmh = facBptr1(i,j,k)*faceAreas(KAXIS,i,j,k)
                 BoxVal(ii+5) = -condkmh*theta*dt/(del(KAXIS))               
                 
              end if
              
              !! i,j,k+1
              if ((k /= blkLimits(HIGH, KAXIS)) .or. (faces(2,KAXIS) == NOT_BOUNDARY)) then                      
                 call AssoMed(facBptr1,facBptrZ,iFactorB) ! facBptr1(:,:,:) => facBptrZ(iFactorB,:,:,:)
                 condkph = facBptr1(i,j,k+1)*faceAreas(KAXIS,i,j,k+1)
                 BoxVal(ii+6) = -condkph*theta*dt/(del(KAXIS))                                        
              end if
              
              BoxVal(ii) =  BoxVal(ii) +  theta*(Condkmh/del(KAXIS)+Condkph/del(KAXIS))*dt !! diag  
              
#endif             
              
#endif               
              ii = ii + nentries
              
           end do
        end do
     end do
     
     
     
     call HYPRE_SStructMatrixSetBoxValues(gr_hypreMatA, mypart, gr_hypreLower(lb,1:NDIM), & 
          gr_hypreUpper(lb,1:NDIM), var, nentries, stencil_indices(1:nentries), BoxVal(:), ierr)
     
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
  
  call Timers_stop("gr_hypreCreateMatrixFcB") 
  
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
end subroutine gr_hypreCreateMatrixFcB
