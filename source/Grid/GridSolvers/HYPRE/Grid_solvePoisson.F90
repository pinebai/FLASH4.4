!!****if* source/Grid/GridSolvers/HYPRE/Grid_solvePoisson
!!
!! NAME
!!
!!  Grid_solvePoisson
!!
!! SYNOPSIS
!!
!!  call Grid_solvePoisson(integer(IN) :: iSoln,
!!                         integer(IN) :: iSrc, 
!!                         integer(IN) :: bcTypes(6),
!!                         real(IN)    :: bcValues(2,6),
!!                         real(INOUT) :: poisfact)
!!
!! DESCRIPTION
!!
!!   Driver routine for poisson solvers in the Grid
!!
!!
!! ARGUMENTS
!!
!!  iSoln    - the index for the solution variable (potential when used for self-gravity)
!!  iSrc     - the index of the source variable (density when used for self-gravity)
!!  bcTypes  - the boundary condition type; only the first entry is used.
!!             Only the first 2*NDIM elements are significant. They are interpreted
!!             in the order (X left, X right, Y left, Y right, Z left, Z right).
!!             Valid values are:
!!               GRID_PDE_BND_PERIODIC (1)
!!               GRID_PDE_BND_DIRICHLET (2) (homogeneous or constant Dirichlet)
!!               GRID_PDE_BND_NEUMANN (3) (homogeneous or constant Neumann)
!!               GRID_PDE_BND_ISOLATED (0)
!!           !DEV: requested bcTypes ignored, always uses GRID_PDE_BND_NEUMANN
!!  bcValues - the values to boundary conditions, currently not used (treated as 0)
!!  poisfact      - scaling factor to be used in calculation
!!
!! NOTES
!!
!!  The symbols listed above for bcTypes are declared as FORTRAN PARAMETERS in
!!  the module Grid_interfaces.  Code using this interface should refer to that
!!  module with a USE statement, like this:
!!
!!    use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
!!       GRID_PDE_BND_ISOLATED, GRID_PDE_BND_DIRICHLET, &
!!       Grid_solvePoisson
!!  
!!  Most implementations only support a limited subset of boundary condition types.
!!  Some implementations require that all significant elements of bcTypes are the same.
!!  (That is always the case for GRID_PDE_BND_ISOLATED.)
!!  Even if an implementation supports combinations of different boundary conditions
!!  on different sides of the domain, the types at left and right sides for the same
!!  axis direction will usually have to be the samme.
!!
!!  Support in some implementations provided with FLASH4:
!!
!!   GridSolvers/Multipole:                  GRID_PDE_BND_ISOLATED                   
!!   GridSolvers/Multigrid (simple):         GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET(hom.),
!!                                            GRID_PDE_BND_NEUMANN(hom.)(?)
!!                                            (same type in all directions),
!!                                            or GRID_PDE_BND_ISOLATED
!!                                           (requires Paramesh as Grid with NBlockX==NBlockY==NBlockZ==1)
!!   GridSolvers/Pfft:                       GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET(hom.), 
!!                                            GRID_PDE_BND_NEUMANN(hom.)
!!                                            in various combinations,
!!                                            depending on GridSolvers/Pfft subdirectory 
!!                                            (i.e. implementation configured in)
!!                                           (requires UG in pencil shape or Paramesh as Grid)
!!   GridSolvers/Multigrid hybrid with Pfft: GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET, 
!!                                            GRID_PDE_BND_NEUMANN in various combinations,
!!                                            or GRID_PDE_BND_ISOLATED
!!                                           (requires Paramesh as Grid)
!!   GridSolvers/HYPRE:                      currently GRID_PDE_BND_NEUMANN only
!!   
!!***

!!REORDER(4): solnVec

subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, bcValues, poisfact)

  use Grid_data,        ONLY : gr_meshMe, gr_meshcomm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_interface,     ONLY : gr_hypreCreateMatrix, gr_hypreComputeB,    &
                               gr_hypreGridStatus
  use gr_hypreLocalInterface, ONLY: gr_hypreExchangeFacB
  use Grid_interface,   ONLY : Grid_fillGuardCells, Grid_getListOfBlocks, &
                               Grid_getBlkPtr, Grid_releaseBlkPtr,        &
                               Grid_getBlkIndexLimits, Grid_getBlkData,   &
                               Grid_getBlkRefineLevel

  
  use gr_hypreData,   ONLY   : gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, gr_hypreVecB, gr_hypreRefineMIN, &
                               gr_hypreUseFloor
  
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
       GRID_PDE_BND_NEUMANN,   &
       GRID_PDE_BND_DIRICHLET



  implicit none

#include "Flash.h"
#include "constants.h"
  
  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes(6)
  real, intent(in)       :: bcValues(2,6)
  real, intent(inout)    :: poisfact !DEV: NOT intent(IN) because some implementation actually changes it? - KW  
  
  integer :: mylevel, mypart, var, ii,i,j,k, ierr
  real, allocatable :: RHSVal(:)
  integer :: datasize(MDIM)
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer :: blockID, lb
  real, allocatable :: cellVolumes(:,:,:)
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  logical :: mask(NUNK_VARS), savedUseFloor

  integer :: iFactorB
  integer :: iFactorA
  real    :: dt, theta

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  integer :: localbcTypes(6) 
!!$    character(len=32) :: matfile

  
  call Timers_start("Grid_solvePoisson")    

  
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  
  !!-----------------------------------------------------------------------
  !!     1.  Do we need to reset HYPRE grid ?, has the underlying AMR
  !!         mesh been modified ?, is this the first call to HYPRE ?
  !!         Only in AMR is this routine meaningful.
  !!-----------------------------------------------------------------------
  call gr_hypreGridStatus (blockCount, blockList)  
  
  !! Needs to have a sensible value.
  iFactorB = 1
  mypart = 0  !! part iterator 
  var    = 0  !! var iterator.
  
  call HYPRE_SStructVectorAssemble(gr_hypreVecB, ierr)  
  call HYPRE_SStructVectorGather(gr_hypreVecB, ierr)

  do lb = 1, blockCount 
     
     blockID = blockList(lb)
     
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     call Grid_getBlkPtr(blockID, solnVec)
     call Grid_getBlkRefineLevel(blockID,mylevel)

     mypart = mylevel - gr_hypreRefineMIN
     
     datasize(1:MDIM) = blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM)+1         
     
     allocate(cellVolumes(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
     
     call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR,          &
          blkLimits(LOW,:), cellVolumes,datasize)      
     
     allocate(RHSVal(product(dataSize(1:NDIM))))     
     
     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
              
              ii = (k - blkLimits(LOW,KAXIS)  + 1)                             +  &
                   (j - blkLimits(LOW,JAXIS))*dataSize(KAXIS)                  +  &
                   (i - blkLimits(LOW,IAXIS))*dataSize(KAXIS)*dataSize(JAXIS)  
              
              RHSVal(ii) = cellVolumes(i,j,k)*solnVec(iSrc,i,j,k)*4.0*3.1416*6.67428E-08
              
           end do
        end do
     end do
     
     solnVec(iFactorB,:,:,:) = 1.0
     solnVec(iSoln,:,:,:)    = 0.0
          
     
     call HYPRE_SStructVectorSetBoxValues(gr_hypreVecB, mypart, gr_hypreLower(lb, 1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), var, RHSVal(:), ierr)  
     
     
     call Grid_releaseBlkPtr(blockID, solnVec)
     
     deallocate (RHSVal)
     deallocate(cellVolumes)
     
  end do


  mask = .false.  
  mask(iSoln) = .true.  
  mask(iSrc)  = .true.
  mask(iFactorB) = .true.  

  
  call Timers_start("Grid_fillGuardCells")   
  call Grid_fillGuardCells(CENTER,ALLDIR,masksize=NUNK_VARS, mask=mask,selectBlockType=LEAF)  
  call Timers_stop("Grid_fillGuardCells")   

  
  !!-----------------------------------------------------------------------
  !!     2. Exchange iFactorB across processors. Needs to be done only in
  !!        PARAMESH / AMR, if UG the function will return without any action.
  !!-----------------------------------------------------------------------
  call gr_hypreExchangeFacB (iFactorB, blockCount, blockList)  
  
  !!-----------------------------------------------------------------------
  !!     3. Set initial guess for solver typically 
  !!        iVar at previous time step is used.
  !!-----------------------------------------------------------------------
  call gr_hypreSetIniGuess (iSoln, blockCount, blockList)      
  
  !!-----------------------------------------------------------------------
  !!     4. Compute the actual A matrix in AX=B
  !!-----------------------------------------------------------------------  


  !! iFactorA, iFactorC not used in routine.
  iFactorA = 1
  dt       = 1.0
  theta    = 1.0 
  
  localbcTypes = GRID_PDE_BND_NEUMANN
  
  call gr_hypreCreateMatrix(iSoln, iFactorB, iFactorA, localbcTypes, bcValues, dt, theta,  &
       blockCount, blockList, .TRUE.)       
  
!!$  call gr_hypreComputeB (blockCount, blockList, iVar, iFactorA, iFactorB, dt, theta, &
!!$       bcTypes, bcValues, iFactorD)  
  
  

  
!!$  matfile = 'ex12f.out'
!!$  matfile(10:10) = char(0)
!!$  call HYPRE_SStructMatrixPrint(matfile, gr_hypreMatA, 0, ierr)
!!$  pause  
  
  !!-----------------------------------------------------------------------
  !!     7. Solve AX = B
  !!-----------------------------------------------------------------------
  call gr_hypreSolve ()
  
  !!-----------------------------------------------------------------------
  !!     8. Update unk variable using X.
  !!-----------------------------------------------------------------------
  savedUseFloor = gr_hypreUseFloor
  gr_hypreUseFloor = .FALSE.    ! Do not floor solution of Poisson problem.
  call gr_hypreUpdateSoln (iSoln, blockCount, blockList)  
  gr_hypreUseFloor = savedUseFloor    ! rstore previous value of flag
  
  
  call Timers_stop("Grid_solvePoisson") 
  
  return
end subroutine Grid_solvePoisson
