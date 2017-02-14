!!****if* source/Grid/GridSolvers/HYPRE/multiScalar/coupled/Grid_dfsvCreateSystem
!!
!! NAME
!!  Grid_dfsvCreateSystem
!!
!! SYNOPSIS
!!
!!  call Grid_dfsvCreateSystem(
!!                            integer(IN), dimension(VARDESC_SIZE) :: baseVarDesc,
!!                            integer(IN)  :: ntotVars,
!!                              integer(IN)         :: iSrc,
!!                                integer(IN)  :: iFactorA,
!!                                integer(IN)  :: bcTypes(6),
!!                                real(IN)     :: bcValues(2,6),
!!                                real(IN)     :: dt,
!!                                real(IN)     :: chi,
!!                                real(IN)     :: scaleFact,
!!                                real(IN)     :: theta,
!!                                logical(IN)  :: solnIsDelta,
!!                        OPTIONAL,integer(IN) :: iFactorC,
!!                        OPTIONAL,integer(IN) :: iFactorD
!!                        OPTIONAL,integer(IN) :: pass,
!!                        OPTIONAL,integer(OUT):: maxChunkSize)
!!  DESCRIPTION 
!!
!!      Create a system of coupled scalar equations for a diffusion problem.
!!
!!      The solver advances a generalized diffusion operator of the form
!!
!!         A*(df/dt) + C*f = div(B*grad(f)) + D ,
!!
!!      where
!!         f = f(x,t) is the  Variable to be diffused (x=1D..3D position);
!!
!!         A,B,C,D are optional given scalar factors/terms that may depend
!!         on position; they are either physcially constant in time, or at
!!         least considered time-independent for the purpose of the operation
!!         implemented here (typically by computing their values from the
!!         solution state reached by the previous time step).
!!
!!      Presently it is used to do heat conduction and multigroup diffusion.
!!
!! ARGUMENTS
!!   baseVarDesc    : A variable on which the diffusion operation is performed, usually
!!                    an extra variable, used to identify the system
!!   ntotVars       : total number or variables to solve for / equations to solve
!!   iFactorA       :| Are factors in the equation with spatial variation.
!!   factorB        :| Factor C,D are optional and are generally used
!!   iFactorC       :| to represent emission/absorption in MGD.
!!   iFactorD       :| iFactorA is needed only for conduction.
!!                   | For this FcB variant of the interface, factorB is passed in
!!                   | allocated scratch arrays, i.e., the implementation will have
!!                   | to call Grid_ascGetBlkPtr to get at it.
!!   bcTypes        : Presently OUTFLOW, VACUUM is supported, DIRICHLET is untested.
!!   bcValues       : Values of iVar,factorB on boundary (DIRICHLET).                        
!!   dt             : The time step.
!!   scaleFact      : Factor by which the end solution is scaled (not used).
!!   chi            : useful for constant diffusion problems (not used).
!!   theta          : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   pass           : Ignored in unsplit solver.
!!                    pass=1 order of directional sweep X-Y-Z, 
!!                    pass=2 order of directional sweep Z-Y-X.
!!   iSrc           : Ignored.
!!   solnIsDelta    : Is the solution only a delta that the caller has to apply to the
!!                    temperature, rather than temperature itself (ignored).
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!  This implementation
!!    * supports: 
!!           1D, 2D, 3D Cartesian PARAMESH (with local refinement)
!!           1D, 2D, 3D Cartesian in UG (3D untested).
!!           1D Spherical in PARAMESH/UG
!!           1D, 2D Cylindrical in PARAMESH/UG (R-Z)
!!
!!    *  uses HYPRE library to solve AX = B
!!
!! SEE ALSO
!!
!!  Grid_dfsvBeginSystem
!!  Grid_dfsvAddToSystem
!!  Grid_dfsvCompleteSystem
!!  Diffuse_solveCoupledScalar
!!
!!***

#include "Flash.h"

subroutine Grid_dfsvCreateSystem (baseVarDesc, ntotVars, dt,theta,thetaC,thetaD, pass, maxChunkSize)
  
  use Grid_data,        ONLY : gr_meshMe, gr_meshcomm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_interface,     ONLY : gr_hypreCreateMatrixFcB, gr_hypreComputeB, &
                               gr_hypreGridStatus
  use gr_hypreLocalInterface, ONLY: gr_hypreExchangeFacBFcB
  use Grid_interface,   ONLY : Grid_fillGuardCells, Grid_getListOfBlocks, &
                               Grid_getBlkPtr, Grid_releaseBlkPtr,        &
                               Grid_getBlkIndexLimits, Grid_getBlkData,   &
                               Grid_getBlkRefineLevel
  
#ifdef FLASH_GRID_PARAMESH
  use physicaldata,     ONLY : nfluxes
#endif
  use gr_hypreData,   ONLY   : gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, gr_hypreVecB, gr_hypreRefineMIN, &
                               gr_hyprePrintSolveInfo, &
                               gr_hypreFloor, gr_hypreRelTol, &
                               gr_hypreAbsTol, gr_hypreSolverAbsTolEff, &
                               gr_hypreSolverAutoAbsTolFact
  use gr_hypreMultiData,ONLY : gr_dfsvBaseVarDesc, &
                               gr_dfsvFirstUnkVar, gr_dfsvLastUnkVar, &
                               gr_dfsvTheta,gr_dfsvThetaC,gr_dfsvThetaD, &
                               gr_dfsvDt, &
                               gr_dfsvPhase

  implicit none

#include "constants.h"
  
  integer, dimension(VARDESC_SIZE), intent(IN):: baseVarDesc
  integer, intent(IN) :: ntotVars
!!$  integer, intent(IN) :: iSrc
!!$  integer, intent(IN) :: iFactorA
  real, intent(IN)    :: dt 
!!$  real, intent(IN)    :: chi
!!$  real, intent(IN)    :: scaleFact
  real, intent(IN)    :: theta, thetaC, thetaD
!!$  logical, intent(IN) :: solnIsDelta
!!$  integer, dimension(6),  intent(IN) :: bcTypes
!!$  real   , dimension(2,6),intent(IN) :: bcValues
  integer, intent(IN), OPTIONAL :: pass
  integer, intent(OUT),OPTIONAL :: maxChunkSize
!!$  integer, intent(IN), OPTIONAL :: iFactorC
!!$  integer, intent(IN), OPTIONAL :: iFactorD   
  
  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  integer, parameter :: fpv = 2**(NDIM-1) ! fluxes per variable
  
!!! BEGIN TEMP

  integer :: iVar

  integer :: mylevel, mypart, var, ii,i,j,k, ierr
  character(len=32) :: matfile

!!! END TEMP

  
  call Timers_start("Grid_dfsvCreateSystem")   
  
!!$  call Timers_start("Diffusion Barrier")
!!$  call MPI_Barrier (gr_meshcomm, ierr)
!!$  call Timers_stop ("Diffusion Barrier")
  
  gr_dfsvPhase = 0

  gr_dfsvBaseVarDesc = baseVarDesc
  gr_dfsvFirstUnkVar = -1
  gr_dfsvLastUnkVar  = -1
  gr_dfsvTheta  = theta
  gr_dfsvThetaC = thetaC
  gr_dfsvThetaD = thetaD

  gr_dfsvDt     = dt

  if (present(maxChunkSize)) then
#ifdef FLASH_GRID_UG
     maxChunkSize = -1
#endif
#ifdef FLASH_GRID_PARAMESH
     maxChunkSize = nfluxes/fpv
#endif
  end if

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  
  !!-----------------------------------------------------------------------
  !!     1.  Do we need to reset HYPRE grid ?, has the underlying AMR
  !!         mesh been modified ?, is this the first call to HYPRE ?
  !!         Calls gr_hypreSetupGrid if necessary.
  !!-----------------------------------------------------------------------
  call gr_hypreGridStatus (blockCount, blockList, nvars=ntotVars)
  

  
  call Timers_stop("Grid_dfsvCreateSystem") 
  
  return
  
end subroutine Grid_dfsvCreateSystem
