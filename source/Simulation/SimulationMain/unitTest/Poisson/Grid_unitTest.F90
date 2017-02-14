!!****if* source/Simulation/SimulationMain/unitTest/Poisson/Grid_unitTest
!!
!! NAME
!!
!!  Grid_unitTest
!!
!! SYNOPSIS
!!
!!  call Grid_unitTest(integer(in):: fileUnit,
!!                     logical(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  This unit test exercises a Poisson solver through the Grid_solvePoisson
!!  interface.
!!
!! ARGUMENTS
!!
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!!
!! NOTES
!!
!!  Additional diagnostic info can be written to the 'fileUnit' if desired,
!!  independent of the test outcome that is expressed in the value
!!  returned in 'perfect'.
!!***

!!REORDER(4): solnData

subroutine Grid_unitTest(fileUnit,perfect)

#include "Flash.h"

  use Grid_data, ONLY: gr_globalMe, gr_globalComm
  use gr_solversTestData, ONLY: gr_testTolL2, gr_testTolLinf

  use Driver_interface, ONLY : Driver_abortFlash

  use Grid_interface,    ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
                                GRID_PDE_BND_DIRICHLET, &
                                Grid_getListOfBlocks, &
                                Grid_updateRefinement, &
                                Grid_getBlkPtr, &
                                Grid_getBlkIndexLimits, &
                                Grid_solvePoisson, &
                                Grid_releaseBlkPtr, &
                                Grid_getBlkRefineLevel, &
                                Grid_getDeltas

  use Grid_interface,    ONLY : Grid_computeVarMean

  use RuntimeParameters_interface, ONLY: RuntimeParameters_get, &
                                         RuntimeParameters_mapStrToInt

  implicit none

#include "constants.h"
 include "Flash_mpi.h"

  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors


  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: solnData

  integer, dimension(6) :: bc_types, gr_PoissonBcTypes
  real, dimension(2,6)  :: bc_values
  real :: poisfact

  integer :: blkpoints, blkpointsaux
  real :: Tvol, Tvolaux, del(MDIM), vcell
  real :: L2_err, L2_erraux, Linf_err, Linf_erraux

  integer :: refineLevel,blockID,lb,ierr

  integer :: TA(2),count_rate
  real ::ET

  real :: meanPHI,meanPHIaux,meanANL,meanANLaux

  integer :: eachBoundary
  character(len=MAX_STRING_LENGTH), dimension(2*MDIM) :: bcTypeStr



  meanPHI = 0.0
  meanPHIaux =0.0
  meanANL = 0.0
  meanANLaux = 0.0


  !Initialise to -1 to help us spot errors.
  bc_types = -1

  ! Retrieve Poisson solution Boundary Conditions for each face:
  call RuntimeParameters_get("xl_boundary_type", bcTypeStr(1))
  call RuntimeParameters_get("xr_boundary_type", bcTypeStr(2))
  if (NDIM >= 2) then
     call RuntimeParameters_get("yl_boundary_type", bcTypeStr(3))
     call RuntimeParameters_get("yr_boundary_type", bcTypeStr(4))
  endif
  if (NDIM == 3) then
     call RuntimeParameters_get("zl_boundary_type", bcTypeStr(5))
     call RuntimeParameters_get("zr_boundary_type", bcTypeStr(6))
  endif

  ! Boundary Conditions:
  do eachBoundary = 1, 2*NDIM

     call RuntimeParameters_mapStrToInt(bcTypeStr(eachBoundary), gr_PoissonBcTypes(eachBoundary))

     select case(gr_PoissonBcTypes(eachBoundary))
     case (PERIODIC)
        bc_types(eachBoundary)  = GRID_PDE_BND_PERIODIC
        if (gr_globalMe .eq. 0) write(*,*) eachBoundary,'GRID_PDE_BND_PERIODIC'
     case (OUTFLOW)
        bc_types(eachBoundary)  = GRID_PDE_BND_NEUMANN
        if (gr_globalMe .eq. 0) write(*,*) eachBoundary,'GRID_PDE_BND_NEUMANN'
     case (DIRICHLET)
        bc_types(eachBoundary)  = GRID_PDE_BND_DIRICHLET
        if (gr_globalMe .eq. 0) write(*,*) eachBoundary,'GRID_PDE_BND_DIRICHLET'
     case default
           write(*,*) 'In DriverEvolveFlash: Poisson Problem Boundary Condition not supported'
           call Driver_abortFlash('BCs unsupported')
     end select
  enddo


  bc_values(1,:) = 0.
  bc_values(2,:) = -1.


  ! Call Poisson Solver
  if (gr_globalMe .eq. 0) write(*,*) 'Into Grid Solve Poisson ..'


  call mpi_barrier(gr_globalComm,ierr)
  if (gr_globalMe .eq. 0) CALL SYSTEM_CLOCK(TA(1),count_rate)


  poisfact = 1.0
  call Grid_solvePoisson(VPHI_VAR, VSRC_VAR, bc_types, bc_values, poisfact)


  call mpi_barrier(gr_globalComm,ierr)
  if (gr_globalMe .eq. 0) then
      CALL SYSTEM_CLOCK(TA(2),count_rate)
      ET=REAL(TA(2)-TA(1))/count_rate
      write(*,*) 'Poisson Solver time =',ET
   endif


  ! Mean Phi
  call Grid_computeVarMean(VPHI_VAR,meanPHI)
  call Grid_computeVarMean(VANL_VAR,meanANL)

  !! Errors.
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  ! Check error in the solution:
  L2_err    = 0.
  blkpoints = 0
  Linf_err  = 0.
  Tvol      = 0.
  do lb = 1,blockCount

     blockID = blockList(lb)

     call Grid_getBlkRefineLevel(blockID,refineLevel)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     call Grid_getDeltas(blockID,del)

     select case (NDIM)
     case(1)
        vcell = del(IAXIS)
     case(2)
        vcell = del(IAXIS)*del(JAXIS)
     case(3)
        vcell = del(IAXIS)*del(JAXIS)*del(KAXIS)
     end select

     solnData(VERR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),           &
                       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           &
                       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   =       &
     abs(  solnData(VPHI_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   - &
           meanPHI                                                         - &
           solnData(VANL_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   )

     blkpoints = blkpoints + &
                 (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1) * &
                 (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1) * &
                 (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)

     Tvol = Tvol + &
     vcell*real( (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1) * &
                 (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1) * &
                 (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1) )

     ! L2 norm of error:
     L2_err = L2_err + &
     vcell*sum(solnData(VERR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
                                 blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                                 blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))**2)

     ! Linf norm of error:
     Linf_err = max(Linf_err,maxval(solnData(VERR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
                                              blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           &
                                              blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) ))

     ! Release Pointer
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  enddo

 ! Sum processors Volumes
 Tvolaux = Tvol
 call MPI_Allreduce(Tvolaux, Tvol, 1, FLASH_REAL,&
      MPI_SUM, gr_globalComm, ierr)

  ! Sum processors cells
  blkpointsaux = blkpoints
  call MPI_Allreduce(blkpointsaux, blkpoints, 1, FLASH_INTEGER,&
                     MPI_SUM, gr_globalComm, ierr)

  ! Sum processors L2 norm of error squared
  L2_erraux = L2_err
  call MPI_Allreduce(L2_erraux, L2_err, 1, FLASH_REAL,&
                     MPI_SUM, gr_globalComm, ierr)

  ! Compute L2 norm for whole domain
  L2_err = sqrt(1./Tvol * L2_err)

  ! Sum processors Linf norm of error
  Linf_erraux = Linf_err
  call MPI_Allreduce(Linf_erraux, Linf_err, 1, FLASH_REAL,&
                     MPI_MAX, gr_globalComm, ierr)


               perfect = (abs(Linf_err) .LE. gr_testTolLinf)
  if (perfect) perfect = (L2_err .LE. gr_testTolL2)

  if ( gr_globalMe .eq. 0) then
     write(*,*) 'Mean Anl, Num=',meanANL,meanPHI
     write(*,*) 'L2 error = ',L2_err
     write(*,*) 'Linf error = ',Linf_err
     write(*,*) '############################################'
  endif


  return

end subroutine Grid_unitTest
