!!****if* source/Simulation/SimulationMain/unitTest/PFFT_PoissonFD/Grid_unitTest
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
!!  This unit test exercises PFFT unit.
!!
!! ARGUMENTS
!!
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!!
!! NOTES
!!
!!***

!!REORDER(4): solnData

subroutine Grid_unitTest(fileUnit,perfect)

  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, &
       Grid_solvePoisson,Grid_getBlkIndexLimits, &
       Grid_getBlkPtr,Grid_releaseBlkPtr, Grid_getListOfBlocks, &
       Grid_getDeltas, Grid_fillGuardCells, Grid_getBlkRefineLevel
  use gr_interface ,ONLY : gr_findMean
  use Grid_data, ONLY : gr_meshMe, gr_meshComm
#include "Flash.h"
#include "constants.h"


  implicit none
#include "Flash_mpi.h"

  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors

  real, pointer, dimension(:,:,:,:) :: solnData

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(2*MDIM) :: bcTypes
  real,dimension(2,2*MDIM) :: bcValues

  logical :: gcMask(NUNK_VARS), evaluateData

  real :: poisfact
  integer,dimension(MAXBLOCKS) :: blkList
  integer :: blockID,blkCount,lb,i,j,k
  real:: del(MDIM)
  real meanASOL,meanPFFT
  integer nx,ny,nz
  real, allocatable, dimension(:,:,:) :: fsrc, usol
  integer blkpoints, blkpointsaux,blkCountaux
  real L2_err, L2_erraux, Linf_err, Linf_erraux, Tvol, Tvolaux, vcell
  integer :: refinelevel
  real, parameter :: tol = 1.e-6
  real, parameter :: tolInf = 4.1e-2
  real, parameter :: tol2  = 2.24e-2

  integer TA(2),count_rate,ierr
  real :: ET

  ! -------------------------------------------------------------------
  bcTypes(:)=GRID_PDE_BND_PERIODIC  !This is always periodic for this problem.
  bcValues(:,:)=0.

  call mpi_barrier(gr_meshComm,ierr)
  if (gr_meshMe .eq. 0) CALL SYSTEM_CLOCK(TA(1),count_rate)  
  poisfact=1.
  call Grid_solvePoisson(PFFT_VAR,DENS_VAR,bcTypes,bcValues,poisfact)
  call mpi_barrier(gr_meshComm,ierr)
  if (gr_meshMe .eq. 0) then
     CALL SYSTEM_CLOCK(TA(2),count_rate)
     ET=REAL(TA(2)-TA(1))/count_rate
     write(*,*) ' ' 
     write(*,*) '3 PERIODIC Poisson Solver time = ',ET,' sec.'
  endif


  ! Get Block list
  call Grid_getListOfBlocks(LEAF,blkList,blkCount)

  ! Check error in the solution:
  L2_err = 0.
  blkpoints = 0
  Linf_err = 0.
  Tvol = 0.
  do lb = 1,blkCount
     blockID = blkList(lb)

     call Grid_getBlkRefineLevel(blockID,refineLevel)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     call Grid_getDeltas(blockID,del)

     select case (NDIM)
     case(1)
        vcell = del(IAXIS)
     case(2)
        vcell = del(IAXIS)*del(JAXIS)
     case(3)
        vcell = del(IAXIS)*del(JAXIS)*del(KAXIS)
     end select

     blkpoints = blkpoints + &
          (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1) * &
          (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1) * &
          (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)

     Tvol = Tvol + vcell*real( (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1) * &
          (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1) * &
          (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1))

     solnData(DIFF_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),           &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   =       &
          abs(  solnData(PFFT_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   - &
          solnData(ASOL_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   )

     ! L2 norm of error:
     L2_err = L2_err + sum( vcell*solnData(DIFF_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),           &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))**2.)

     ! Linf norm of error:
     Linf_err = max(Linf_err,maxval(solnData(DIFF_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) ))

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  enddo


  ! Sum processors Volumes
  Tvolaux = Tvol
  call MPI_Allreduce(Tvolaux, Tvol, 1, FLASH_REAL,&
       MPI_SUM, gr_meshComm, ierr)

  ! Sum processors points
  blkpointsaux = blkpoints
  call MPI_Allreduce(blkpointsaux, blkpoints, 1, FLASH_INTEGER,&
       MPI_SUM, gr_meshComm, ierr)

  ! Sum processors L2 norm of error squared
  L2_erraux = L2_err
  call MPI_Allreduce(L2_erraux, L2_err, 1, FLASH_REAL,&
       MPI_SUM, gr_meshComm, ierr)


  ! Compute L2 norm for whole domain
  L2_err = sqrt(L2_err/Tvol)

  ! Sum processors Linf norm of error
  Linf_erraux = Linf_err
  call MPI_Allreduce(Linf_erraux, Linf_err, 1, FLASH_REAL,&
       MPI_MAX, gr_meshComm, ierr)


  ! Fill GuardCells:
  gcMask = .TRUE.                         
  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS,mask=gcMask)          

  ! Export to Tecplot:
  !call outtotecplot(gr_meshMe,0.0,1.,1,0,0.0,blkList,blkCount,0)


  call gr_findMean(ASOL_VAR,2,.false.,meanASOL)
  call gr_findMean(PFFT_VAR,2,.false.,meanPFFT)

  !Unit test gives a mean analytical solution of zero.  Ensure the absolute
  !value of the mean numerical solution is between 0 and the tolerance value.
  perfect = abs(meanPFFT) < MAX(TINY(1.), tol)

  if (perfect) perfect = (abs(Linf_err) .LE. tolInf)
  if (perfect) perfect = (L2_err .LE. tol2)

  ! Sum processors points
  call MPI_Allreduce(blkCount, blkCountaux, 1, FLASH_INTEGER,&
       MPI_SUM, gr_meshComm, ierr)


  if (gr_meshMe .eq. 0) then
     write(*,*) ' ' 
     write(*,'(A,2g16.8)') ' Mean Analytical, Numerical Sol=',meanASOL,meanPFFT
     write(*,'(A,1g16.8)') " ||Phi - PhiAnalytical||inf =" ,Linf_err
     write(*,'(A,1g16.8)') " ||Phi - PhiAnalytical||2   =" ,L2_err
     write(*,*) " Total Volume =",Tvol
     write(*,*) " Total Number of Leaf Blocks=", blkCountaux
     write(*,*) ' ' 
!!$  print*,"Processor", gr_meshMe, "the Scalar result is Linf err=", err2
!!$  print*,"Processor", gr_meshMe, "Pfft against Scalar result is Linf err=", err3
  endif
  return

end subroutine Grid_unitTest
