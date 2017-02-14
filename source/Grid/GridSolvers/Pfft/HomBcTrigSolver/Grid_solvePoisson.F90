!!****if* source/Grid/GridSolvers/Pfft/HomBcTrigSolver/Grid_solvePoisson
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
!!   Poisson solver routine.  This implementation provides the
!!   FFT based method for periodic, homogeneous Dirichlet, and
!!   homogeneous Neumann problems on a uniform grid.
!!   Mixed boundary problems, in which different boundary types
!!   apply at different boundary faces, are also supported.
!!   
!!   Isolated problems are not supported.
!!
!!
!! ARGUMENTS
!!
!!  iSoln -  index to variable containing potential
!!  iSrc - index to variable containing density
!!  bcTypes - boundary types along various faces,
!!          valid values are: (although only some are implemented)
!!          GRID_PDE_BND_PERIODIC (1) (supported)
!!          GRID_PDE_BND_DIRICHLET (2) (homogeneous Dirichlet supported)
!!          GRID_PDE_BND_NEUMANN (3) (homogeneous Neumann supported)
!!          GRID_PDE_BND_ISOLATED (0) (not supported in this implementation)
!!  bcValues - the values to boundary conditions, currently not used
!!  poisfact -  factor to be used in calculation
!!
!! NOTES
!!
!!  The symbols listed above for bcTypes are declared as FORTRAN PARAMETERS in
!!  the module Grid_interfaces.  Code using this interface should refer to that
!!  module with a USE statement, like this:
!!
!!    use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
!!       GRID_PDE_BND_DIRICHLET, &
!!       Grid_solvePoisson
!!
!!***
!#define DEBUG_MAPPING
!#define PRINT_DEBUG_MAPPING

!!REORDER(4): solnData

subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, bcValues, poisfact)

  use gr_pfftData, ONLY : pfft_inLen,pfft_outLen,pfft_setupOnce,pfft_usableProc,pfft_localLimits

  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
       Grid_pfftMapToInput,Grid_getGlobalIndexLimits,&
       Grid_getBlkIndexLimits,Grid_pfftInit, Grid_pfft,Grid_pfftMapFromOutput, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_pfftFinalize
  use gr_pfftInterface, ONLY : gr_pfftDerivs, gr_pfftSpecifyTransform
  use gr_interface, ONLY:  gr_findMean
  use Grid_data, ONLY : gr_meshNumProcs, gr_meshMe, gr_meshComm
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Pfft.h"  

  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes(2*MDIM)
  real, intent(in)       :: bcValues(2,2*MDIM)
  real, intent(inout)    :: poisfact !DEV: NOT intent(IN) because some implementation actually changes it? - KW

  !--------------------------------------------------------------------------
  real, allocatable, dimension(:) :: inArray,tranArray,outArray
  real :: meanDensity, normScale
  integer, dimension(MDIM) :: localSize,globalSize,transformType
  integer, dimension(0:MDIM) :: baseDatType
  integer :: inSize,tranSize
  logical :: needMap

  !=========================================================================================


  if(.not.pfft_setupOnce) then
     call Grid_getGlobalIndexLimits(globalSize)

     call gr_pfftSpecifyTransform(transformType, baseDatType, bcTypes)
     needMap=.true.

     call Grid_pfftInit(NDIM,needMap,globalSize,&
          localSize,transformType, baseDatType)
  end if

  !Important.  Tests that this processor should be doing work
  if(.not.pfft_usableProc) return   

  !This is a temporary test used to check the mapping to and 
  !from PFFT grid independently of the Tranpose & FFT.
#ifdef DEBUG_MAPPING
  call CheckMapping(iSrc)
#endif

  inSize=pfft_inLen(IAXIS)*pfft_inLen(JAXIS)*pfft_inLen(KAXIS)
  tranSize=2*pfft_outLen(IAXIS)*pfft_outLen(JAXIS)*pfft_outLen(KAXIS)

  allocate(inArray(inSize+2))
  allocate(outArray(inSize+2))
  allocate(tranArray(tranSize))


  !! Here's the real work of the fft
  ! Converts to uniform mesh (on output, inArray contains uniformly mapped density)
  call Grid_pfftMapToInput(iSrc,inArray) 

  ! Figure out the mean of the density and subtract
  call gr_findMean(iSrc,2,.false.,meanDensity)
  !print *, "meanDensity = ",meanDensity
  if (ALL(bcTypes(1:2*NDIM)==GRID_PDE_BND_PERIODIC .OR. &
       bcTypes(1:2*NDIM)==GRID_PDE_BND_NEUMANN)) then
     inArray(1:inSize) = inArray(1:inSize) - meanDensity !!was 1.5e+7
  end if

  ! Forward transform of density 
  call Grid_pfft(PFFT_FORWARD,inArray,tranArray)
!!$  print*,'the local limits',pfft_localLimits

  ! Calculates the transform of iSoln = GPOT
  !  which is the transform of the delSquared(u) = rho in Poisson equation
  call gr_pfftDerivs(tranArray)

  ! Inverse transform of GPOT
  call Grid_pfft(PFFT_INVERSE,tranArray,outArray)

  ! Now multiply by the poisson factor
   outArray(1:inSize) = outArray(1:inSize)*poisfact

   ! We do not multiply with other scaling factors here, the right scaling
   ! should already be applied in Grid_pfft and in gr_pfftDerivs.  And if
   ! other factors are needed, the caller could pass them in as part
   ! of poisfact.

  ! We do not have to add back the subtracted mean, it is
  ! not necessary because we never affected the actual density array.

  ! Map back to the non-uniform mesh
  call Grid_pfftMapFromOutput(iSoln,outArray)

  


  deallocate(inArray)
  deallocate(tranArray)
  deallocate(outArray)
  if(.not.pfft_setupOnce) then
     call Grid_pfftFinalize()
  end if

  return
end subroutine Grid_solvePoisson


!!-------------------------------------------------------------------------------------------------------------
!This is used to test only the mapping to and from PFFT grid for 
!UG and fixed refinement PARAMESH grids only.
!Each leaf block on a processor is filled with a globally unique 
!block identifier. We then perform the mapping to and from PFFT grid, and
!check if we are left with the original global block identifier on each block.
subroutine CheckMapping(iSrc)

  use Grid_interface, ONLY : Grid_pfftMapToInput,Grid_getGlobalIndexLimits,&
       Grid_getBlkIndexLimits,Grid_pfftMapFromOutput, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getListOfBlocks, &
       Grid_getBlkCornerID
  use gr_pfftData, ONLY : pfft_inLen, pfft_comm
  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs, gr_meshComm
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Flash_mpi.h"

  integer, intent(in)    :: iSrc

  real, dimension(:,:,:,:), pointer :: solnData
  real, allocatable, dimension(:) :: inArray
  character (len=*), PARAMETER :: baseName3="SolutionArray",  baseName2="IntermediateArray", &
       baseName="InitialArray"
  character (len=200) :: fileName, fileName2, fileName3  !** NOTE 200 characters!!!
  character (len=200) :: uniqueName, uniqueName2, uniqueName3 !** NOTE 200 characters!!!
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(MAXBLOCKS) :: listOfBlocks
  integer,dimension(IAXIS:KAXIS) :: stride,cornerID,startPos,endPos,blkSize 
  integer, dimension(:), allocatable :: globalLeafCount, globalLeafCounttmp
  integer :: inSize, ierr, i, expectedValue, error
  integer :: localNumbLeafBlocks, startBlockID, blockID


  if (gr_meshComm /= pfft_comm(IAXIS)) then
     call Driver_abortFlash("Use fewer processors in this test!")
  end if


  !Define an arbitrary global block identifier for each leaf block on a processor.
  !---------------------------------------------------------------------
  call Grid_getListOfBlocks(LEAF,listofBlocks,localNumbLeafBlocks)
  allocate(globalLeafCount(gr_meshNumProcs), globalLeafCounttmp(gr_meshNumProcs))
  globalLeafCount = 0; globalLeafCounttmp = 0;
  globalLeafCounttmp(gr_meshMe+1) = localNumbLeafBlocks
  call MPI_ALLREDUCE(globalLeafCounttmp,globalLeafCount,gr_meshNumProcs,FLASH_INTEGER,MPI_SUM,&
       gr_meshComm,ierr)  !MPI_IN_PLACE is an MPI-2 construct, hence we use 2 arrays.
  if (gr_meshMe == 0) then
     startBlockID = 0
  else
     startBlockID = sum(globalLeafCount(1:gr_meshMe))
  end if


  !Loop over every single leaf block, and label the internal grid 
  !points with the global block ID.
  !---------------------------------------------------------------------
  do i = 1, localNumbLeafBlocks
     blockID = listofBlocks(i)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     expectedValue = startBlockID + i - 1
     solnData(iSrc,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = real(expectedValue)

#ifdef PRINT_DEBUG_MAPPING
     write(uniqueName, fmt='(i3.3,''.dat'')') startBlockID + i - 1
     fileName = trim(baseName) // uniqueName
     open(unit=20, file=filename)
     write(20,*) int(solnData(iSrc,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
     close(20)
#endif

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  end do


  !Allocate PFFT array, and fill with garbage.
  !---------------------------------------------------------------------
  inSize=pfft_inLen(IAXIS)*pfft_inLen(JAXIS)*pfft_inLen(KAXIS) 
  allocate(inArray(inSize+2), STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("Severe error: Memory cannot be allocated!")
  end if
  inArray(:) = -8.0


  !Map from grid to PFFT array.  Print out PFFT array (this should 
  !be a mix of many different integer block IDs).
  !---------------------------------------------------------------------
  call Grid_pfftMapToInput(iSrc,inArray)
#ifdef PRINT_DEBUG_MAPPING
  write(uniqueName2, fmt='(i3.3,''.dat'')') gr_meshMe
  fileName2 = trim(baseName2) // uniqueName2
  open(unit=20, file=filename2)
  write(20,*) int(inArray(1:inSize))
  close(20)
#endif


  !Now all the data is in PFFT array. Fill the leaf blocks with garbage.
  !---------------------------------------------------------------------
  do i = 1, localNumbLeafBlocks
     blockID = listofBlocks(i)
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     solnData(iSrc,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = -9.0
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  end do


  !Map from PFFT array back to leaf blocks.
  !---------------------------------------------------------------------
  call Grid_pfftMapFromOutput(iSrc,inArray)


  !If the map back and forth went OK, the leaf blocks will 
  !contain just their original global block identifier.
  !---------------------------------------------------------------------
  do i = 1, localNumbLeafBlocks
     blockID = listofBlocks(i)
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     expectedValue = startBlockID + i - 1

     if(all(int(solnData(iSrc,&
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))) == expectedValue)) then
        print *, "Success!!! - Processor", gr_meshMe, "predicted global ID.", expectedValue
     else
        print *, "Failure!!! - Processor", gr_meshMe, "messed up during the map."
        call Driver_abortFlash("FAILURE DURING MAP")
     end if

#ifdef PRINT_DEBUG_MAPPING
     write(uniqueName3, fmt='(i3.3,''.dat'')') expectedValue
     fileName3 = trim(baseName3) // uniqueName3
     open(unit=20, file=filename3)
     write(20,*) int(solnData(iSrc,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     close(20)
#endif

  end do

  deallocate(inArray, globalLeafCount, globalLeafCounttmp, STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("Severe error: Memory cannot be deallocated!")
  end if
  call MPI_Finalize(ierr)
  stop

end subroutine CheckMapping
