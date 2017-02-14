!!****if* source/IO/IOMain/io_getVarExtrema
!!
!! NAME
!!
!!  io_getVarExtrema
!!
!!
!! SYNOPSIS
!!
!!  io_getVarExtrema(integer(in) :: nvars,
!!                   real(out)   :: globalVarMin(nvars),
!!                   real(out)   :: globalVarMax(nvars),
!!                   integer(in) :: gridDataStruct)
!!          
!!
!!
!!
!! DESCRIPTION
!!
!!  This function gets the maximum and minimum values for each variable 
!!  stored in the kind of array that is indicated by gridDataStruct.
!!
!! ARGUMENTS
!! 
!!  nvars - the number of mesh variables in gridDataStruct
!!  globalVarMin - array holding max value for each variable in the data structure
!!  globalVarMax - array holding min value for each variable in the data structure
!!  gridDataStruct - one of CENTER (for UNK), SCRATCH, or FACE{X,Y,Z}
!!
!! NOTES
!!
!!***


!!REORDER(4): solnData

subroutine io_getVarExtrema(nvars, globalVarMin, globalVarMax, gridDataStruct)
      use Grid_interface, ONLY :  Grid_getListOfBlocks, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkIndexLimits
      use Driver_interface, ONLY : Driver_abortFlash
      use IO_data, only: io_unkToGlobal, io_globalComm
      implicit none
      
#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
#ifdef Grid_releaseBlkPtr
! doing Grid_releaseBlkPtr macro expansion by hand because otherwise
! it generates too long of a line for the fortran compiler, see: drift
#undef Grid_releaseBlkPtr
#endif

      integer, intent(in) :: nvars, gridDataStruct
      real, DIMENSION(nvars), INTENT(out) :: globalVarMin, globalVarMax

      real, allocatable :: varMin(:), varMax(:)
      real, dimension(:,:,:,:), pointer :: solnData

      integer :: i, j, k, n, numLeafBlocks, lb
      integer :: listOfBlocks(MAXBLOCKS)

      integer :: blockID, ierr
      integer, dimension(2,MDIM) :: blkLmts, blkLmtsGC

      if(nvars > 0) then
            allocate(varMin(nvars))
            allocate(varMax(nvars))

            ! start by finding the extrema locally on each processor
            varMin(:) = huge(varMin)
            varMax(:) = -huge(varMax)

            call Grid_getListOfBlocks(LEAF, listOfBlocks, numLeafBlocks)

            do lb = 1, numLeafBlocks
                  call Grid_getBlkPtr(listOfBlocks(lb), solnData, gridDataStruct)
                  call Grid_getBlkIndexLimits(listOfBlocks(lb), blkLmts, blkLmtsGC, gridDataStruct)

                  select case(gridDataStruct)
                  case(CENTER)
                     do k = blkLmts(1,KAXIS), blkLmts(2,KAXIS)
                              do j = blkLmts(1,JAXIS), blkLmts(2,JAXIS)
                                    do i = blkLmts(1,IAXIS), blkLmts(2,IAXIS)
                                          do n = UNK_VARS_BEGIN, UNK_VARS_END
                                                if(io_unkToGlobal(n) > 0) then
                                                varMin(io_unkToGlobal(n)) = min(varMin(io_unkToGlobal(n)), solnData(n,i,j,k))
                                                varMax(io_unkToGlobal(n)) = max(varMax(io_unkToGlobal(n)), solnData(n,i,j,k))
                                       end if
                                          enddo
                                    enddo
                              enddo
                        enddo
                  case(SCRATCH, FACEX, FACEY, FACEZ)
                        do k = blkLmts(1,KAXIS), blkLmts(2,KAXIS)
                              do j = blkLmts(1,JAXIS), blkLmts(2,JAXIS)
                                    do i = blkLmts(1,IAXIS), blkLmts(2,IAXIS)
                                          do n = 1, nvars
                                                !if(.not. io_unkActive(n)) cycle
                                                varMin(n) = min(varMin(n), solnData(n,i,j,k))
                                                varMax(n) = max(varMax(n), solnData(n,i,j,k))
                                          enddo
                                    enddo
                              enddo
                        enddo
                  case default
                        call Driver_abortFlash("io_getVarExtrema: dataStruct not implemented")
                  end select
                  
                  ! doing Grid_releaseBlkPtr expansion by hand, see: drift
                  call Driver_driftSetSrcLoc(__FILE__,__LINE__)
                  call Grid_releaseBlkPtr(listOfBlocks(lb), solnData, gridDataStruct)
                  
            enddo

            ! now do a global minimization or maximization
            call MPI_AllReduce(varMin, globalVarMin, nvars, FLASH_REAL, &
                  MPI_MIN, io_globalComm, ierr)

            call MPI_AllReduce(varMax, globalVarMax, nvars, FLASH_REAL, &
                  MPI_MAX, io_globalComm, ierr)
            
            deallocate(varMin)
            deallocate(varMax)
      endif
      return
end subroutine io_getVarExtrema

