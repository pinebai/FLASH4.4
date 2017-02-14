!!****f* source/flashUtilities/system/memoryUsage/UnitTest2/ut_sysMemUnitTest2
!!
!! NAME
!!  ut_sysMemUnitTest2
!!
!! SYNOPSIS
!!  ./ut_sysMemUnitTest2 (standalone program)
!!
!! DESCRIPTION
!!  This is a standalone unit test independent of FLASH.
!!
!!***

#include "ut_sysMem.h"

Program ut_sysMemUnitTest2
  implicit none
  include 'mpif.h'
  integer, parameter :: XSIZE=10000
  integer, parameter :: YSIZE=10000
  real, dimension(:,:), allocatable  :: matrix
  integer :: ierr
  
  call MPI_Init (ierr)

  call print_memory_stats("START")

  allocate(matrix(XSIZE,YSIZE))
  matrix = -1.0
  call print_memory_stats("Just allocated a big heap array")

  deallocate(matrix)
  call print_memory_stats("Just deallocated a big heap array")

  call MPI_Finalize (ierr)

End Program ut_sysMemUnitTest2


subroutine print_memory_stats(msg)
  use ut_sysMemData, ONLY : memsummary_t
  use ut_sysMemInterface, ONLY  : ut_sysMemSummaryStats
  implicit none
  include 'mpif.h'
  character(len=*), intent(IN) :: msg

  integer, parameter :: comm = MPI_COMM_WORLD, verbosity = 1, &
       memsummarySize = 20, memorySampler = UT_SYSMEM_AUTO
  type(memsummary_t), dimension(memsummarySize) :: memsummary
  integer :: numStats, i, rank, ierr
  character(len=12) :: minMem_string, maxMem_string, avgMem_string
  
  call MPI_Comm_rank(comm, rank, ierr)
  call ut_sysMemSummaryStats(comm, verbosity, memorySampler, &
       memsummary, numStats)

  if (rank == 0) then
     write(*,*) msg
     do i = 1, numStats
        write (minMem_string, '(f12.2)') memsummary(i) % min
        write (maxMem_string, '(f12.2)') memsummary(i) % max
        write (avgMem_string, '(f12.2)') memsummary(i) % avg

        write (*,*) trim(memsummary(i) % description) // & 
             trim(minMem_string) // ' (min)  ' // &
             trim(maxMem_string) // ' (max)  ' // &
             trim(avgMem_string) // ' (avg) '
     end do
     write(*,*) ""
  end if
end subroutine print_memory_stats
