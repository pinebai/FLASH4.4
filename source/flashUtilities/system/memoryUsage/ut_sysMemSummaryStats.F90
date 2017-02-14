!!****if* source/flashUtilities/system/memoryUsage/ut_sysMemSummaryStats
!!
!! NAME
!!  ut_sysMemSummaryStats
!!
!! SYNOPSIS
!!
!!  ut_sysMemSummaryStats(integer(IN)                :: comm,
!!                        integer(IN)                :: verbosity,
!!                        integer(IN)                :: memorySampler,
!!                        type(memsummary_t)(:)(OUT) :: memsummary
!!                        integer(OUT)               :: numStats)
!!  
!! DESCRIPTION
!!  
!!  Returns a memory usage summary of all MPI ranks in a given MPI
!!  communicator.  The summary includes the minimum, maximum and
!!  mean memory usage for each memory statistic provided by a
!!  given memory sampler.
!!
!! ARGUMENTS
!!
!!  comm : The MPI communicator
!!
!!  verbosity : The verbosity of the memory usage information.  There are
!!              currently two levels
!!               0 = regular
!!               1 = more verbose
!!
!!  memorySampler : The method for sampling memory.  The following are valid
!!                   UT_SYSMEM_PROC = use /proc/self/stat
!!                   UT_SYSMEM_RUSAGE = use rusage
!!                   UT_SYSMEM_MALLINFO = use mallinfo (f2003)
!!                   UT_SYSMEM_BGKERNEL = use Kernel_GetMemorySize (f2003)
!!                   UT_SYSMEM_AUTO = use the default sampler for this machine
!!                   UT_SYSMEM_ALL = use all samplers supported by this machine
!!                  (f2003) indicates that the sampler is only supported
!!                  by the memory usage code in directory f2003.  This is
!!                  only included in a FLASH application if you setup with
!!                  useFortran2003=True.
!!
!!  memsummary : An array containing the returned memory summary.
!!
!!  numStats : The number of memory statistics memsummary.
!!
!!
!! NOTES
!!
!!  This subroutine performs collective MPI communication and thus
!!  must be called by all MPI ranks in the passed communicator.
!!
!!***

#include "constants.h"

subroutine ut_sysMemSummaryStats(comm, verbosity, memorySampler, &
     memsummary, numStats)
  use ut_sysMemInterface, ONLY : ut_sysMemStats
  use ut_sysMemData, ONLY : memsummary_t, meminfo_t
  implicit none
  include "Flash_mpi.h"
  integer, intent(IN) :: comm, verbosity, memorySampler
  type(memsummary_t), dimension(:), intent(OUT) :: memsummary
  integer, intent(OUT) :: numStats

  type(meminfo_t), dimension(size(memsummary,1)) :: meminfo
  real, dimension(2,size(memsummary,1)) :: minMem, maxMem, memAndRank
  real, dimension(size(memsummary,1)) :: totMem, mem
  integer :: commSize, commRank, ierr, i

  call MPI_Comm_rank(comm, commRank, ierr)
  call MPI_Comm_size(comm, commSize, ierr)

  call ut_sysMemStats(verbosity, memorySampler, meminfo, numStats)
  do i = 1, numStats
     memsummary(i) % description = meminfo(i) % description
     memAndRank(1,i) = meminfo(i) % measurement
     memAndRank(2,i) = real(commRank)
     mem(i) = memAndRank(1,i)
  end do

  call MPI_AllReduce(memAndRank, minMem, numStats, FLASH_2REAL, &
       MPI_MINLOC, comm, ierr)
  call MPI_AllReduce(memAndRank, maxMem, numStats, FLASH_2REAL, &
       MPI_MAXLOC, comm, ierr)
  call MPI_AllReduce(mem, totMem, numStats, FLASH_REAL, &
       MPI_SUM, comm, ierr)

  do i = 1, numStats
     memsummary(i) % measurement = mem(i)
     memsummary(i) % min = minMem(1,i)
     memsummary(i) % max = maxMem(1,i)
     memsummary(i) % avg = totMem(i) / commSize
     memsummary(i) % min_rank = int(minMem(2,i))
     memsummary(i) % max_rank = int(maxMem(2,i))
  end do
end subroutine ut_sysMemSummaryStats
