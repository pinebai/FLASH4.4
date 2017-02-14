!!****if* source/Driver/DriverMain/Driver_initParallel
!!
!! NAME
!!
!!  Driver_initParallel
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!  Initialize the parallel message-passing interface,
!!  the number of processors in a run and each processing
!!  element
!!
!!
!!  ARGUMENTS
!!    myPE : current processor
!!    numProcs : number of processors
!!  
!!
!!
!!***

!In mpif.h MPI_VERSION is an integer (and thus can't be used to conditionally
!compile code) and so we allow the user to define FLASH_MPI1 to indicate
!that they have an MPI-1 implementation.

#ifdef _OPENMP
#ifndef FLASH_MPI1
#define FLASH_MPI2_OPENMP
#endif
#endif

subroutine Driver_initParallel ()

  use Driver_data, ONLY : dr_globalMe, dr_globalNumProcs, dr_globalComm, &
       dr_mpiThreadSupport
  !$ use omp_lib
  
  implicit none             

  include "Flash_mpi.h"
  integer :: error, iprovided, errcode

#ifdef _OPENMP
#ifdef FLASH_MPI2_OPENMP
  integer, parameter :: MPI_thread_level = MPI_THREAD_SERIALIZED
#endif
#ifdef __INTEL_COMPILER
  integer(kind=kmp_size_t_kind) :: stksize
#endif
#endif
  logical :: mpiThreadSupport
  mpiThreadSupport = .false.

  !We should use MPI_Init_thread rather than MPI_Init when using multiple
  !threads so that we get a guaranteed level of thread support.

#ifdef FLASH_MPI2_OPENMP
  !We have some OpenMP parallel regions spanning MPI calls - any such
  !MPI calls are currently contained in $omp single sections and so
  !we use MPI_THREAD_SERIALIZED to give us exactly the thread support we need
  !to operate safely.  I print a warning message to the screen when your
  !MPI installation is not providing this level of thread support - it
  !is up to you whether you are happy with this risk.

  !Support Levels                     Description
  !MPI_THREAD_SINGLE     Only one thread will execute.
  !MPI_THREAD_FUNNELED   Process may be multi-threaded, but only main
  !                      thread will make MPI calls (calls are funneled to
  !                      main thread). "Default"
  !MPI_THREAD_SERIALIZED Process may be multi-threaded, any thread can
  !                      make MPI calls, but threads cannot execute MPI
  !                      calls concurrently (MPI calls are serialized).
  !MPI_THREAD_MULTIPLE   Multiple threads may call MPI, no restrictions.
   
  !The MPI standard says that "a call to MPI_INIT has the same effect as
  !a call to MPI_INIT_THREAD with a required = MPI_THREAD_SINGLE".

  call MPI_Init_thread(MPI_thread_level, iprovided, error)
  if (error /= MPI_SUCCESS) then
     print *, "Error from MPI_Init_thread"
     stop
  end if
#else
  call MPI_Init (error)
#endif


  dr_globalComm=FLASH_COMM
  call MPI_Comm_Rank (dr_globalComm, dr_globalMe, error)
  call MPI_Comm_Size (dr_globalComm, dr_globalNumProcs, error)


#ifdef _OPENMP
  if (dr_globalMe == 0) then

# ifdef FLASH_MPI2_OPENMP
     !The default thread support in Open-MPI (in the versions I have used) is
     !MPI_THREAD_SINGLE unless you configure Open-MPI with --enable-mpi-threads.

     !On Cray systems the MPI environment is limited to MPI_THREAD_SINGLE
     !by default.  This can be changed with the environmental variable
     !MPICH_MAX_THREAD_SAFETY - it has possible values of "single", "funneled",
     !"serialized" or "multiple".  To obtain MPI_THREAD_MULTIPLE thread level:
     !1) Set MPICH_MAX_THREAD_SAFETY to multiple in job submission script:
     !   export MPICH_MAX_THREAD_SAFETY="multiple"
     !2) link FLASH against a special MPI library:
     !   -lmpich_threadm.
     write(6,'(a,i3,a,i3)') " [Driver_initParallel]: "//&
          "Called MPI_Init_thread - requested level ", MPI_thread_level, &
          ", given level ", iprovided
     mpiThreadSupport = (iprovided >= MPI_thread_level);
# endif

     if (.not.mpiThreadSupport) then
        write(6,"(/ a /)") " [Driver_initParallel]: WARNING! We do not have "//&
             "a safe level of MPI thread support! (see Driver_initParallel.F90)"
        !write(6,*) "[Driver_initParalllel]: ERROR! MPI thread support too limited"
        !call MPI_Abort (dr_globalComm, errcode, error)
        !stop
     end if
  end if

  !$omp parallel
  if (dr_globalMe == 0) then
     if (omp_get_thread_num() == 0) then
        write(6,'(a,i3)') " [Driver_initParallel]: "//&
             "Number of OpenMP threads in each parallel region", &
             omp_get_num_threads()

        !Add Intel compiler specific code.  It is possible to overflow the
        !stack of the spawned OpenMP threads (e.g. WD_def 3d with block list
        !threading).  The default value for intel software stack on
        !code.uchicago.edu is 4MB (it is useful to print this information).
        !I recommend increasing this to 16MB:
        !export OMP_STACKSIZE="16M".
# ifdef __INTEL_COMPILER
        stksize = kmp_get_stacksize_s() / (1024*1024)
        write(6,'(a,i8,a)') " OpenMP thread stack size:", stksize, " MB"
# endif

        !Add Absoft compiler specific code.  The same loop iteration is
        !executed by multiple threads in parallel do loops that have 1 loop
        !iteration!  This bug happens when compiling the following test
        !problem with Absoft 64-bit Pro Fortran 11.1.4 on code.uchicago.edu.
        !
        !./setup unitTest/Multipole -auto -geometry=cartesian -3d -maxblocks=1 \
        !  +newMpole +noio threadBlockList=True -nxb=64 -nyb=64 -nzb=64
        !
        ! Set lrefine_min = lrefine_max = 1 in the flash.par.
# ifdef __ABSOFT__
        print *, ""
        print *, "WARNING!!!! Absoft compiler OpenMP bug!!!!"
        print *, "A parallel do loop with 1 loop iteration will be executed incorrectly"
        print *, ""
# endif

     end if
  end if
# ifdef DEBUG_THREADING
  write(6,'(a,i3,a,i3)') " [Driver_initParallel]: MPI rank ", dr_globalMe, &
       " has a team that includes thread ", omp_get_thread_num()
# endif
  !$omp end parallel
#endif

  dr_mpiThreadSupport = mpiThreadSupport

end subroutine Driver_initParallel
