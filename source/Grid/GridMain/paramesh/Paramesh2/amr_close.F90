      subroutine amr_close


! $RCSfile: amr_close.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $

      implicit none

      include 'mpif.h'
      
      integer ierr

! This subroutine closes the amr package. 

!----------------------------------------------------------------

! Call the machine/software environment specific closure routine.
! Different versions of comm_finish are provided for use with machines
! which run shmem or mpi. Make sure to compile with the appropriate version
! for your environment.
      call MPI_FINALIZE(ierr)

      return
      end
