!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      include 'mpif.h'

      integer ierror,mype,nprocs

   
      call mpi_init(ierror)

      call mpi_comm_rank(amr_mpi_meshComm, mype)
      call mpi_comm_size(amr_mpi_meshComm, nprocs)

      write(*,*) 'proc ',mype,' nprocs ',nprocs

      call mpi_barrier(amr_mpi_meshComm, ierror)

      write(*,*) 'proc ',mype,' nprocs ',nprocs
      call mpi_finalize(ierror)

      end
