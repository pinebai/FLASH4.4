!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_checkpoint_wr
!! NAME
!!
!!   amr_checkpoint_wr
!!
!! SYNOPSIS
!!
!!   call amr_checkpoint_wr(file_num)
!!   call amr_checkpoint_wr(file_num, l_with_guardcells, check_format, 
!!                          user_attr1, user_attr2, user_attr3, 
!!                          user_attr4, user_attr5)
!!
!!   call amr_checkpoint_wr(integer, optional logical, optional char*80,
!!                          optional real, optional real, optional real, 
!!                          optional real, optional real)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: file_num
!!     An integer number which be appended to the end of the file name.
!!
!!   optional, logical, intent(in) :: l_with_guardcells
!!     If true, then guardcells are included in the checkpoint file.  Otherwise 
!!     (the default) they are not included.
!!   
!!   optional, character(len=80), intent(in) :: check_format
!!     Argument describing what type of output to use.  Currently supports
!!     fortran binary(the default), parallel hdf5 output, or parallel mpiio
!!     in native binary.
!!     To produce parallel hdf5 output:
!!       character(len=80) :: checkf
!!       checkf = 'hdf5'
!!       call amr_checkpoint_wr(...., check_format=checkf, ....) 
!!
!!   optional, real, intent(in) :: user_attr1(2,3,4,5)
!!     Arguments which allow the user to add some extra information to the file.  
!!     Currently only 5 real numbers can be added.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_checkpoint_wr_default
!!   amr_checkpoint_wr_hdf5
!!   amr_checkpoint_wr_mpiio
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit a checkpoint file has been written.
!!
!! DESCRIPTION
!! 
!!  Subroutine to checkpoint runs using PARAMESH.
!!  It writes out the tree data structure and data stored in PARAMESH blocks.
!!  Optionally, a user may add a small amout of attribute data to the files written.
!!  If the default bevaviour is selected, writes are done serially by processor 0. 
!!  I.e. data is collected from other processors and then written out.  
!!  The default behaviour USES UNFORMATTED DIRECT I/O.
!!  Optionally, the user can specify if the file will be written using the HDF5
!!  portable data format.  If supported on the system, selecting hdf5 output will 
!!  also result in the file will be written in parallel.
!!
!!  The files produced will have names of the form 'paramesh_chk_######' or 
!!  'paramesh_chk_######.hdf5'. where '######' is the file_num argument passed into
!!  this routine.
!!
!! AUTHORS
!!
!!   Peter MacNeice (1997) with modifications by Kevin Olson for parallel hdf5 and 
!!   mpiio (2004-2005).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      subroutine amr_checkpoint_wr(file_num, & 
     &                             l_with_guardcells, & 
     &                             check_format, & 
     &                             user_attr_1, & 
     &                             user_attr_2, & 
     &                             user_attr_3, & 
     &                             user_attr_4, & 
     &                             user_attr_5)

      use paramesh_interfaces, only : amr_checkpoint_wr_default, & 
     &                                amr_checkpoint_wr_hdf5, & 
     &                                amr_checkpoint_wr_mpiio
      use paramesh_comm_data

      implicit none

      include 'mpif.h'

      integer, intent(in)                      :: file_num
      logical, optional, intent(in)            :: l_with_guardcells
      character (len=80), optional, intent(in) :: check_format
      real, optional, intent(in)               :: user_attr_1, & 
     &                                            user_attr_2, & 
     &                                            user_attr_3, & 
     &                                            user_attr_4, & 
     &                                            user_attr_5

      character (len=80) :: check_format_in
      integer            :: mype, ierr
      
      call MPI_COMM_RANK (amr_mpi_meshComm, mype, ierr)

      check_format_in = 'default'
      if (present(check_format)) check_format_in = check_format
      check_format_in = trim(check_format_in)

      if (check_format_in(1:7) == 'default') then
         call amr_checkpoint_wr_default (file_num, l_with_guardcells, & 
     &     user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5)
      elseif (check_format_in(1:4) == 'hdf5') then
         call amr_checkpoint_wr_hdf5 (file_num, l_with_guardcells, & 
     &     user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5)
      elseif (check_format_in(1:5) == 'mpiio') then
         call amr_checkpoint_wr_mpiio (file_num, l_with_guardcells, & 
     &     user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5)
      else
         if (mype == 0) then
            print *,' UNRECOGNIZED I/O FORMAT, CHECKPOINT NOT WRITTEN '
         end if
      end if

      end subroutine amr_checkpoint_wr



!-----------------------------------------------------------------------
      

!!****f* mpi_source/amr_checkpoint_re
!! NAME
!!
!!   amr_checkpoint_re
!!
!! SYNOPSIS
!!
!!   call amr_checkpoint_re(file_num)
!!   call amr_checkpoint_re(file_num, l_with_guardcells, check_format, 
!!                          user_attr1, user_attr2, user_attr3, 
!!                          user_attr4, user_attr5)
!!
!!   call amr_checkpoint_re(integer, optional logical, optional char*80,
!!                          optional real, optional real, optional real, 
!!                          optional real, optional real)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: file_num
!!     An integer number which be appended to the end of the file name.
!!
!!   optional, logical, intent(in) :: l_with_guardcells
!!     If true, then guardcells are included in the checkpoint file.  Otherwise 
!!     (the default) they are not included.
!!   
!!   optional, character(len=80), intent(in) :: check_format
!!     Argument describing what type of output has been used.  Currently only supports
!!     fortran binary(the default), parallel hdf5 output, or parallel mpiio in
!!     in native binary.
!!     To read and hdf5 checkpoint file:
!!       character(len=80) :: checkf
!!       checkf = 'hdf5'
!!       call amr_checkpoint_re(...., check_format=checkf, ....) 
!!
!!   optional, real, intent(out) :: user_attr1(2,3,4,5)
!!     Arguments which allow the user to add some extra information to the file.  
!!     Currently only 5 real numbers can be added. Note also that the intent for these
!!     arguments is different than for the amr_checkpoint_wr routine.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_checkpoint_re_default
!!   amr_checkpoint_re_hdf5
!!   amr_checkpoint_re_mpiio
!!   amr_abort
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit a checkpoint file has been read in.
!!
!! DESCRIPTION
!! 
!!  Subroutine to read checkpoint files written by amr_checkpoint_wr routine.
!!  It read in the tree data structure and data stored in the PARAMESH blocks.
!!  Optionally, a user may read a small amout of attribute data fime the file.
!!  If the default bevaviour is selected, reads are done serially by processor 0. 
!!  The data is sent to the other processors from processor 0.  
!!  The default behaviour USES UNFORMATTED DIRECT I/O.
!!  Optionally, the user can specify if the file will be read in using the HDF5
!!  portable data format.  If supported on the system, selecting hdf5 will 
!!  also result in the file will be read in parallel.
!!
!!  The files read in must have names of the form 'paramesh_chk_######' or 
!!  'paramesh_chk_######.hdf5'. where '######' is the file_num argument passed into
!!  this routine.
!!
!! AUTHORS
!!
!!   Peter MacNeice (1997) with modifications by Kevin Olson for parallel hdf5 and
!!   mpiio (2004-2005).
!!
!!***

      subroutine amr_checkpoint_re (file_num, & 
     &                              l_with_guardcells, & 
     &                              check_format, & 
     &                              user_attr_1, & 
     &                              user_attr_2, & 
     &                              user_attr_3, & 
     &                              user_attr_4, & 
     &                              user_attr_5)


      Use paramesh_comm_data

      use paramesh_interfaces, only : amr_checkpoint_re_default, & 
     &                                amr_checkpoint_re_hdf5, & 
     &                                amr_checkpoint_re_mpiio, & 
     &                                amr_abort

      implicit none

      include 'mpif.h'

      integer, intent(in)                      :: file_num
      logical, optional, intent(in)            :: l_with_guardcells
      character (len=80), optional, intent(in) :: check_format
      real, optional, intent(out)              :: user_attr_1, & 
     &                                            user_attr_2, & 
     &                                            user_attr_3, & 
     &                                            user_attr_4, & 
     &                                            user_attr_5

      character (len=80) :: check_format_in
      integer            :: mype, ierr
      
      call MPI_COMM_RANK (amr_mpi_meshComm, mype, ierr)

      check_format_in = 'default'
      if (present(check_format)) check_format_in = check_format
      check_format_in = trim(check_format_in)

      if (check_format_in(1:7) == 'default') then
         call amr_checkpoint_re_default (file_num, l_with_guardcells, & 
     &     user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5)
      elseif (check_format_in(1:4) == 'hdf5') then
         call amr_checkpoint_re_hdf5    (file_num, l_with_guardcells, & 
     &     user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5)
      elseif (check_format_in(1:5) == 'mpiio') then
         call amr_checkpoint_re_mpiio    (file_num, l_with_guardcells, & 
     &     user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5)
      else
         if (mype == 0) then
            print *,' UNRECOGNIZED I/O FORMAT,  & 
     & UNABLE TO READ CHECKPOINT FILE '
         end if
         call amr_abort()
      end if

      ! broadcast user_attr's to all processors
      if (present(user_attr_1)) then
         call MPI_BCAST(user_attr_1,1,amr_mpi_real,0, & 
     &                  amr_mpi_meshComm, & 
     &                  ierr)
      endif
      if (present(user_attr_2)) then
         call MPI_BCAST(user_attr_2,1,amr_mpi_real,0, & 
     &                  amr_mpi_meshComm, & 
     &                  ierr)
      endif
      if (present(user_attr_3)) then
         call MPI_BCAST(user_attr_3,1,amr_mpi_real,0, & 
     &                  amr_mpi_meshComm, & 
     &                  ierr)
      endif
      if (present(user_attr_4)) then
         call MPI_BCAST(user_attr_4,1,amr_mpi_real,0, & 
     &                  amr_mpi_meshComm, & 
     &                  ierr)
      endif
      if (present(user_attr_5)) then
         call MPI_BCAST(user_attr_5,1,amr_mpi_real,0, & 
     &                  amr_mpi_meshComm, & 
     &                  ierr)
      endif

      end subroutine amr_checkpoint_re


