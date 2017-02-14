!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      subroutine amr_plotfile_chombo (file_num)
      use paramesh_comm_data

      implicit none
      include 'mpif.h'

      integer, intent(in) :: file_num
      integer :: mype, ierr

      call MPI_COMM_RANK (amr_mpi_meshComm,mype,ierr)

      if (mype == 0) then
         print *,' WARNING: you are calling amr_plotfile_chombo '
         print *,'          but your version of paramesh is not '
         print *,'          yet configured to do this.          '
         print *,'          Go to utilities/io/plotting/chombovis '
         print *,'          in the main paramesh directory, run  '
         print *,'          the INSTALL script, and recompile !!!'
      end if

      return
      end subroutine amr_plotfile_chombo

!----------------------------------------------------------------------

